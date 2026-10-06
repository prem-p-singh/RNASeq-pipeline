#!/usr/bin/env python3
"""Durable SSH bridge to submit.sh. Requires POSIX file locks and Python 3.9+."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[2]


def save(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def lock(path):
    handle = path.open('a')
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        raise RuntimeError('Another controller is running for this project')
    return handle


def status(run):
    state = json.loads((run / 'status.json').read_text())
    if state['status'] == 'running':
        try:
            handle = lock(run / 'run.lock')
        except RuntimeError:
            pass
        else:
            handle.close()
            state.update(status='interrupted', error='Controller stopped without a final status. Inspect jobs and logs before resuming.')
    return state


def source_identity(repo):
    """Hash executable/configuration inputs, independent of Git or archive layout."""
    repo = Path(repo).resolve()
    if not (repo / 'submit.sh').is_file():
        raise ValueError('pipeline_repo must contain submit.sh')
    digest = hashlib.sha256()
    paths = [repo / 'submit.sh', repo / 'Snakefile']
    for folder in ('scripts', 'workflow', 'config', 'environments', 'profiles'):
        paths.extend(p for p in (repo / folder).rglob('*')
                     if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc')
    for path in sorted(paths):
        if path.is_file():
            digest.update(str(path.relative_to(repo)).encode() + b'\0')
            digest.update(hashlib.sha256(path.read_bytes()).digest())
    launcher = (repo / 'submit.sh').read_text()
    questions = repo / 'config/intake_questions.yaml'
    schema = re.search(r'^schema_version:\s*(\d+)', questions.read_text(), re.M) if questions.is_file() else None
    return {'intake_schema': int(schema[1]) if schema else None, 'pipeline_repo': str(repo), 'source_sha256': digest.hexdigest(),
            'workbook_revisions': '--accept-revision' in launcher,
            'tagseq': (repo / 'scripts/tagseq.py').is_file(),
            'small_rna': (repo / 'workflow/rules/small_rna.smk').is_file(),
            'modules': [name for name in ('star', 'screen', 'srna', 'umi')
                        if (repo / f'environments/{name}-linux-64.explicit.txt').is_file()]}


def request(path):
    data = json.loads(path.read_text())
    allowed = {'project', 'conda_bin', 'executor', 'cores', 'profile', 'mode', 'intake',
               'intake_sheet', 'env_prefix', 'star_env_prefix', 'screen_env_prefix',
               'srna_env_prefix', 'umi_env_prefix', 'pipeline_repo', 'expected_source_sha256',
               'accept_revision'}
    if set(data) - allowed:
        raise ValueError('Unknown request fields: ' + ', '.join(sorted(set(data) - allowed)))
    for key in ('project', 'conda_bin', 'pipeline_repo'):
        if not isinstance(data.get(key), str) or not Path(data[key]).is_absolute():
            raise ValueError(key + ' must be an absolute path')
    project = Path(data['project']).resolve()
    repo = Path(data['pipeline_repo']).resolve()
    if project == repo or repo in project.parents or project in repo.parents:
        raise ValueError('Project must be separate from the source checkout')
    if not (Path(data['conda_bin']) / 'conda').is_file():
        raise ValueError('conda_bin must contain the conda executable')
    if data.get('executor', 'local') not in ('local', 'slurm'):
        raise ValueError('executor must be local or slurm')
    if data.get('mode', 'dry_run') not in ('dry_run', 'run', 'raw_qc'):
        raise ValueError('mode must be dry_run, run or raw_qc')
    cores = data.get('cores', 4)
    if type(cores) is not int or not 1 <= cores <= 4096:
        raise ValueError('cores must be an integer from 1 to 4096')
    if data.get('profile') not in (None, 'small', 'medium', 'large'):
        raise ValueError('profile must be small, medium or large')
    for key in ('intake', 'env_prefix', 'star_env_prefix', 'screen_env_prefix', 'srna_env_prefix', 'umi_env_prefix'):
        if key in data and (not isinstance(data[key], str) or not Path(data[key]).is_absolute()):
            raise ValueError(key + ' must be an absolute path')
    if data.get('intake_sheet') and not data.get('intake'):
        raise ValueError('intake_sheet requires intake')
    if data.get('intake'):
        if not Path(data['intake']).is_file():
            raise ValueError('Intake file does not exist')
    elif not (project / 'config/config.yaml').is_file():
        raise ValueError('Provide an intake or an existing project config/config.yaml')
    if type(data.get('accept_revision', False)) is not bool:
        raise ValueError('accept_revision must be true or false')
    if data.get('accept_revision') and not data.get('intake'):
        raise ValueError('accept_revision requires intake')
    if data.get('intake_sheet') and not isinstance(data['intake_sheet'], str):
        raise ValueError('intake_sheet must be text')
    identity = source_identity(repo)
    expected = data.get('expected_source_sha256')
    if not isinstance(expected, str) or not re.fullmatch('[0-9a-f]{64}', expected):
        raise ValueError('expected_source_sha256 must pin the reviewed pipeline source; use inspect')
    if identity['source_sha256'] != expected:
        raise ValueError('Pipeline source differs from expected_source_sha256; inspect the checkout before updating the request')
    if data.get('accept_revision') and not identity['workbook_revisions']:
        raise ValueError('Selected pipeline does not support workbook revisions')
    if data.get('mode') == 'raw_qc':
        if data.get('intake') or data.get('accept_revision'):
            raise ValueError('raw_qc needs an existing canonical project; omit intake and accept_revision')
        if data.get('executor', 'local') != 'local':
            raise ValueError('Standalone raw_qc runs locally; use a worker allocation or dedicated host')
        if cores > 16:
            raise ValueError('raw_qc supports at most 16 cores')
    data['project'] = str(project)
    data['pipeline_repo'] = str(repo)
    return data


def command(data):
    repo = Path(data['pipeline_repo'])
    cmd = ['bash', str(repo / 'submit.sh'), '-d', data['project'],
           '--executor', data.get('executor', 'local'), '--cores', str(data.get('cores', 4))]
    for field, flag in [('profile', '--profile'), ('intake', '--intake'), ('intake_sheet', '--intake-sheet')]:
        if data.get(field):
            cmd.extend([flag, data[field]])
    if data.get('accept_revision'):
        cmd.append('--accept-revision')
    if data.get('mode', 'dry_run') == 'dry_run':
        cmd.append('--dry-run')
    cmd.extend(['--', '--rerun-incomplete'])
    return cmd


def work(run, descriptors):
    # The project and run locks stay inherited by the launcher until it exits.
    data = json.loads((run / 'request.json').read_text())
    state = json.loads((run / 'status.json').read_text())
    env = os.environ.copy()
    # Re-activating the same Conda prefix may leave PATH unchanged. Keep its
    # interpreter ahead of the base Conda bin directory on allocated workers.
    active_bin = str(Path(env['CONDA_PREFIX']) / 'bin') if env.get('CONDA_PREFIX') else ''
    env['PATH'] = os.pathsep.join(p for p in (active_bin, data['conda_bin'], env.get('PATH', '')) if p)
    for key, name in [('env_prefix', 'RNASEQ_ENV_PREFIX'), ('star_env_prefix', 'RNASEQ_STAR_ENV_PREFIX'),
                      ('screen_env_prefix', 'RNASEQ_SCREEN_ENV_PREFIX'),
                      ('srna_env_prefix', 'RNASEQ_SRNA_ENV_PREFIX'), ('umi_env_prefix', 'RNASEQ_UMI_ENV_PREFIX')]:
        if key in data:
            env[name] = data[key]
    try:
        # Recheck immediately before execution, including after a delayed child start.
        if source_identity(data['pipeline_repo'])['source_sha256'] != data['expected_source_sha256']:
            raise ValueError('Pipeline source changed after launch')
        cmd = command(data)
        if data.get('mode') == 'raw_qc':
            # Positional parameters preserve paths; no request text is evaluated as shell code.
            script = ('set -euo pipefail; prefix=$(bash "$1/scripts/bootstrap.sh"); '
                      'source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$prefix"; '
                      'exec python3 "$1/scripts/raw_qc.py" -d "$2" --out "$3" --threads "$4"')
            cmd = ['bash', '-c', script, 'raw-qc', data['pipeline_repo'], data['project'],
                   str(run / 'raw_qc'), str(data.get('cores', 4))]
        state['command'] = cmd
        save(run / 'status.json', state)
        result = subprocess.run(cmd, env=env, stdin=subprocess.DEVNULL, pass_fds=descriptors)
        state.update(exit_code=result.returncode, status='succeeded' if result.returncode == 0 else 'failed')
    except Exception as exc:
        state.update(status='failed', error=str(exc))
    # Only expose current-run gates; never relabel older gates as this run's results.
    state['reports'] = {}
    for name in ('recommendation', 'preflight_plan', 'resource_plan'):
        path = Path(data['project']) / 'gates' / (name + '.json')
        if data.get('mode') != 'raw_qc' and path.is_file() and path.stat().st_mtime >= state['started_at']:
            try:
                state['reports'][name] = json.loads(path.read_text())
            except (ValueError, OSError) as exc:
                state.setdefault('report_warnings', []).append(str(exc))
    if data.get('mode') == 'raw_qc':
        state['raw_qc_report'] = str(run / 'raw_qc/index.html')
    state['finished_at'] = time.time()
    save(run / 'status.json', state)


def start(data, run_id):
    root = Path(data['project']) / '.n8n'
    root.mkdir(parents=True, exist_ok=True)
    run = root / 'runs' / run_id
    if (run / 'status.json').exists():
        if json.loads((run / 'request.json').read_text()) != data:
            raise ValueError('Run ID already belongs to a different request; use a new run_id')
        return status(run)
    with lock(root / 'controller.lock') as project_lock:
        # Recheck after locking in case another start won the race.
        if run.exists():
            raise ValueError('Run directory exists; inspect it and choose a new run_id')
        run.mkdir(parents=True)
        with lock(run / 'run.lock') as run_lock:
            save(run / 'request.json', data)
            state = {'run_id': run_id, 'status': 'running', 'mode': data.get('mode', 'dry_run'),
                     'project': data['project'], 'started_at': time.time(),
                     'pipeline': source_identity(data['pipeline_repo']),
                     'log': str(run / 'controller.log'), 'command': command(data),
                     'review': ['gates/preflight_issues.tsv', 'gates/recommendation.json',
                                'gates/resource_plan.json', 'config/config.yaml'],
                     'result_note': 'Inspect project.output_dir and optional-analysis manifests. A dry run produces no scientific results.'}
            save(run / 'status.json', state)
            descriptors = (project_lock.fileno(), run_lock.fileno())
            try:
                with (run / 'controller.log').open('ab') as output:
                    subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '_work',
                                      '--run', str(run), '--fds', ','.join(map(str, descriptors))],
                                     stdin=subprocess.DEVNULL, stdout=output, stderr=subprocess.STDOUT,
                                     start_new_session=True, pass_fds=descriptors)
            except Exception as exc:
                state.update(status='failed', error=str(exc))
                save(run / 'status.json', state)
                raise
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['start', 'status', 'inspect', '_work'])
    parser.add_argument('--request', type=Path)
    parser.add_argument('--project', type=Path)
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--run-id')
    parser.add_argument('--run', type=Path)
    parser.add_argument('--fds')
    args = parser.parse_args()
    if args.action == 'inspect':
        if args.repo is None or not args.repo.is_absolute():
            parser.error('inspect requires an absolute --repo')
        print(json.dumps(source_identity(args.repo)))
        return
    if args.action == '_work':
        work(args.run, tuple(map(int, args.fds.split(','))))
        return
    if not args.run_id or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}', args.run_id):
        parser.error('run-id must be 1–80 letters, digits, underscores or hyphens')
    if args.action == 'start':
        if args.request is None:
            parser.error('start requires --request')
        data = request(args.request)
        if args.project is None or args.project.resolve() != Path(data['project']):
            parser.error('--project must match the request project')
        answer = start(data, args.run_id)
    else:
        if args.project is None or not args.project.is_absolute():
            parser.error('status requires an absolute --project')
        answer = status(args.project.resolve() / '.n8n/runs' / args.run_id)
    print(json.dumps(answer))


if __name__ == '__main__':
    main()
