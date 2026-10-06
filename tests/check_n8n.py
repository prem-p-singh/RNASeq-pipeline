#!/usr/bin/env python3
"""Exercise the SSH bridge with a fake launcher; no scientific jobs are submitted."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

SOURCE = Path(__file__).resolve().parents[1]


def main():
    # This fixture controls PATH to test the bridge's Conda precedence. A host
    # BASH_ENV can replace PATH again on shell startup (including agent hooks).
    # Isolate only this test process; retain the exact interpreter assertion.
    os.environ.pop('BASH_ENV', None)
    workflow = json.loads((SOURCE / 'integrations/n8n/rnaseq.workflow.json').read_text())
    names = {node['name'] for node in workflow['nodes']}
    assert len(names) == len(workflow['nodes']) and not workflow['active']
    for source, links in workflow['connections'].items():
        assert source in names
        for branch in links['main']:
            assert all(link['node'] in names for link in branch)
    with tempfile.TemporaryDirectory(prefix='rnaseq n8n ') as tmp:
        root = Path(tmp)
        repo = root / 'repo'
        runner = repo / 'integrations/n8n/runner.py'
        runner.parent.mkdir(parents=True)
        shutil.copy(SOURCE / 'integrations/n8n/runner.py', runner)
        (repo / 'submit.sh').write_text('#!/bin/bash\nsleep 0.5\ncommand -v python3\nprintf "%s\\n" "$@"\nif [ -f "' + str(root / 'fail') + '" ]; then exit 7; fi\n')
        spec = importlib.util.spec_from_file_location('bridge', runner)
        bridge = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(bridge)
        project = root / "study ' with spaces"
        (project / 'config').mkdir(parents=True)
        (project / 'config/config.yaml').write_text('project: {}\n')
        conda = root / 'conda/bin'
        conda.mkdir(parents=True)
        (conda / 'conda').touch()
        data = {'project': str(project), 'conda_bin': str(conda), 'mode': 'dry_run', 'pipeline_repo': str(repo),
                'expected_source_sha256': bridge.source_identity(repo)['source_sha256']}
        request = root / 'request.json'
        request.write_text(json.dumps(data))
        data = bridge.request(request)

        def finish(run_id):
            run = project / '.n8n/runs' / run_id
            for _ in range(100):
                result = bridge.status(run)
                if result['status'] != 'running':
                    # Terminal status is published just before the worker exits.
                    try:
                        released = bridge.lock(project / '.n8n/controller.lock')
                    except RuntimeError:
                        time.sleep(0.05)
                        continue
                    released.close()
                    return result
                time.sleep(0.05)
            raise AssertionError('Worker did not finish')

        bad_project = subprocess.run([sys.executable, str(runner), 'start', '--request', str(request),
                                      '--project', str(root / 'wrong'), '--run-id', 'wrong'], capture_output=True)
        assert bad_project.returncode != 0 and not (project / '.n8n').exists()
        book = root / 'intake.xlsx'; book.touch()
        assert '--accept-revision' in bridge.command({**data, 'intake': str(book), 'accept_revision': True})
        active = root / 'active/bin'; active.mkdir(parents=True)
        (active / 'python3').symlink_to(sys.executable)
        old_prefix = os.environ.get('CONDA_PREFIX')
        os.environ['CONDA_PREFIX'] = str(active.parent)
        try:
            first = bridge.start(data, 'check1')
        finally:
            if old_prefix is None:
                os.environ.pop('CONDA_PREFIX')
            else:
                os.environ['CONDA_PREFIX'] = old_prefix
        assert first['status'] == 'running'
        assert bridge.start(data, 'check1')['run_id'] == 'check1'
        try:
            bridge.start(data, 'conflict')
        except RuntimeError:
            pass
        else:
            raise AssertionError('Concurrent controller accepted')
        done = finish('check1')
        assert done['status'] == 'succeeded' and done['exit_code'] == 0
        log = Path(done['log']).read_text()
        assert log.splitlines()[0] == str(active / 'python3'), log
        assert '--dry-run\n' in log and str(project) + '\n' in log and '--rerun-incomplete' in log
        assert bridge.start(data, 'check1')['status'] == 'succeeded'
        try:
            bridge.start({**data, 'mode': 'run'}, 'check1')
        except ValueError:
            pass
        else:
            raise AssertionError('Reused run ID accepted changed request')
        (root / 'fail').touch()
        bridge.start({**data, 'mode': 'run'}, 'failure')
        assert finish('failure')['exit_code'] == 7
        run = project / '.n8n/runs/check1'
        state = json.loads((run / 'status.json').read_text())
        state['status'] = 'running'
        bridge.save(run / 'status.json', state)
        assert bridge.status(run)['status'] == 'interrupted'
        invalid = subprocess.run([sys.executable, str(runner), 'status', '--project', str(project),
                                  '--run-id', '../escape'], capture_output=True)
        assert invalid.returncode != 0
        for change in ({'expected_source_sha256': '0'*64}, {'accept_revision': 'yes'}, {'accept_revision': True}, {'cores': True}, {'executor': 'anything'}, {'extra': True}, {'project': str(repo)}):
            request.write_text(json.dumps({**data, **change}))
            try:
                bridge.request(request)
            except ValueError:
                pass
            else:
                raise AssertionError('Invalid request accepted: ' + repr(change))
    print('PASS: graph, durable completion/failure, idempotency, concurrency, interruption, quoting and validation')


if __name__ == '__main__':
    main()
