#!/usr/bin/env python3
"""Linux integration smoke: real release launcher/DAG through the n8n bridge.

All assay routes get dry runs. Animal small RNA and QuantSeq UMI also execute
on synthetic reads. This checks orchestration, not new scientific qualification.
Run within a worker allocation, with RNASEQ_*_ENV_PREFIX set for locked modules.
"""
import argparse
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import platform
import tempfile
import time


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--gate', action='store_true', help='Check deployment prerequisites; exit 77 when unavailable')
    parser.add_argument('--raw-qc-only', action='store_true')
    parser.add_argument('--only', choices=['FWD','REV','FWD-UMI','animal','plant','bulk-salmon','bulk-star'])
    args = parser.parse_args()
    if args.gate:
        missing = []
        if platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'AMD64'):
            missing.append('Linux x86_64 worker allocation')
        if importlib.util.find_spec('yaml') is None:
            missing.append('PyYAML')
        source = os.environ.get('RNASEQ_N8N_SOURCE')
        if args.repo is None and source:
            args.repo = Path(source)
        if args.repo is None or not (args.repo / 'submit.sh').is_file():
            missing.append('RNASEQ_N8N_SOURCE pointing to the reviewed schema-7 bundle')
        for name in ('CONDA_EXE', 'RNASEQ_ENV_PREFIX', 'RNASEQ_STAR_ENV_PREFIX',
                     'RNASEQ_SCREEN_ENV_PREFIX', 'RNASEQ_SRNA_ENV_PREFIX', 'RNASEQ_UMI_ENV_PREFIX'):
            value = os.environ.get(name)
            if not value or not (Path(value).is_file() if name == 'CONDA_EXE'
                                 else (Path(value) / 'conda-meta/history').is_file()):
                missing.append(name + ' (installed locked runtime)')
        if missing:
            print('SKIP: ' + '; '.join(missing))
            return 77
        # An explicitly supplied but wrong source is a failure, never a skip.
        local = Path(__file__).resolve().parents[1] / 'integrations/n8n'
        reviewed = json.loads((local / 'source.review.json').read_text())
        identity = load(local / 'runner.py', 'gate_bridge').source_identity(args.repo)
        assert identity['intake_schema'] == 7, identity
        assert identity['source_sha256'] == reviewed['source_sha256'], identity
        assert (args.repo / 'integrations/n8n/runner.py').read_bytes() == (local / 'runner.py').read_bytes(), 'Rebuild bundle: bridge differs from gate checkout'
        if args.out is None:
            args.out = Path(tempfile.mkdtemp(prefix='rnaseq-n8n-release-')) / 'results'
    if args.repo is None or args.out is None:
        parser.error('--repo and --out are required unless --gate supplies them')
    import yaml
    repo, out = args.repo.resolve(), args.out.resolve()
    print('Release smoke output: ' + str(out), flush=True)
    out.mkdir(parents=True, exist_ok=False)
    bridge = load(repo / 'integrations/n8n/runner.py', 'bridge')
    identity = bridge.source_identity(repo)
    assert identity['tagseq'] and identity['small_rna'] and identity['workbook_revisions']
    if args.raw_qc_only:
        sys.path.insert(0, str(repo / 'scripts'))
        from metadata import write_tables
        project = out / 'raw-qc-project'
        (project / 'config').mkdir(parents=True)
        (project / 'config/config.yaml').write_text('samples:\n  metadata_dir: metadata\n  sheet: metadata/samples.tsv\n  seq_type: rnaseq_single\n')
        read = project / 'sample.fastq.gz'
        with gzip.open(read, 'wt') as f:
            f.write('@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n')
        original = read.read_bytes()
        write_tables({'samples':[{'sample_id':'001'}],
            'libraries':[{'sample_id':'001','library_id':'L1','layout':'single','assay_family':'bulk','umi':'no','strandedness':'unknown'}],
            'reads':[{'read_unit_id':'R1','library_id':'L1','role':'R1','uri':str(read)}]},project/'metadata')
        data={'pipeline_repo':str(repo),'expected_source_sha256':identity['source_sha256'],
              'project':str(project),'conda_bin':str(Path(os.environ['CONDA_EXE']).parent),
              'mode':'raw_qc','executor':'local','cores':2,'env_prefix':os.environ['RNASEQ_ENV_PREFIX']}
        req=out/'request.json';req.write_text(json.dumps(data))
        state=bridge.start(bridge.request(req),'raw-qc')
        deadline=time.monotonic()+300
        while state['status']=='running' and time.monotonic()<deadline:
            time.sleep(1);state=bridge.status(project/'.n8n/runs/raw-qc')
        assert state['status']=='succeeded',state
        assert Path(state['raw_qc_report']).stat().st_size>0
        assert read.read_bytes()==original
        (out/'summary.json').write_text(json.dumps(state,indent=2)+'\n')
        print('PASS: independent raw QC through bridge, real fastp report and unchanged source reads',flush=True)
        return
    srna = load(repo / 'tests/fixtures/small_rna_workflow_project.py', 'srna_fixture')
    bulk = load(repo / 'tests/fixtures/fastq_project.py', 'bulk_fixture')
    # Reuse the release's own trusted synthetic kit builder, without running
    # the surrounding shell test or duplicating kit geometry in this test.
    tag_test = (repo / 'tests/check_tagseq_workflow.sh').read_text()
    tag_builder = tag_test.split("<<'PY'\n", 1)[1].split('\nPY\n}', 1)[0]
    assert 'expected.json' in tag_builder
    projects = []
    for kit in ('FWD', 'REV', 'FWD-UMI'):
        project = out / kit
        subprocess.run([sys.executable, '-c', tag_builder, str(project), 'QuantSeq '+kit, str(repo)], check=True)
        projects.append((project, 'star_tagseq'))
    for profile in ('animal', 'plant'):
        project = out / profile
        srna.build(project, repo)
        if profile == 'plant':
            cfg_path = project / 'config/config.yaml'
            cfg = yaml.safe_load(cfg_path.read_text())
            cfg['analysis']['small_rna_profile'] = 'plant'
            cfg['reference'].pop('mature_mirna_url')
            cfg['organism'].pop('mirbase_prefix')
            genome = project / 'inputs/genome.fa'
            genome.write_text('>chr1\n' + 'ACGT'*10000 + '\n')
            cfg['reference']['genome_fasta_url'] = genome.as_uri()
            cfg_path.write_text(yaml.safe_dump(cfg))
        projects.append((project, 'small_rna_'+profile))
    project = out / 'bulk-salmon'
    bulk.build(project, repo, 'single')
    projects.append((project, 'salmon_bulk'))
    project = out / 'bulk-star'
    subprocess.run([sys.executable, '-c', tag_builder, str(project), 'QuantSeq FWD', str(repo)], check=True)
    cfg_path = project / 'config/config.yaml'
    cfg = yaml.safe_load(cfg_path.read_text())
    cfg['samples']['seq_type'] = 'rnaseq_single'
    cfg['analysis'].pop('kit')
    cfg['analysis']['quantifier'] = 'star'
    cfg['samples']['expected_libtype'] = 'SF'
    cfg_path.write_text(yaml.safe_dump(cfg))
    projects.append((project, 'star_counts'))
    if args.only:
        projects = [(p,r) for p,r in projects if p.name == args.only]
    results=[]
    for project, route in projects:
        modes = ['dry_run', 'run'] if project.name in ('animal', 'FWD-UMI') else ['dry_run']
        for mode in modes:
            request = {'pipeline_repo': str(repo), 'expected_source_sha256':identity['source_sha256'],
                       'project':str(project), 'conda_bin':str(Path(os.environ['CONDA_EXE']).parent),
                       'mode':mode, 'executor':'local', 'cores':4}
            for module, key in [('','env_prefix'),('STAR_','star_env_prefix'),('SCREEN_','screen_env_prefix'),
                                ('SRNA_','srna_env_prefix'),('UMI_','umi_env_prefix')]:
                request[key]=os.environ['RNASEQ_'+module+'ENV_PREFIX']
            request_path = project / (mode+'.json')
            request_path.write_text(json.dumps(request))
            result = subprocess.run([sys.executable,str(repo/'integrations/n8n/runner.py'),'start',
                '--request',str(request_path),'--project',str(project),'--run-id',mode],capture_output=True,text=True,check=True)
            state=json.loads(result.stdout)
            deadline=time.monotonic()+1200
            while state['status']=='running' and time.monotonic()<deadline:
                time.sleep(1)
                state=bridge.status(project/'.n8n/runs'/mode)
            if state['status']!='succeeded':
                print(Path(state['log']).read_text()[-10000:],flush=True)
                raise AssertionError(state)
            assert state['reports']['recommendation']['route']==route,state
            results.append({'project':project.name,'mode':mode,'route':route,'status':state['status']})
            (out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
            print('PASS', project.name, mode, route, flush=True)
    # Check actual outputs, not only exit codes.
    for project in ('animal','FWD-UMI') if not args.only else [args.only] if args.only in ('animal','FWD-UMI') else []:
        assert (out/project/'results/counts.tsv').stat().st_size>0
    if not args.only or args.only == 'animal':
        with (out/'animal/results/counts.tsv').open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        assert {r['gene_id'] for r in rows}=={'fix-miR-1','fix-miR-3'}
    print('PASS: requested assay/kit DAGs and selected executions through the bridge', flush=True)


if __name__=='__main__':
    sys.exit(main())
