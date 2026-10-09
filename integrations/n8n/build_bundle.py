#!/usr/bin/env python3
"""Package the reviewed release source and its explicit integration patch."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True, help='Git checkout containing the reviewed commit')
    p.add_argument('--out', type=Path, required=True)
    args=p.parse_args()
    review=json.loads((HERE/'source.review.json').read_text())
    raw=subprocess.check_output(['git','-C',str(args.source),'archive',review['reviewed_commit']])
    with tempfile.TemporaryDirectory(prefix='rnaseq-n8n-bundle-') as temporary:
        folder=Path(temporary)
        # This archive is generated from the explicit local source commit above.
        with tarfile.open(fileobj=io.BytesIO(raw)) as archive:
            archive.extractall(folder)
        for patch in review['integration_patches']:
            subprocess.run(['git','apply',str(HERE/patch)],cwd=folder,check=True)
        shutil.copytree(HERE,folder/'integrations/n8n',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for name in ('check_n8n.py','check_n8n_nodes.js','check_n8n_release.py'):
            shutil.copy2(ROOT/'tests'/name,folder/'tests'/name)
        spec=importlib.util.spec_from_file_location('bridge',HERE/'runner.py')
        bridge=importlib.util.module_from_spec(spec);spec.loader.exec_module(bridge)
        actual=bridge.source_identity(folder)['source_sha256']
        if actual!=review['source_sha256']:
            raise ValueError('Bundle source hash does not match the reviewed integration')
        readme=folder/'README.md'
        readme.write_text(readme.read_text()+'\n## n8n integration\n\nSee [setup, coverage and validation](integrations/n8n/README.md). This bundle includes the reviewed deployment patch; it is not a new published release.\n')
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with tarfile.open(args.out,'w:gz') as archive:
            archive.add(folder,arcname='.')
    print(json.dumps({'archive':str(args.out.resolve()),'source_sha256':actual}))


if __name__=='__main__':
    main()
