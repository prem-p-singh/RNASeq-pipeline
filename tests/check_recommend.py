#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from recommend import recommend, handbook_backend

cfg = {"samples": {"seq_type": "rnaseq_paired"},
       "model": {"fixed_effects": "~ condition", "random_effects": None}}
a, b = recommend(cfg, 6), recommend(cfg, 300)
assert a["route"] == b["route"] == "salmon_bulk"
assert a["backend"] == b["backend"] == "limma_voom"
assert a["handbook_ids"] and a["policy_sha256"] and a["rule_ids"]
for backend in ("limma_voom", "edger_ql", "deseq2"):
    cfg["analysis"] = {"backend": backend}
    result = recommend(cfg)
    assert result["status"] == "selected" and result["backend"] == backend
cfg.pop("analysis")
cfg["model"]["random_effects"] = "(1|subject)"
assert recommend(cfg)["backend"] == "dream"
for analysis in ({"objectives": ["splicing"]}, {"umi": True}, {"backend": "deseq2"}):
    cfg["analysis"] = analysis
    assert recommend(cfg)["status"] == "blocked"
cfg.pop("analysis")
for seq in ("small_rna", "single_cell", "spatial", "long_read", "unknown"):
    cfg["samples"]["seq_type"] = seq
    assert recommend(cfg)["status"] == "blocked"
print("check_recommend.py: scope, dependence, traceability and sample-count invariance passed")
star_cfg = {'samples': {'seq_type': 'rnaseq_paired', 'expected_libtype': 'IU'},
            'reference': {'genome_fasta_url': 'file:///genome.fa'},
            'analysis': {'objectives': ['alignment']}}
assert recommend(star_cfg)['route'] == 'star_counts'
assert recommend(star_cfg)['backend'] is None
assert recommend(star_cfg)['count_treatment'] == 'raw_gene_counts'
for strand in ('A', 'U', None):
    star_cfg['samples']['expected_libtype'] = strand
    assert recommend(star_cfg)['status'] == 'blocked'
star_cfg['samples']['expected_libtype'] = 'IU'
star_cfg['analysis']['quantifier'] = 'salmon'
assert recommend(star_cfg)['status'] == 'blocked'
star_cfg['analysis'] = {'objectives': ['qc']}
star_cfg.pop('reference')
assert recommend(star_cfg)['route'] == 'fastp_qc'
assert recommend(star_cfg)['count_treatment'] is None
print('STAR objective selection, layout/strand validation and QC independence passed')
# Handbook default for new independent designs; 2 units needs an explicit choice.
assert [handbook_backend(n) for n in (None, 1, 2, 3, 12, 13)] == [None, None, None, 'deseq2', 'deseq2', 'limma_voom']
print('handbook backend defaults passed')
srna_cfg = {'samples': {'seq_type': 'srnaseq'},
            'analysis': {'objectives': ['gene_expression'], 'small_rna_profile': 'animal', 'small_rna_kit': 'TruSeq Small RNA'},
            'reference': {'mature_mirna_url': 'file:///mature.fa'}, 'organism': {'mirbase_prefix': 'hsa'}}
r = recommend(srna_cfg)
assert (r['status'], r['route'], r['count_treatment']) == ('selected', 'small_rna_animal', 'raw_small_rna_counts'), r
for key in ('small_rna_profile', 'small_rna_kit'):
    assert recommend({**srna_cfg, 'analysis': {**srna_cfg['analysis'], key: None}})['status'] == 'blocked'
assert recommend({**srna_cfg, 'organism': {}})['status'] == 'blocked'
assert recommend({**srna_cfg, 'reference': {}})['status'] == 'blocked'
print('animal small-RNA route selection passed')
assert r['qualification'] == 'experimental; no public-data check passed', r['qualification']
for kit, word in (('QIAseq miRNA', 'UMI'), ('No Such Kit', 'Unknown small-RNA kit')):
    b = recommend({**srna_cfg, 'analysis': {**srna_cfg['analysis'], 'small_rna_kit': kit}})
    assert b['status'] == 'blocked' and any(word in i for i in b['issues']), b['issues']
b = recommend({**srna_cfg, 'analysis': {**srna_cfg['analysis'], 'objectives': ['enrichment']},
               'downstream': {'run_go': True}})
assert b['status'] == 'blocked' and any('miRNA/locus identifiers' in i for i in b['issues']), b['issues']
print('small-RNA kit, qualification and enrichment blocks passed')
plant = {**srna_cfg, 'analysis': {**srna_cfg['analysis'], 'small_rna_profile': 'plant'},
         'reference': {'genome_fasta_url': 'file:///g.fa'}, 'organism': {}}
r = recommend(plant)
assert (r['status'], r['route'], r['count_treatment']) == ('selected', 'small_rna_plant', 'raw_small_rna_counts'), r
assert 'ShortStack' in r['reasons'][0], r['reasons']
assert recommend({**plant, 'reference': {**plant['reference'], 'mature_mirna_url': 'file:///m.fa'},
                  'organism': {'mirbase_prefix': 'ath'}})['status'] == 'selected'
b = recommend({**plant, 'reference': {}})
assert b['status'] == 'blocked' and any('genome_fasta_url' in i for i in b['issues']), b['issues']
b = recommend({**plant, 'reference': {**plant['reference'], 'mature_mirna_url': 'file:///m.fa'}})
assert b['status'] == 'blocked' and any('mirbase_prefix' in i for i in b['issues']), b['issues']
print('plant small-RNA route selection passed')
