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
assert r['qualification'] == 'validated; public-data check passed', r['qualification']
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
tag = {'samples': {'seq_type': 'tagseq', 'expected_libtype': None},
       'analysis': {'objectives': ['gene_expression'], 'quantifier': 'auto', 'kit': 'QuantSeq FWD-UMI'},
       'reference': {'genome_fasta_url': 'file:///g.fa', 'gtf_url': 'file:///a.gtf'}}
r = recommend(tag)
assert r['route'] == 'star_tagseq' and r['count_treatment'] == 'umi_gene_counts' and not r['issues'], r
tag['analysis']['kit'] = 'QuantSeq REV'
assert recommend(tag)['count_treatment'] == 'raw_gene_counts'
tag['analysis']['kit'] = 'Mystery'
assert any('Unknown TAG-seq kit' in i for i in recommend(tag)['issues'])
tag['analysis']['kit'] = None
assert recommend(tag)['route'] == 'salmon_end_tag_provisional', 'legacy tagseq must not change'
print('TAG-seq route selection passed')
tag['analysis'].update(kit='QuantSeq FWD', objectives=['alignment', 'gene_expression'])
r = recommend(tag)
assert (r['status'], r['route']) == ('selected', 'star_tagseq'), r['issues']
for bad, word in (({**tag, 'reference': {'genome_fasta_url': 'file:///g.fa'}}, 'gtf_url'),
                  ({**tag, 'samples': {'seq_type': 'rnaseq_single'}}, 'seq_type tagseq only'),
                  ({**tag, 'analysis': {**tag['analysis'], 'quantifier': 'salmon'}}, 'STAR route')):
    assert any(word in i for i in recommend(bad)['issues']), recommend(bad)['issues']
print('TAG-seq kit route blocks passed')
dual = {'samples': {'seq_type': 'rnaseq_paired'},
        'dual_organism': {'enabled': True,
                          'host': {'label': 'grape', 'genome_fasta_url': 'file:///h.fa', 'gtf_url': 'file:///h.gtf'},
                          'pathogen': {'label': 'GRBV', 'genome_fasta_url': 'file:///p.fa', 'gtf_url': 'file:///p.gtf'}}}
# Each organism runs as its own project afterwards, so the parent refuses whatever a child would refuse.
star = {**dual, 'analysis': {'quantifier': 'star'}, 'samples': {'seq_type': 'rnaseq_paired', 'expected_libtype': 'ISR'}}
salmon = {**dual, 'dual_organism': {**dual['dual_organism'], **{
    role: {**dual['dual_organism'][role], 'transcriptome_fasta_url': f'file:///{role}.rna.fa'} for role in ('host', 'pathogen')}}}
for valid, qualified in ((star, 'validated; public-data check passed'), (salmon, 'validated; public-data check passed')):
    r = recommend(valid)
    assert (r['status'], r['route'], r['backend'], r['count_treatment']) == ('selected', 'dual_assignment', None, None), r
    assert not r['issues'] and r['qualification'] == qualified, r
enrich = {**star, 'dual_organism': {**dual['dual_organism'], 'pathogen': {**dual['dual_organism']['pathogen'], 'objectives': ['gene_expression', 'enrichment']}}}
coexpr = {**star, 'dual_organism': {**dual['dual_organism'], 'pathogen': {**dual['dual_organism']['pathogen'], 'objectives': ['gene_expression', 'coexpression']}}}
nogtf = {**star, 'dual_organism': {**dual['dual_organism'], 'pathogen': {'label': 'GRBV', 'genome_fasta_url': 'file:///p.fa'}}}
half = {**salmon, 'dual_organism': {**salmon['dual_organism'], 'pathogen': dual['dual_organism']['pathogen']}}
for bad, words in (({**star, 'samples': {'seq_type': 'rnaseq_paired'}}, ('host: ', 'pathogen: ', 'expected_libtype')),
                   (half, ('dual_organism.pathogen.transcriptome_fasta_url',)),
                   (nogtf, ('dual_organism.pathogen.gtf_url',)),
                   (enrich, ('pathogen: Enrichment objective',)), (coexpr, ('pathogen: Coexpression objective',))):
    b = recommend(bad)
    assert b['status'] == 'blocked' and b['route'] == 'dual_assignment', b
    assert all(any(w in i for i in b['issues']) for w in words), (words, b['issues'])
assert not any(i.startswith('host: ') for i in recommend(enrich)['issues']), recommend(enrich)['issues']
b = recommend({**dual, 'samples': {'seq_type': 'tagseq'}})
assert b['status'] == 'blocked' and any('tagseq' in i for i in b['issues']), b['issues']
off = recommend({**dual, 'dual_organism': {**dual['dual_organism'], 'enabled': False}})
assert (off['route'], off['backend']) == ('salmon_bulk', 'limma_voom'), off
for bad, word in (({**dual, 'samples': {'seq_type': 'rnaseq_pared'}}, 'No implemented route'),
                  ({**dual, 'analysis': {'umi': True}}, 'UMI'),
                  ({**dual, 'analysis': {'kit': 'QuantSeq FWD'}}, 'analysis.kit'),
                  ({**dual, 'analysis': {'objectives': ['splicing']}}, 'splicing'),
                  ({**dual, 'analysis': {'objectives': 'qc'}}, 'nonempty list')):
    b = recommend(bad)
    assert b['status'] == 'blocked' and b['route'] == 'dual_assignment' and any(word in i for i in b['issues']), (word, b['issues'])
b = recommend({**star, 'analysis': {'quantifier': 'star', 'objectives': ['qc']}})
assert b['status'] == 'blocked' and any('quantification of both organisms' in i for i in b['issues']), b['issues']
r = recommend(star)
assert r['status'] == 'selected' and not r['issues'], r
print('dual-organism route selection passed')
