import pysam

from bespoke_bio.stage1.precomputed import coding_regions, phase_identity


def record(*, phased, explicit_ps=False, ps=None):
    header = pysam.VariantHeader()
    header.contigs.add('chr1', length=100)
    header.formats.add('GT', 1, 'String', 'Genotype')
    header.formats.add('PS', 1, 'Integer', 'Phase set')
    header.add_sample('HG002')
    row = header.new_record(contig='chr1', start=9, alleles=('A', 'C'))
    row.samples['HG002']['GT'] = (0, 1)
    row.samples['HG002'].phased = phased
    if explicit_ps:
        row.samples['HG002']['PS'] = ps
    return row


def test_import_keeps_declared_phase_without_inventing_unphased_links():
    assert phase_identity(record(phased=False), 'HG002', 'a'*64) is None
    implicit = phase_identity(record(phased=True), 'HG002', 'a'*64)
    assert implicit.endswith(':HG002:chr1:vcf-default')
    assert implicit != phase_identity(record(phased=True), 'HG002', 'b'*64)
    assert phase_identity(record(phased=True, explicit_ps=True), 'HG002', 'a'*64) is None
    assert phase_identity(record(phased=True, explicit_ps=True, ps=7), 'HG002', 'a'*64).endswith(':PS:7')


def test_coding_annotation_regions_include_flanks_and_preserve_full_source(tmp_path):
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chr1\n' + 'A'*100 + '\n')
    pysam.faidx(str(fasta))
    gtf = tmp_path / 'genes.gtf'
    gtf.write_text(
        'chr1\ttest\tCDS\t2\t9\t.\t+\t0\ttranscript_type "protein_coding";\n'
        'chr1\ttest\tstop_codon\t10\t12\t.\t+\t0\ttranscript_type "protein_coding";\n'
        'chr1\ttest\tCDS\t50\t60\t.\t+\t0\ttranscript_type "nonsense_mediated_decay";\n'
    )
    assert coding_regions(gtf, fasta, 5) == {'chr1': [(0, 17)]}
