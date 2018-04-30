import hail as hl
from hail.typecheck import *
from hail.expr.expressions.expression_typecheck import *

@typecheck(a=expr_str)
def flip_strand(a):
    alleles = hl.literal({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})
    return alleles.get(a)

@typecheck(reference_ref=expr_str,
           reference_alt=expr_str,
           ds_a1=expr_str,
           ds_a2=expr_str)
def add_strand_flip_annotation(reference_ref, reference_alt, ds_a1, ds_a2):
    """ Document me here :)
    """
    is_strand_ambig = hl.is_strand_ambiguous(ds_a1, ds_a2)
    ds_a1_flipped   = flip_strand(ds_a1)
    ds_a2_flipped   = flip_strand(ds_a2)
    is_snp = hl.is_snp(ds_a1, ds_a2)
    null = hl.null(hl.tbool)

    return (hl.case()
       .when((ds_a1 == reference_alt) & (ds_a2 == reference_ref), 
             hl.cond(is_strand_ambig, 
                     [hl.struct(swap = True,  flip = True), hl.struct(swap = False, flip = False)], 
                     [hl.struct(swap = False, flip = False)]))
       .when((ds_a1 == reference_ref) & (ds_a2 == reference_alt), 
             hl.cond(is_strand_ambig, 
                     [hl.struct(swap = True,  flip = False), hl.struct(swap = False, flip = True)],
                     [hl.struct(swap = True,  flip = False)]))
       .when((ds_a1_flipped == reference_alt) & (ds_a2_flipped == reference_ref) & is_snp, 
                     [hl.struct(swap = False, flip = True)])
       .when((ds_a1_flipped == reference_ref) & (ds_a2_flipped == reference_alt) & is_snp, 
                     [hl.struct(swap = True,  flip = True)])
       .default(hl.empty_array(hl.tstruct(swap=hl.tbool, flip=hl.tbool))))  

@typecheck(gwas_row_key=oneof(expr_locus(), expr_struct(), expr_str),
           gwas_alleles=expr_array(expr_str),
           reference_row_key=oneof(expr_locus(), expr_struct(), expr_str),
           reference_alleles=expr_array(expr_str))
def match_variants(gwas_row_key, gwas_alleles, reference_row_key, reference_alleles):
    """
    Assumes gwas row key and reference row key are same type and are keys of datasets already
    """
    assert(gwas_row_key.dtype == reference_row_key.dtype) # FIXME: Allow locus to be matched with tstruct(contig, pos)
    assert(gwas_row_key._indices == gwas_alleles._indices)
    assert(reference_row_key._indices == reference_alleles._indices)
    # FIXME: assert all row indices

    reference_alleles_fn = reference_alleles._ast.name

    gwas = gwas_row_key._indices.source 
    reference = reference_row_key._indices.source

    reference = reference.collect_by_key()

    matched_t = gwas.annotate(matches=hl.map(lambda x: x.annotate(match_alleles=add_strand_flip_annotation(
        x[reference_alleles_fn][0], 
        x[reference_alleles_fn][1], 
        gwas_alleles[0], 
        gwas_alleles[1])), reference[gwas_row_key].values))

    # t = matched_t.annotate(matches=hl.map(lambda x: x.annotate(flip=add_strand_flip_annotation(x[???], x[???], matched_t[???], matched_t[???])) , matched_t.matches))

    # t = t.select(t.locus, t.A1, t.A2, matches = hl.map(lambda x: x, t.matches))

    # Table(gwas_row_key, 
    #       matches[struct{
    #                       ref_index = int, 
    #                       ref_locus = locus, 
    #                       ref_rsid = str, 
    #                       ref_alleles = (ref,alt), 
    #                       allele_swap = bool, 
    #                       strand_flip = bool
    #                   }, ...
    #               ]
    #       )
    return matched_t
