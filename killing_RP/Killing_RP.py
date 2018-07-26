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

def match_variants(gwas, reference):
    """
    Groups `gwas` and `reference` by the row key (e.g. locus) and then
    compares the allele fields to determine whether a strand flip has occurred.

    Assumes `gwas` row key and `reference` row key are the same type.
    Assumes both `gwas` and `reference` have a row field `alleles` of
    type <array<str>>. `alleles` cannot be a row key of the dataset.
    """
    def has_field_of_type(source, name, dtype):
        return name in source.row and source[name].dtype == dtype

    if not has_field_of_type(gwas, 'alleles', hl.tarray(hl.tstr)):
        raise TypeError("'gwas' must have a row field 'alleles' with type <array<str>>")

    if not has_field_of_type(reference, 'alleles', hl.tarray(hl.tstr)):
        raise TypeError("'reference' must have a row field 'alleles' with type <array<str>>")

    if 'alleles' in gwas.key:
        raise TypeError("'alleles' cannot be a row key in 'gwas'.")

    if 'alleles' in reference.key:
        raise TypeError("'alleles' cannot be a row key in 'reference'.")        

    reference = reference.collect_by_key()

    matched = gwas.annotate(matches=hl.map(lambda x: x.annotate(match_alleles=add_strand_flip_annotation(
        x.alleles[0], 
        x.alleles[1], 
        gwas.alleles[0], 
        gwas.alleles[1])), reference[gwas.key].values))

    return matched
