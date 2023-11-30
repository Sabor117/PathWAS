#' @export
flip.allele = function(ref_ld, alt_ld, ref_gen, alt_gen) {

  ### Initialise variable for each row

  output = 0

  if (ref_ld == ref_gen & alt_ld == alt_gen){

    ### If alleles in both files match, return 1

    output = 1

  }	else if (ref_ld == alt_gen & alt_ld == ref_gen){

    ### If alleles in both files are exact opposites, return -1

    output = -1

  }	else{

    ### Make sure that files aren't using alternate strands
    ### Switch alleles to opposing strand using switch_strand function

    ref_ld_switched = switch.strand(ref_ld)
    alt_ld_switched = switch.strand(alt_ld)

    if (ref_ld_switched == ref_gen & alt_ld_switched == alt_gen){

      ### If new switched alleles match, return 1

      output = 1

    }	else if (ref_ld_switched == alt_gen & alt_ld_switched == ref_gen){

      ### If new switched alleles are opposites, return -1

      output = -1

    }	else {

      ### If the alleles do not match then return 0 for QC

      output = 0

    }
  }

  return(output)
}
