switch.strand = function(allele){

  allele = toupper(allele)

  intermediate = str_replace_all(allele, c("A" = "H", "T" = "J", "G" = "K", "C" = "L"))

  out_allele = str_replace_all(intermediate, c("H" = "T", "J" = "A", "K" = "C", "L" = "G"))

  return(out_allele)

}
