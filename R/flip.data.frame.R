flip.data.frame = function(df) {

  apply(df[2:5], 1, function(x) flip.allele(x[1], x[2], x[3], x[4]))

}
