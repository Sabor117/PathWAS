lmp = function (modelobject) {

  if (class(modelobject) != "lm"){

    stop("Not an object of class 'lm'. Womp womp.")

  }

  f = summary(modelobject)$fstatistic
  p = pf(f[1], f[2], f[3],
         lower.tail = F)

  attributes(p) = NULL
  return(p)

}
