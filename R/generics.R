fitjoint <- function(model, form.mean, form.disp, data, family.mean = gaussian,
                     family.disp = Gamma(link='log'), eps = 10^(-6),
                     iter.max = 100, maxit = 50, reml = TRUE){

  ##This function will fit the Joint Model through the framework describe
  ##in [McCullagh1989]
  ##By default, it fit a gaussian distribution for the mean component
  ##and a Gamma distribution with log link function for the dispersion
  ##component

  ##form.mean     : a R ``formula'' put in a character type way which
  ##                specify the linear model for the mean
  ##form.disp     : a R ``formula'' put in a character type way which
  ##                specify the linear model for the dispersion
  ##data          : a ``data.frame'' corresponding to the data.
  ##                The first column correspond to the response
  ##                observations, while others to the explicative
  ##                variables.
  ##family.mean   : The ``family'' for the mean component.
  ##family.disp   : The ``family'' for the dispersion component.
  ##eps           : Optional Numeric. It specifies the precision
  ##                for the convergence test.
  ##iter.max      : The maximum number of iteration to fit the Joint
  ##                Model.
  ##maxit         : The maximum number of iteration to fit a generic
  ##                GLM - see function ``glm.fit''
  ##reml          : Logical. If TRUE, the Restricted Maximum Likelihood
  ##                Estimation is used.

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  if (!(model %in% c("glm", "gam")))
    stop("``model'' should be one of ``glm'' or ``gam''.")
  
  if (model == "glm")
    mod <- joint.glm(form.mean, form.disp, data, family.mean = family.mean,
                     family.disp = family.disp, eps = eps,
                     iter.max = iter.max, maxit = maxit, reml = reml)

  if (model == "gam")
    mod <- joint.gam(form.mean, form.disp, data, family.mean = family.mean,
                     family.disp = family.disp, eps = eps,
                     iter.max = iter.max, maxit = maxit)

  class(mod) <- c("joint", model)
  return(mod)
}

print.joint <- function(x, ...){
  ##A method function for printing ``joint'' objects.
  x.mean <- x$mod.mean
  x.disp <- x$mod.disp

  cat("EQL:", x$eql,"\n")
  cat("\n\tThe Mean Component\n")
  print(x.mean)
  cat("\n\tThe Dispersion Component\n")
  print(x.disp)
}
  
summary.joint <- function(object, ...){
  ##A method function to summarize ``joint'' objects.
  x.mean <- x$mod.mean
  x.disp <- x$mod.disp

  cat("EQL:", x$eql,"\n")
  cat("\n\tThe Mean Component\n")
  summary(x.mean, ...)
  cat("\n\tThe Dispersion Component\n")
  summary(x.disp, ...)
}
