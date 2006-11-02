obs.vs.model <- function(x.joint, plot.disp = FALSE, ...){

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  ##Plot observations functions of predicted values

  ##x.joint  : A list with two component. Both are objects of class
  ##             ``glm'' corresponding to the mean and dispersion
  ##             respectively.
  ##plot.disp : Logical. If ``TRUE'' error bar - i.e. +/- standard
  ##            deviation - are displayed on selected points.
  ##...       : Optional arguments to be passed to ``plot'' and
  #             ``abline'' functions.

  mod.mean <- x.joint$mod.mean
  mod.disp <- x.joint$mod.disp
  
  Model <- mod.mean$fitted
  Obs <- mod.mean$y
  
  plot(Model, Obs, ...)
  abline(a=0, b=1, ...)

  if (plot.disp){
    fct.var <- mod.mean$family$variance      #The variance function
    choice <- identify(x = Model, y = Obs, plot = FALSE)
    ecart.type <- sqrt(mod.disp$fitted * fct.var(mod.mean$fitted)) 
    arrows(Model[choice], Obs[choice], Model[choice] +
           ecart.type[choice], Obs[choice],
           angle= 90, length= 0.05)
    arrows(Model[choice], Obs[choice], Model[choice] -
           ecart.type[choice], Obs[choice],
           angle= 90, length= 0.05)
  }
}


rstand.vs.linpred <- function(x, smooth = TRUE, ...){

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  ##Plot standardized residuals functions of  the linear
  ##predictor

  ##x    : A fitted joint model object
  ##smooth : Logical. If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull
  ##         to detect departure from the model.
  ##...    : Optional arguments to be passed to the
  ##         ``plot'' and ``abline'' functions.

  res.stand <- residuals(x, 'deviance') / sqrt(1 - hat.glm(x))
  lin.pred <- x$linear.predictors
  plot(lin.pred, res.stand, ...)
  abline( h =0, ...)

  if (smooth){
    lines(lowess(x = lin.pred, y = res.stand),
          col = 2, ...)
  }
}

res.vs.explvar <- function(x, var, res = 'standard',
                           smooth = TRUE, ...){

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  ##Plot residuals functions of an explanatory variable

  ##x    : Object of class ``glm''
  ##var    : Character. The name of the explanatory variable
  ##res    : Should be 'standard', 'student', 'brut'
  ##         The residual type chosen.
  ##smooth : Logical. If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull to
  ##         detect departure from the model.
  ##...    : Optional arguments to be passed to 
  ##         ``plot'' and ``abline'' functions.

  if (!any(names(x$data) == var))
    stop("Invalid explanatory variable !")
  else
    var <- as.vector(x$data[[var]])

  res <- switch(res, 'standard' = residuals(x,'deviance') /
                sqrt(1 - hat.glm(x)),
                'student' = residuals(x, 'deviance') /
                sqrt(summary(x)$dispersion*(1-hat.glm(x))),
                'brut' = residuals(x, 'deviance'))
  plot(var, res, ...)
  abline( h = 0, ...)

  if ( smooth ){
    lines(lowess(x = var, y = res), col = 2, ...)
  }
}
  
absres.vs.fitted <- function(x, res = 'standard',
                             smooth = TRUE, ...){

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  ##Plot absolute residuals functions of predicted values

  ##x    : Object of class``glm''
  ##res    : Should be 'standard', 'student', 'brut'
  ##         The residual type chosen.
  ##smooth : Logical. . If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull to
  ##         detect departure from the model.
  ##...    : Optional arguments to be passed to the ``plot''
  ##         function.


  res <- switch(res, 'standard' = residuals(x,'deviance') /
                sqrt(1 - hat.glm(x)),
                'student' = residuals(x, 'deviance') /
                sqrt(summary(x)$dispersion*(1-hat.glm(x))),
                'brut' = residuals(x, 'deviance'))

  fitted <- x$fitted

  plot(fitted, abs(res), ...)

  if (smooth){
    lines(lowess(x = fitted, y = abs(res)), col = 2, ...)
  }
}

adjvar.vs.linpred <- function(x, smooth = TRUE, ...){

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  ##Plot the adjusted dependent variable functions of
  ##the linear predictor.

  ##x    : Object of class``glm''
  ##smooth : Logical. . If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull to
  ##         detect departure from the model.
  ##...    : Optional arguments to be passed to the ``plot''
  ##         function.

  link <- x$family$linkfun
  y <- x$y
  lin.pred <- x$linear.predictors
  z <- link(y)

  link <- x$family$link
  switch(link,
         'identity' = { z <- z + y - x$fitted },
         'log' = { z <- z + (y - x$fitted) / x$fitted },
         'inverse' = { z <- z - (y - x$fitted) / x$fitted^2 },
         'logit' = { z <- z + (y - x$fitted) / (1 - x$fitted)^2 },
         'sqrt' = { z <- z + 2 * (y - x$fitted) / sqrt(x$fitted) },
         '1/mu^2' = { z <- z - 2 * (y - x$fitted) / x$fitted^3 },
         'probit' = stop('Not implemented yet !!!'),
         'cloglog' = stop('Not implemented yet !!!')
         )
         
         
  plot(lin.pred, z, ...)
  abline(a=0, b=1, ...)

  if (smooth){
    lines(lowess(x = lin.pred, y = z), col = 2, ...)
  }

}

plot.joint <- function(x, comp = "mean", var = NULL,
                       res = 'standard', which = 1:5,
                       ask = nb.fig < length(which) &&
                       dev.interactive(), smooth = TRUE, ...){

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  if (!is.numeric(which) || any(which < 1) || any(which > 5)) 
        stop("`which' must be in 1:5")
  if (is.null(var) && any(which == 2))
    stop("Select an explanatory variable please.")

  if (comp == "mean")
    x <- x$mod.mean

  else
    x <- x$mod.disp
  
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (show[1]) {
    rstand.vs.linpred(x, smooth = smooth, ...)
  }
  if (show[2]) {
    res.vs.explvar(x, var, res, smooth = smooth, ...)
  }
  if (show[3]) {
    absres.vs.fitted(x, res, smooth = smooth, ...)
  }
  if (show[4]) {
    adjvar.vs.linpred(x, smooth = smooth, ...)
  }
  if (show[5]) {
    qqglm(x, ...)
  }
}

qqglm <- function(x, ...){

  warning("Package JointGLM is outdated. Please use the JointModeling packgage.\n")
  
  ##Produce a QQ-plot for the studentized residuals

  ##x    : A fitted joint model object
  ##...  : Optional arguments to be passed to the ``qqnorm''
  ##       function.

  res <- residuals(x, 'deviance') /
                sqrt((1-hat.glm(x)) * summary(x)$dispersion)
                
  qqnorm(res, ...)
  abline(a = 0, b = 1, ...)

}
  
