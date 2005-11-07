obs.vs.model <- function(glm.joint, plot.disp = FALSE, ...){
  ##Plot observations functions of predicted values

  ##glm.joint  : A list with two component. Both are objects of class
  ##             ``glm'' corresponding to the mean and dispersion
  ##             respectively.
  ##plot.disp : Logical. If ``TRUE'' error bar - i.e. +/- standard
  ##            deviation - are displayed on selected points.
  ##...       : Optional arguments to be passed to ``plot'' and
  #             ``abline'' functions.

  glm.mean <- glm.joint$glm.mean
  glm.disp <- glm.joint$glm.disp
  
  Model <- glm.mean$fitted
  Obs <- glm.mean$y
  
  plot(Model, Obs, ...)
  abline(a=0, b=1, ...)

  if (plot.disp){
    fct.var <- glm.mean$family$variance      #The variance function
    choice <- identify(x = Model, y = Obs, plot = FALSE)
    ecart.type <- sqrt(glm.disp$fitted * fct.var(glm.mean$fitted)) 
    arrows(Model[choice], Obs[choice], Model[choice] +
           ecart.type[choice], Obs[choice],
           angle= 90, length= 0.05)
    arrows(Model[choice], Obs[choice], Model[choice] -
           ecart.type[choice], Obs[choice],
           angle= 90, length= 0.05)
  }
}


rstand.vs.linpred <- function(glm, smooth = TRUE, ...){
  ##Plot standardized residuals functions of  the linear
  ##predictor

  ##glm    : Object of class ``glm''
  ##smooth : Logical. If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull
  ##         to detect departure from the model.
  ##...    : Optional arguments to be passed to the
  ##         ``plot'' and ``abline'' functions.

  res.stand <- residuals(glm, 'deviance') / sqrt(1 - hat.glm(glm))
  lin.pred <- glm$linear.predictors
  plot(lin.pred, res.stand, ...)
  abline( h =0, ...)

  if (smooth){
    lines(lowess(x = lin.pred, y = res.stand),
          col = 2, ...)
  }
}

res.vs.explvar <- function(glm, var, res = 'standard',
                           smooth = TRUE, ...){
  ##Plot residuals functions of an explanatory variable

  ##glm    : Object of class ``glm''
  ##var    : Character. The name of the explanatory variable
  ##res    : Should be 'standard', 'student', 'brut'
  ##         The residual type chosen.
  ##smooth : Logical. If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull to
  ##         detect departure from the model.
  ##...    : Optional arguments to be passed to 
  ##         ``plot'' and ``abline'' functions.

  if (!any(names(glm$data) == var))
    stop("Invalid explanatory variable !")
  else
    var <- as.vector(glm$data[[var]])

  res <- switch(res, 'standard' = residuals(glm,'deviance') /
                sqrt(1 - hat.glm(glm)),
                'student' = residuals(glm, 'deviance') /
                sqrt(summary(glm)$dispersion*(1-hat.glm(glm))),
                'brut' = residuals(glm, 'deviance'))
  plot(var, res, ...)
  abline( h = 0, ...)

  if ( smooth ){
    lines(lowess(x = var, y = res), col = 2, ...)
  }
}
  
absres.vs.fitted <- function(glm, res = 'standard',
                             smooth = TRUE, ...){
  ##Plot absolute residuals functions of predicted values

  ##glm    : Object of class``glm''
  ##res    : Should be 'standard', 'student', 'brut'
  ##         The residual type chosen.
  ##smooth : Logical. . If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull to
  ##         detect departure from the model.
  ##...    : Optional arguments to be passed to the ``plot''
  ##         function.


  res <- switch(res, 'standard' = residuals(glm,'deviance') /
                sqrt(1 - hat.glm(glm)),
                'student' = residuals(glm, 'deviance') /
                sqrt(summary(glm)$dispersion*(1-hat.glm(glm))),
                'brut' = residuals(glm, 'deviance'))

  fitted <- glm$fitted

  plot(fitted, abs(res), ...)

  if (smooth){
    lines(lowess(x = fitted, y = abs(res)), col = 2, ...)
  }
}

adjvar.vs.linpred <- function(glm, smooth = TRUE, ...){
  ##Plot the adjusted dependent variable functions of
  ##the linear predictor.

  ##glm    : Object of class``glm''
  ##smooth : Logical. . If ``TRUE'' response surface is estimated
  ##         thanks to the ``lowess'' function. It is helpfull to
  ##         detect departure from the model.
  ##...    : Optional arguments to be passed to the ``plot''
  ##         function.

  link <- glm$family$linkfun
  y <- glm$y
  lin.pred <- glm$linear.predictors
  z <- link(y)

  link <- glm$family$link
  switch(link,
         'identity' = { z <- z + y - glm$fitted },
         'log' = { z <- z + (y - glm$fitted) / glm$fitted },
         'inverse' = { z <- z - (y - glm$fitted) / glm$fitted^2 },
         'logit' = { z <- z + (y - glm$fitted) / (1 - glm$fitted)^2 },
         'sqrt' = { z <- z + 2 * (y - glm$fitted) / sqrt(glm$fitted) },
         '1/mu^2' = { z <- z - 2 * (y - glm$fitted) / glm$fitted^3 },
         'probit' = stop('Pas encore implémenté !!!'),
         'cloglog' = stop('Pas encore implémenté !!!')
         )
         
         
  plot(lin.pred, z, ...)
  abline(a=0, b=1, ...)

  if (smooth){
    lines(lowess(x = lin.pred, y = z), col = 2, ...)
  }

}

summplot.glm <- function(glm, var = NULL, res = 'standard',
                         which = 1:5,
                         ask = nb.fig < length(which) &&
                         dev.interactive(), smooth = TRUE,
                         ...){
  if (!is.numeric(which) || any(which < 1) || any(which > 5)) 
        stop("`which' must be in 1:5")
  if (is.null(var) && any(which == 2))
    stop("Select an explanatory variable please.")
  
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (show[1]) {
    rstand.vs.linpred(glm, smooth = smooth, ...)
  }
  if (show[2]) {
    res.vs.explvar(glm, var, res, smooth = smooth, ...)
  }
  if (show[3]) {
    absres.vs.fitted(glm, res, smooth = smooth, ...)
  }
  if (show[4]) {
    adjvar.vs.linpred(glm, smooth = smooth, ...)
  }
  if (show[5]) {
    qqglm(glm, ...)
  }
}

qqglm <- function(glm, ...){
  ##Produce a QQ-plot for the studentized residuals

  ##glm    : Object of class``glm''
  ##...    : Optional arguments to be passed to the ``qqnorm''
  ##         function.

  res <- residuals(glm, 'deviance') /
                sqrt((1-hat.glm(glm)) * summary(glm)$dispersion)
                
  qqnorm(res, ...)
  abline(a = 0, b = 1, ...)

}
  
