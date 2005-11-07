joint.glm <- function(form.mean, form.disp, data,
                      family.mean = gaussian, family.disp = Gamma(link='log'),
                      eps = 10^(-6), iter.max = 100, maxit = 50, reml = TRUE){
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

  ##First convert formulas in ``real'' formula
  form.mean <- as.formula(form.mean)
  form.disp <- as.formula(form.disp)
  
  ##Set the initial weigth vector
  weights <- rep(1, length(data[,1]))
  ##Set an arbitrary value for the EQL beacause
  ##difference between two iterations will be a criteria
  ##for our test of convergence.
  dev.new <- Inf
  
  ##Introduce a logical object ``flag'' which specifie
  ##if convergence is reached or not.
  ##   TRUE    ->   convergence test OK
  ##   FALSE   ->   convergence test non OK
  ##Initialize it to FALSE to enter in the ``while'' loop
  flag <- FALSE

  ##Introduce also a counter for the number of
  ##iterations
  n.iter <- 0

  ##First, fit the GLM associated to the mean
  glm.mean <- glm(form.mean, family = family.mean, data = data,
                  maxit = maxit, weights = weights)
  
  while ( !flag ){

    ##We can then compute the dispersion response
    ##i.e. the deviance contribution
    glm.mean.alias <- glm.mean          #craation of an alias of
                                        #``glm.mean'' object
    glm.mean.alias$prior.weights <- rep(1, length(data[,1]))
                                        #Set the prior weights
                                        #to  1 to evaluate the deviance
                                        #contribution
    d <- residuals(glm.mean.alias, 'deviance')^2
    
    if (reml)
      d <- d / (1 - hat.glm(glm.mean))
    
    data.disp <- cbind(d = d, data[,-1])
    data.disp <- as.data.frame(data.disp)
    names(data.disp) <- c('d', names(data)[-1])

    weights.disp <- rep(1, length(d))
    if (reml)
      weights.disp <- 1 - hat.glm(glm.mean)
    
    glm.disp <- glm(form.disp, family = family.disp,
                    data = data.disp, weights = weights.disp,
                    maxit = maxit)
    
    gam <- glm.disp$coefficients
    phi <- glm.disp$fitted

    ##Update for the mean weights
    weights <- 1 / phi
    glm.mean <- glm(form.mean, family = family.mean, data = data,
                    maxit = maxit, weights = weights)

    dev.old <- dev.new
            
    dev.new <- eqd.joint(glm.mean,glm.disp)
    cat("Percentage of EQL variation is : ",
        abs( (dev.old - dev.new) / dev.old * 100),"\n")
    flag <- ( iter.max < n.iter) || abs( (dev.old - dev.new) / dev.new * 100) < eps
   
    n.iter <- n.iter + 1
    
  }

  return(list(glm.mean = glm.mean, glm.disp = glm.disp,
              iterations = n.iter, eqd = dev.new))
}

eqd.joint <- function(glm.mean, glm.disp, df.adj = TRUE){
  ##A function which compute the Extended Quasi-Deviance
  ##(EQD) for the Joint Model

  ##glm.mean : Object of class ``glm'' corresponding to
  ##           the mean component
  ##glm.disp : Object of class ``glm'' corresponding to
  ##           the dispersion component
  ##df.adj   : Logique. If ``TRUE'' adjustment on degree
  ##           of freedom is applied

  dispersion <- glm.disp$fitted              #dispersion parameter
                                             #estimate
  dev.res <- glm.disp$y
  
  fct.var <- glm.mean$family$variance        #Variance function

  if(df.adj)
    eqllik <- - dev.res / (2 * dispersion)  -
      1/2 * glm.mean$df.residual / length(dispersion) *
        log(2*pi*dispersion*fct.var(glm.mean$y))
  else
    eqllik <- -dev.res / (2 * dispersion) -
      1/2 * log(2*pi*dispersion*fct.var(glm.mean$y))
  
  
  eqd <- -2 * sum(dispersion * eqllik)

  return(eqd)
}

hat.glm <- function(glm){
  X <- as.matrix(glm$data[,-1])
  W <- diag(glm$weights)
  H <- sqrt(W) %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% sqrt(W)
  h <- diag(H)

  return(h)
}
