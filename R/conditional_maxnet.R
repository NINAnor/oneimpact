#' @export
conditional_maxnet <-
  function(p, strata, data, f=maxnet::maxnet.formula(p, data), regmult=1.0,
           regfun=maxnet::maxnet.default.regularization, ...)
  {
    if (anyNA(data)) stop("NA values in data table. Please remove them and rerun.")

    mm <- model.matrix(f, data)
    reg <- regfun(p,mm) * regmult
    weights <- rep(1, length(p))
    glmnet::glmnet.control(pmin=1.0e-8, fdev=0)
    p_strat <- glmnet::stratifySurv(survival::Surv(rep(1, length(p)), p), strata)
    model <- glmnet::glmnet(x=mm, y=p_strat, family="cox", standardize=F, penalty.factor=reg, lambda=10^(seq(4,0,length.out=200))*sum(reg)/length(reg)*sum(p)/sum(weights), ...)
    class(model) <- c("maxnet", class(model))
    if (length(model$lambda) < 200) {
      msg <- "Error: glmnet failed to complete regularization path.  Model may be infeasible."
      if (!addsamplestobackground)
        msg <- paste(msg, " Try re-running with addsamplestobackground=T.")
      stop(msg)
    }
    bb <- model$beta[,200]
    model$betas <- bb[bb!=0]
    model$alpha <- 0
    ##rr <- maxnet::predict.maxnet(model, data[p==0, , drop = FALSE], type="exponent", clamp=F) #check this part...
    rr <- predict.maxnet(model, data[p==0, , drop = FALSE], type="exponent", clamp=F) #check this part...
    raw <- rr / sum(rr)
    model$entropy <- -sum(raw * log(raw))
    model$alpha <- -log(sum(rr))
    model$penalty.factor <- reg
    model$featuremins <- apply(mm, 2, min)
    model$featuremaxs <- apply(mm, 2, max)
    vv <- (sapply(data, class)!="factor")
    model$varmin <- apply(data[,vv, drop = FALSE], 2, min)
    model$varmax <- apply(data[,vv, drop = FALSE], 2, max)
    means <- apply(data[p==1,vv, drop = FALSE], 2, mean)
    majorities <- sapply(names(data)[!vv],
                         function(n) which.max(table(data[p==1,n, drop = FALSE])))
    names(majorities) <- names(data)[!vv]
    model$samplemeans <- unlist(c(means, majorities))
    model$levels <- lapply(data, levels)
    model
  }
