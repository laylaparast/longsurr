

estimate_surrogate_value <- function(y_t, y_c, X_t, X_c, method = c('gam', 'linear', 'kernel'), k = 3, var = FALSE, bootstrap_samples = 50, alpha = 0.05) {
  if (method == 'linear') {
    Deltahat_S <- estimate_linear(y_t, y_c, X_t, X_c)
  } else if (method == 'gam') {
    Deltahat_S <- estimate_gam(y_t, y_c, X_t, X_c)
  } else if (method == 'kernel') {
    Deltahat_S <- estimate_kernel(y_t, y_c, X_t, X_c, k)
  }
  Deltahat <- mean(y_t) - mean(y_c)
  if (var) {
    boot_ests <- purrr::map(1:bootstrap_samples, boot_fn, method, k, y_t, y_c, X_t, X_c)
    boot_ests <- dplyr::bind_rows(boot_ests)
    boot_se <- summarise_all(boot_ests, sd, na.rm = TRUE)
    boot_ci_l <- summarise_all(boot_ests, quantile, alpha/2, na.rm = TRUE)
    boot_ci_h <- summarise_all(boot_ests, quantile, 1-alpha/2, na.rm = TRUE)
  } else {
    boot_se <- boot_ci_l <- boot_ci_h <- 
      tibble(Deltahat = NA,
             Deltahat_S = NA,
             R = NA,
             Deltahat_S_se = NA,
             Deltahat_S_ci_l = NA,
             Deltahat_S_ci_h = NA,
             R_se = NA,
             R_ci_l = NA,
             R_ci_h = NA)
  }
  if(var) {tibble::tibble(
    Deltahat = Deltahat,
    Deltahat_S = Deltahat_S,
    R = 1 - Deltahat_S/Deltahat,
    Deltahat_S_se = boot_se$Deltahat_S,
    Deltahat_S_ci_l = boot_ci_l$Deltahat_S,
    Deltahat_S_ci_h = boot_ci_h$Deltahat_S,
    R_se = boot_se$R,
    R_ci_l = boot_ci_l$R,
    R_ci_h = boot_ci_h$R
  )
  }
 if(!var) {tibble::tibble(
    Deltahat = Deltahat,
    Deltahat_S = Deltahat_S,
    R = 1 - Deltahat_S/Deltahat,
  )
  }
}
  
boot_fn <- function(b, method, k, y_t, y_c, X_t, X_c) {
    n1 <- length(y_t)
    n0 <- length(y_c)
    
    ind_t <- sample(1:n1, replace = TRUE)
    ind_c <- sample(1:n0, replace = TRUE)
    
    boot_yt <- y_t[ind_t]
    boot_yc <- y_c[ind_c]
    boot_Xt <- X_t[ind_t,]
    boot_Xc <- X_c[ind_c,]
    estimate_surrogate_value(
      boot_yt, 
      boot_yc, 
      boot_Xt, 
      boot_Xc, 
      method, 
      k, 
      bootstrap_samples = 0)
}




estimate_linear <- function(y_t, y_c, X_t, X_c, verbose = FALSE) {
  stopifnot(length(y_t) == nrow(X_t))
  stopifnot(length(y_c) == nrow(X_c))
  stopifnot(ncol(X_t) == ncol(X_c))
  
  lin_fit <- pfr(y_t ~ lf(X_t))
  lin_yhat = predict(lin_fit, newdata = list(X_t = X_c), type = 'response')
  
  k_lin_fit <- Rsurrogate::R.s.estimate(sone = predict(lin_fit, type = 'response'),
                            szero = lin_yhat,
                            yone = y_t,
                            yzero = y_c, var = FALSE, warn.support = TRUE)
  if (verbose) {
    unsm_delta_s <- mean(lin_yhat) - mean(y_c)
    print(glue::glue("Unsmoothed linear result is {unsm_delta_s}; smoothed linear results is {k_lin_fit$delta.s}."))
  }
  k_lin_fit$delta.s
}

estimate_kernel <- function(y_t, y_c, X_t, X_c, k = 3) {
  stopifnot(length(y_t) == nrow(X_t))
  stopifnot(length(y_c) == nrow(X_c))
  stopifnot(ncol(X_t) == ncol(X_c))
  
  fdX_t <- fdata(X_t)
  fdX_c <- fdata(X_c)
  
par.fda.usc<-list()
par.fda.usc$verbose <- FALSE
par.fda.usc$trace <- FALSE
par.fda.usc$warning <- FALSE
par.fda.usc$ncores <- 1 
par.fda.usc$int.method <- "TRAPZ"
par.fda.usc$eps <- as.double(.Machine[[1]]*10)
fda.usc::ops.fda.usc(ncores=1)

kernel_fit <- fregre.np.cv(fdataobj = fdX_t, y = y_t, 
                           metric = make_semimetric_pca(k), type.CV = "CV.S") 
  kernel_yhat = predict(kernel_fit, fdX_c)
  kernel_deltahat_s <- mean(kernel_yhat) - mean(y_c)
  kernel_deltahat_s
}

fpca <- function(ds, ycol, tcol, idcol, options = list(plot = FALSE, methodBwCov = 'GCV')) {
  # browser()
  list_data <- dataframe_to_list(ds, ycol, tcol, idcol)
  fpca_result <- with(list_data, fdapace::FPCA(y_list, t_list, options))
  score_ds <- data.frame(fpca_result$xiEst)
  colnames(score_ds) <- stringr::str_c('xi.', 1:ncol(score_ds))
  score_ds$id <- names(list_data$t_list)
  xi <- fpca_result$xiEst
  phi <- fpca_result$phi
  mu <- fpca_result$mu
  yh <- t(t(xi %*% t(phi)) + mu)
  yh_ds <- data.frame(yh)
  colnames(yh_ds) <- stringr::str_c('yhat.', 1:ncol(yh_ds))
  yh_ds$id <- names(list_data$t_list)
  # phi_pl <- plot_phi(fpca_result)
  # mu_pl <- plot_mu(fpca_result)
  phi_hat <- approx_phi(fpca_result$phi, fpca_result$workGrid, ds %>% select(tcol) %>% unlist)
  muhat <- approx_mu(fpca_result$mu, fpca_result$workGrid, ds %>% select(tcol) %>% unlist)
  phi_ds <- data.frame(id = ds %>% select(idcol), tt = ds %>% select(tcol), mu = muhat, phi = phi_hat)

  Lambda <- ifelse(length(fpca_result$lambda) == 1,
                    fpca_result$lambda,
                    diag(fpca_result$lambda))
  ids <- unique(ds$id)
  nn <- length(ids)
  Sigma <- fpca_result$smoothedCov
  diag(Sigma) <- ifelse(diag(Sigma) < 0, 1e-6, diag(Sigma))
  # var_l <- list()
  # for (i in 1:nn) {
  #   id_i <- ids[i]
  #   ds_i <- ds %>% filter(id == id_i)
  #   sigma2 <- fpca_result$sigma2
  #   Sigma_i <- approx_cov(Sigma,
  #                         fpca_result$workGrid,
  #                         ds_i %>% select(tcol) %>% unlist) +
  #     diag(nrow(ds_i))*sigma2
  #   phi_i <- phi_ds %>%
  #     filter(id == id_i) %>%
  #     select(contains('phi')) %>%
  #     as.matrix
  #   H <- phi_i %*% Lambda
  #   Omega <- Lambda - t(H) %*% solve(Sigma_i) %*% H
  #   # yh_var <- phi_i %*% Omega %*% t(phi_i)
  #   var_l[[i]] <- data.frame(id = as.character(i), t(diag(Omega)))
  #   colnames(var_l[[i]])[-1] <- paste0('xi.var.', 1:ncol(phi_i))
  # }
  # browser()
  # var_ds <- bind_rows(var_l)
  score_ds <- score_ds #%>% inner_join(var_ds)


  list(score_ds = score_ds, phi_ds = phi_ds, yh_ds = yh_ds, fpca_result = fpca_result)
}


presmooth_data <- function(obs_data, ...) {
  a = id = tt = tp = X = type = t_n = NULL
  treatment_arm <- obs_data %>%
    filter(a == 1) %>%
    arrange(id, tt)
  control_arm <- obs_data %>%
    filter(a == 0) %>%
    arrange(id, tt)
  n_trt <- treatment_arm %>%
    count(id) %>%
    nrow
  n_ctrl <- control_arm %>%
    count(id) %>%
    nrow
  
  
  trt_fpc_fit <- fpca(ds = treatment_arm, ycol = 'x', tcol = 'tt', idcol = 'id', ...)
  ctrl_fpc_fit <- fpca(ds = control_arm, ycol = 'x', tcol = 'tt', idcol = 'id', ...)
  
  
  times_1 <- tibble(tt = trt_fpc_fit$fpca_result$workGrid,
                    t_n = rank(tt))
    
  times_0 <- tibble(tt = ctrl_fpc_fit$fpca_result$workGrid,
                    t_n = rank(tt))
  
# browser()
trt_yh <- trt_fpc_fit$yh_ds %>%
  gather(tp, X, -id)
ctrl_yh <- ctrl_fpc_fit$yh_ds %>%
  gather(tp, X, -id)
  trt_xhat <- trt_yh %>%
    mutate(id = as.integer(id),
           # tt = rep(seq(-1, 1, length = 51), each = n_trt),
           t_n = as.numeric(stringr::str_remove(tp, 'yhat\\.')),
           type = 'estimated') %>%
    inner_join(times_1)

  ctrl_xhat <- ctrl_yh %>%
    mutate(id = as.integer(id),
           t_n = as.numeric(stringr::str_remove(tp, 'yhat\\.')),
           type = 'estimated') %>%
    inner_join(times_0)
  # browser()
  trt_xhat_wide <- trt_xhat %>%
    dplyr::select(-tp, -type, -t_n) %>%
    spread(tt, X) %>%
    dplyr::select(-id) %>%
    as.matrix
  colnames(trt_xhat_wide) <- colnames(trt_xhat_wide)
  rownames(trt_xhat_wide) <- trt_xhat$id %>% unique

  ctrl_xhat_wide <- ctrl_xhat %>%
    dplyr::select(-tp, -type, -t_n) %>%
    spread(tt, X) %>%
    dplyr::select(-id) %>%
    as.matrix
  colnames(ctrl_xhat_wide) <- colnames(ctrl_xhat_wide)
  rownames(ctrl_xhat_wide) <- ctrl_xhat$id %>% unique

  list(X_t = trt_xhat_wide, X_c = ctrl_xhat_wide)
}


approx_phi <- function(phi_mat, t_grid, new_t) {
  apply(phi_mat, 2, function(p) {
    approx(t_grid, p, xout = new_t, rule = 2)$y
  })
}

approx_mu <- function(mu, t_grid, new_t) {
  approx(t_grid, mu, xout = new_t, rule = 2)$y
}


estimate_gam <- function(y_t, y_c, X_t, X_c, verbose = FALSE) {
  stopifnot(length(y_t) == nrow(X_t))
  stopifnot(length(y_c) == nrow(X_c))
  stopifnot(ncol(X_t) == ncol(X_c))
  
  fgam_fit <- refund::fgam(y_t ~ af(X_t))
  fgam_yhat = predict(fgam_fit, newdata = list(X_t = X_c), type = 'response')  
  k_fgam_fit <- Rsurrogate::R.s.estimate(sone = predict(fgam_fit, type ='response'),
                             szero = fgam_yhat,
                             yone = y_t,
                             yzero = y_c, var = FALSE, warn.support = TRUE)
  if (verbose) {
    unsm_delta_s <- mean(fgam_yhat) - mean(y_c)
    print(glue::glue("Unsmoothed GAM result is {unsm_delta_s}; smoothed GAM results is {k_fgam_fit$delta.s}."))
  }
  k_fgam_fit$delta.s
}



dataframe_to_list <- function(ds, ycol, tcol, idcol) {
  id <- dplyr::select(ds, idcol)
  y_list <- plyr::dlply(ds, .variables = idcol, function(x) {
    unlist(dplyr::select(x, ycol))
  })
  t_list <- plyr::dlply(ds, .variables = idcol, function(x) {
    unlist(dplyr::select(x, tcol))
  })
  list(y_list = y_list, t_list = t_list, id = id)
}

get_fpc_scores <- function(ds, ycol, tcol, idcol, options = list(plot = TRUE), return_eigenfunction = FALSE) {
  list_data <- dataframe_to_list(ds, ycol, tcol, idcol)
  fpca_result <- with(list_data, fdapace::FPCA(y_list, t_list, options))
  score_ds <- data.frame(fpca_result$xiEst)
  colnames(score_ds) <- stringr::str_c('xi.', 1:ncol(score_ds))
  score_ds$id <- names(list_data$t_list)
  
  if (return_eigenfunction) {
    phi <- fpca_result$phi
    phi_hat = apply(fpca_result$phi, 2, function(p) {
    approx(fpca_result$workGrid, p, xout = ds %>% select(tcol) %>% unlist, rule = 2)$y})
    phi_ds <- data.frame(id = ds %>% select(idcol), tt = ds %>% select(tcol), phi = phi_hat)
    list(score_ds = score_ds,
         eigenfunction_ds = phi_ds,
         fpca_result = fpca_result
    )
  } else {
    list(score_ds = score_ds)
  }
}

make_semimetric_pca <- function(k) {
  function(x, y) semimetric.pca(x, y, q = k)
}


predict.fgam <- function (object, newdata, type = "response", se.fit = FALSE,
                          terms = NULL, PredOutOfRange = FALSE, ...)
{
  # browser()
  call <- match.call()
  string <- NULL
  if (!missing(newdata)) {
    nobs <- nrow(as.matrix(newdata[[1]]))
    
    stopifnot(length(unique(sapply(newdata, function(x) ifelse(is.matrix(x),
                                                               nrow(x), length(x))))) == 1)
    gamdata <- list()
    varmap <- sapply(names(object$fgam$labelmap), function(x) all.vars(formula(paste("~",
                                                                                     x))))
    for (cov in names(newdata)) {
      trms <- which(sapply(varmap, function(x) any(grep(paste("^",
                                                              cov, "$", sep = ""), x))))
      J <- ncol(as.matrix(newdata[[cov]]))
      if (length(trms) != 0) {
        for (trm in trms) {
          is.af <- trm %in% object$fgam$where$where.af
          is.lf <- trm %in% object$fgam$where$where.lf
          if (is.af) {
            af <- object$fgam$ft[[grep(paste(cov, "[,\\)]",
                                             sep = ""), names(object$fgam$ft))]]
            
            if(J!=length(af$xind) & type!='lpmatrix'){
              stop(paste('Provided data for functional covariate',cov,'does not have same observation times as original data',sep=''))
            }
            L <- matrix(af$L[1, ], nobs, J, byrow = T)
            tmat <- matrix(af$xind,nobs,J,byrow=TRUE)
            if (grepl(paste(cov, "\\.[ot]mat", sep = ""),
                      deparse(af$call$x))) {
              if (length(attr(newdata, "L"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new L matrix of quadrature weights only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "L")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied L matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                    
                  }else{
                    L <- as.vector(attr(newdata, "L"))
                  }
                }
              }
              if (length(attr(newdata, "tmat"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new tmat matrix of observation times only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "tmat")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied tmat matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                  }else{
                    tmat <- as.vector(attr(newdata, "tmat"))
                  }
                }
              }
              if (PredOutOfRange) {
                newdata[[cov]][newdata[[cov]]>af$Xrange[2]]  <- af$Xrange[2]
                newdata[[cov]][newdata[[cov]]<af$Xrange[1]]  <- af$Xrange[1]
              }
              if (!is.null(af$presmooth)) {
                if(type=='lpmatrix' & J!=length(af$xind)){
                  warning('Presmoothing of new functional covariates is only implemented for when new covariates observed at same time points as original data. No presmoothing of new covariates done.')
                } else if (is.logical(af$presmooth)) {
                  # af_old term
                  if (af$presmooth) {
                    newXfd <- fd(tcrossprod(af$Xfd$y2cMap,
                                            newdata[[cov]]), af$Xfd$basis)
                    newdata[[cov]] <- t(eval.fd(af$xind, newXfd))
                  }
                } else {
                  # af term
                  newdata[[cov]] <- af$prep.func(newX = newdata[[cov]])$processed
                }
              }
              # if (af$Qtransform) {
              #   if(type=='lpmatrix' & J!=length(af$xind)){
              #     stop('Prediction with quantile transformation only implemented for when new data observation times match original data observation times')
              #   }
              #   for (i in 1:nobs) {
              #     newdata[[cov]][i, ] <- mapply(function(tecdf,
              #                                            x) {
              #       tecdf(x)
              #     }, tecdf = af$ecdflist, x = newdata[[cov]][i,
              #                                                ])
              #   }
              # }
              if(type=='lpmatrix') newdata[[cov]] <- as.vector(newdata[[cov]])
              gamdata[[paste(cov, ".omat", sep = "")]] <- newdata[[cov]]
              gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
              gamdata[[paste("L.", cov, sep = "")]] <- L
            }
          }
          if (is.lf) {
            lf <- object$fgam$ft[[grep(paste(cov, "[,\\)]",
                                             sep = ""), names(object$fgam$ft))]]
            if(J!=length(lf$xind) & type!='lpmatrix'){
              stop(paste('Provided data for functional covariate',cov,'does not have same observation times as original data',sep=''))
            }
            L <- matrix(lf$L[1, ], nobs, J, byrow = T)
            tmat <- matrix(lf$xind,nobs,J,byrow=TRUE)
            if (grepl(paste(cov, "\\.[t]mat", sep = ""),
                      deparse(lf$call$x))) {
              if (length(attr(newdata, "L"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new L matrix of quadrature weights only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "L")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied L matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                    
                  }else{
                    L <- as.vector(attr(newdata, "L"))
                  }
                }
              }
              if (length(attr(newdata, "tmat"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new tmat matrix of observation times only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "tmat")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied tmat matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                  }else{
                    tmat <- as.vector(attr(newdata, "tmat"))
                  }
                }
              }
              if (!is.null(lf$presmooth)) {
                if(type=='lpmatrix' & J!=length(lf$xind)){
                  warning('Presmoothing of new functional covariates is only implemented for when new covariates observed at same time points as original data. No presmoothing of new covariates done.')
                } else if (is.logical(lf$presmooth)) {
                  # lf_old term
                  if (lf$presmooth) {
                    newXfd <- fd(tcrossprod(lf$Xfd$y2cMap,
                                            newdata[[cov]]), lf$Xfd$basis)
                    newdata[[cov]] <- t(eval.fd(lf$xind, newXfd))
                  }
                } else {
                  # lf() term
                  newdata[[cov]] <- lf$prep.func(newX = newdata[[cov]])$processed  
                }
              }
              if(type=='lpmatrix') newdata[[cov]] <- as.vector(newdata[[cov]])
              gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
              gamdata[[paste("L.", cov, sep = "")]] <- L *
                newdata[[cov]]
            }
          }
          if (!(is.af || is.lf)) {
            gamdata[[cov]] <- drop(newdata[[cov]])
          }
        }
      }
    }
    gamdata <- list2df(gamdata)
    call[["newdata"]] <- gamdata
  }
  else {
    call$newdata <- eval(call$newdata)
    nobs <- object$fgam$nobs
  }
  if (PredOutOfRange) {
    suppressMessages(trace(splines::spline.des, at = 2, quote({
      outer.ok <- TRUE
    }), print = FALSE))
    on.exit(suppressMessages(try(untrace(splines::spline.des),
                                 silent = TRUE)))
  }
  oterms <- terms
  if (type %in% c("terms", "iterms") & !all(terms %in% c(names(object$smooth),
                                                         attr(object$pterms, "term.labels")))) {
    if (terms %in% unlist(varmap)) {
      tnames <- c(names(object$smooth), attr(object$pterms,
                                             "term.labels"))
      terms <- tnames[grep(paste(string, "[,\\.\\)]", sep = ""),
                           tnames)]
    }
    else {
      stop("Invalid terms specified")
    }
  }
  call[[1]] <- mgcv::predict.gam
  call$object <- as.name("object")
  call$terms <- terms
  res <- eval(call)
  if (type %in% c("terms", "iterms"))
    colnames(res) <- oterms
  return(res)
}

list2df <- function(l){
  # make a list into a dataframe -- matrices are left as matrices!
  nrows <- sapply(l, function(x) nrow(as.matrix(x)))
  stopifnot(length(unique(nrows)) == 1)
  ret <- data.frame(rep(NA, nrows[1]))
  for(i in 1:length(l)) ret[[i]] <- l[[i]]
  names(ret) <- names(l)
  return(ret)
  
  
  
}



###############################################################################



resam<- function(v,X, Time, Delta, obsT, Y){
  	n=length(X)
	dN=(obsT>0); 
	tn=max(rowSums(dN));
	m=rowSums(dN)

	array=as.vector(obsT)
	sortIndex = sort(array, index.return=TRUE) 
	I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
	T1=array[I]
	J=sapply(1:length(array),function(k){min(which(T1==array[k]))})

	n=nrow(obsT);tn=ncol(obsT)
	Y_=matrix(0,n,length(T1));dN_=Y_;
	for (i in 1:n){
  		tryCatch(
    		{ ID=J[seq(i,(i+n*(m[i]-1)),by=n)];
    		Y_[i,ID]=Y[i,1:m[i]];
    		dN_[i,ID]=1; }
    		, error = function(e){})
	}
	
  dat=data.frame('time'=Time,'status'=Delta,'X'=X)
  
  ## Cox
  fit=coxph(Surv(time, status) ~ X,data = dat,weights = v) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((Time[k]-cumhaz.tt$time)[(Time[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[as.numeric(ind)])
  
  ## Estimate beta
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt))
  dGbar=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt))
  dGbar2=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= Time[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        temp = as.numeric(v)*Y_[,k]*dN_[,k]*psai
        temp2=as.numeric(v)*psai
        dYbar[i,k]=sum(temp)/sum(temp2)
        
        temp = as.numeric(v)*X*dN_[,k]*psai
        dGbar[i,k]=sum(temp)/sum(temp2)
        
        temp = as.numeric(v)*X*psai
        Xbar[i,k]=sum(temp)/sum(temp2)
        
        temp = as.numeric(v)*X*dN_[,k]*psai*T1[k]
        dGbar2[i,k]=sum(temp)/sum(temp2)
        
        temp = as.numeric(v)*X*psai*T1[k]
        Xbar2[i,k]=sum(temp)/sum(temp2)
      }
    }
  }
  
  nu=rep(NA,length(tt))
  nu2=rep(NA,length(tt))
  sigma=matrix(0,2,2)
  for (k in 1:length(tt)){
  	temp = as.numeric(v)*(X-Xbar[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu[k]=sum(na.omit(temp))
    temp = as.numeric(v)*(X*tt[k]-Xbar2[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu2[k]=sum(na.omit(temp))
    tmp1=rbind(as.numeric(v)*as.numeric((X-Xbar[,k])*(Time>=tt[k])),
               as.numeric(v)*as.numeric((X*tt[k]-Xbar2[,k])*(Time>=tt[k]))); tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  
  est=ginv(sigma)%*%c(sum(nu),sum(nu2))   
  
  out=c((eta_hat),est[1],est[2],ginv(sigma[2,2])%*%c(sum(nu2)))
  
}


resam_nonlinear<- function(v,X, Time, Delta, obsT, Y, gap_time){
  #data setup
	n=length(X)
	dN=(obsT>0); 
	tn=max(rowSums(dN));
	m=rowSums(dN)

	array=as.vector(obsT)
	sortIndex = sort(array, index.return=TRUE) 
	I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
	T1=array[I]
	J=sapply(1:length(array),function(k){min(which(T1==array[k]))})

	n=nrow(obsT);tn=ncol(obsT)
	Y_=matrix(0,n,length(T1));dN_=Y_;
	for (i in 1:n){
  		tryCatch(
    		{ ID=J[seq(i,(i+n*(m[i]-1)),by=n)];
    		Y_[i,ID]=Y[i,1:m[i]];
    		dN_[i,ID]=1; }
    		, error = function(e){})
	}

	#coxph	
  dat=data.frame('time'=Time,'status'=Delta,'X'=X)
 
  fit=coxph(Surv(time, status) ~ X,data = dat,weights = v) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((Time[k]-cumhaz.tt$time)[(Time[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[as.numeric(ind)])
  

  ## nonlinear: splines
 end_time = max(Time[Delta == 1])
  nn=length(seq(0,end_time,gap_time))

  phi = data.frame(bs(tt,df=5))
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt));dGbar=matrix(NA,n,length(tt))
  Xbar1=matrix(NA,n,length(tt));dGbar1=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt));dGbar2=matrix(NA,n,length(tt))
  Xbar3=matrix(NA,n,length(tt));dGbar3=matrix(NA,n,length(tt))
  Xbar4=matrix(NA,n,length(tt));dGbar4=matrix(NA,n,length(tt))
  Xbar5=matrix(NA,n,length(tt));dGbar5=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= Time[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(as.numeric(v)*Y_[,k]*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        dGbar[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        Xbar[i,k]=sum(as.numeric(v)*X*psai)/sum(as.numeric(v)*psai)
        dGbar1[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,1])/sum(as.numeric(v)*psai)
        Xbar1[i,k]=sum(as.numeric(v)*X*psai*phi[k,1])/sum(as.numeric(v)*psai)
        dGbar2[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,2])/sum(as.numeric(v)*psai)
        Xbar2[i,k]=sum(as.numeric(v)*X*psai*phi[k,2])/sum(as.numeric(v)*psai)
        dGbar3[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,3])/sum(as.numeric(v)*psai)
        Xbar3[i,k]=sum(as.numeric(v)*X*psai*phi[k,3])/sum(as.numeric(v)*psai)
        dGbar4[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,4])/sum(as.numeric(v)*psai)
        Xbar4[i,k]=sum(as.numeric(v)*X*psai*phi[k,4])/sum(as.numeric(v)*psai)
        dGbar5[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,5])/sum(as.numeric(v)*psai)
        Xbar5[i,k]=sum(as.numeric(v)*X*psai*phi[k,5])/sum(as.numeric(v)*psai)
      }
    }
  }

  nu=rep(NA,length(tt));nu1=rep(NA,length(tt))
  nu2=rep(NA,length(tt));nu3=rep(NA,length(tt))
  nu4=rep(NA,length(tt));nu5=rep(NA,length(tt))
  sigma=matrix(0,6,6)
  for (k in 1:length(tt)){
  	temp = as.numeric(v)*(X-Xbar[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu[k]=sum(na.omit(temp))
    temp = as.numeric(v)*(X*phi[k,1]-Xbar1[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu1[k]=sum(na.omit(temp))
    temp = as.numeric(v)*(X*phi[k,2]-Xbar2[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu2[k]=sum(na.omit(temp))
    temp = as.numeric(v)*(X*phi[k,3]-Xbar3[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu3[k]=sum(na.omit(temp))
    temp = as.numeric(v)*(X*phi[k,4]-Xbar4[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu4[k]=sum(na.omit(temp))
    temp = as.numeric(v)*(X*phi[k,5]-Xbar5[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu5[k]=sum(na.omit(temp))
    tmp1=rbind(as.numeric(v)*(X-Xbar[,k])*(Time>=tt[k]),as.numeric(v)*(X*phi[k,1]-Xbar1[,k])*(Time>=tt[k]),
               as.numeric(v)*(X*phi[k,2]-Xbar2[,k])*(Time>=tt[k]),
               as.numeric(v)*(X*phi[k,3]-Xbar3[,k])*(Time>=tt[k]),as.numeric(v)*(X*phi[k,4]-Xbar4[,k])*(Time>=tt[k]),
               as.numeric(v)*(X*phi[k,5]-Xbar5[,k])*(Time>=tt[k]))
    tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*phi[k,1]*dN_[,k]-dGbar1[,k]),(X*phi[k,2]*dN_[,k]-dGbar2[,k]),
               (X*phi[k,3]*dN_[,k]-dGbar3[,k]),(X*phi[k,4]*dN_[,k]-dGbar4[,k]),(X*phi[k,5]*dN_[,k]-dGbar5[,k]))
    tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  est.n=ginv(sigma)%*%c(sum(nu),sum(nu1),sum(nu2),sum(nu3),sum(nu4),sum(nu5))

  delta.s=(as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1)/tt)
  t.grid <- seq(0, end_time, by = gap_time)
  delta.s.hat <- approx(x = tt, y = delta.s,
                      xout = t.grid, rule = 2)$y
  
   out=c(eta_hat, delta.s.hat)
  
}

#only works for treatment, not set up for actual covariates

sjm_linear_estimate=function(X, Time, Delta, obsT, Y, n.resample=100, var=FALSE){
	#data setup
	n=length(X)
	dN=(obsT>0); 
	tn=max(rowSums(dN));
	m=rowSums(dN)

	array=as.vector(obsT)
	sortIndex = sort(array, index.return=TRUE) 
	I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
	T1=array[I]
	J=sapply(1:length(array),function(k){min(which(T1==array[k]))})

	n=nrow(obsT);tn=ncol(obsT)
	Y_=matrix(0,n,length(T1));dN_=Y_;
	for (i in 1:n){
  		tryCatch(
    		{ ID=J[seq(i,(i+n*(m[i]-1)),by=n)];
    		Y_[i,ID]=Y[i,1:m[i]];
    		dN_[i,ID]=1; }
    		, error = function(e){})
	}
	
  ## Cox
  dat=data.frame('time'=Time,'status'=Delta,'X'=X)
  fit=coxph(Surv(time, status) ~ X,data = dat) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[as.numeric(ind)])
  
  ## Estimate beta
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt))
  dGbar=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt))
  dGbar2=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= Time[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(Y_[,k]*dN_[,k]*psai)/sum(psai)
        dGbar[i,k]=sum(X*dN_[,k]*psai)/sum(psai)
        Xbar[i,k]=sum(X*psai)/sum(psai)
        dGbar2[i,k]=sum(X*dN_[,k]*psai*T1[k])/sum(psai)
        Xbar2[i,k]=sum(X*psai*T1[k])/sum(psai)
      }
    }
  }
  
  nu=rep(NA,length(tt))
  nu2=rep(NA,length(tt))
  sigma=matrix(0,2,2)
  for (k in 1:length(tt)){
  	temp = (X-Xbar[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu[k]= sum(na.omit(temp))
    temp2 = (X*tt[k]-Xbar2[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu2[k]=sum(na.omit(temp2))
    tmp1=rbind(as.numeric((X-Xbar[,k])*(Time>=tt[k])),as.numeric((X*tt[k]-Xbar2[,k])*(Time>=tt[k]))); tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  
  est=c((eta_hat),ginv(sigma)%*%c(sum(nu),sum(nu2)))
  out=list("est" = as.vector(est))
  
  if(var) {
  	#### variance estimation
  	v=matrix(rexp(n*n.resample),nrow=n)
  	temp_list <- lapply(seq_len(ncol(v)), function(j) {
  resam(v[, j], X = X, Time = Time, Delta = Delta, obsT = obsT, Y = Y)
})
	temp <- do.call(cbind, temp_list)
  	est.se=apply(temp,1,sd)[1:3]
  	ci.lower=apply(temp,1,quantile,0.025,na.rm=TRUE)[1:3]
  	ci.upper=apply(temp,1,quantile,0.975,na.rm=TRUE)[1:3]
  	out = list("est" = as.vector(est), "SE" = as.vector(est.se),"CI_lower"= as.vector(ci.lower),"CI_upper"= as.vector(ci.upper))
   }
  return(out)

}

sjm_nl_estimate=function(X, Time, Delta, obsT, Y, gap_time=0.1, n.resample=100, var=FALSE){
	
	#data setup
	n=length(X)
	dN=(obsT>0); 
	tn=max(rowSums(dN));
	m=rowSums(dN)

	array=as.vector(obsT)
	sortIndex = sort(array, index.return=TRUE) 
	I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
	T1=array[I]
	J=sapply(1:length(array),function(k){min(which(T1==array[k]))})

	n=nrow(obsT);tn=ncol(obsT)
	Y_=matrix(0,n,length(T1));dN_=Y_;
	for (i in 1:n){
  		tryCatch(
    		{ ID=J[seq(i,(i+n*(m[i]-1)),by=n)];
    		Y_[i,ID]=Y[i,1:m[i]];
    		dN_[i,ID]=1; }
    		, error = function(e){})
	}

	#coxph	
  dat=data.frame('time'=Time,'status'=Delta,'X'=X)
  fit=coxph(Surv(time, status) ~ X,data = dat) 
  
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[as.numeric(ind)])
  
    
  ## nonlinear: splines
  end_time =  end_time = max(Time[Delta == 1])
  nn=length(seq(0,end_time,gap_time))

  phi = data.frame(bs(tt,df=5))
  dat.spline=data.frame(tt=tt,phi=phi)
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt));dGbar=matrix(NA,n,length(tt))
  Xbar1=matrix(NA,n,length(tt));dGbar1=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt));dGbar2=matrix(NA,n,length(tt))
  Xbar3=matrix(NA,n,length(tt));dGbar3=matrix(NA,n,length(tt))
  Xbar4=matrix(NA,n,length(tt));dGbar4=matrix(NA,n,length(tt))
  Xbar5=matrix(NA,n,length(tt));dGbar5=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= Time[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(Y_[,k]*dN_[,k]*psai)/sum(psai)
        dGbar[i,k]=sum(X*dN_[,k]*psai)/sum(psai)
        Xbar[i,k]=sum(X*psai)/sum(psai)
        dGbar1[i,k]=sum(X*dN_[,k]*psai*phi[k,1])/sum(psai)
        Xbar1[i,k]=sum(X*psai*phi[k,1])/sum(psai)
        dGbar2[i,k]=sum(X*dN_[,k]*psai*phi[k,2])/sum(psai)
        Xbar2[i,k]=sum(X*psai*phi[k,2])/sum(psai)
        dGbar3[i,k]=sum(X*dN_[,k]*psai*phi[k,3])/sum(psai)
        Xbar3[i,k]=sum(X*psai*phi[k,3])/sum(psai)
        dGbar4[i,k]=sum(X*dN_[,k]*psai*phi[k,4])/sum(psai)
        Xbar4[i,k]=sum(X*psai*phi[k,4])/sum(psai)
        dGbar5[i,k]=sum(X*dN_[,k]*psai*phi[k,5])/sum(psai)
        Xbar5[i,k]=sum(X*psai*phi[k,5])/sum(psai)
      }
    }
  }

  nu=rep(NA,length(tt));nu1=rep(NA,length(tt))
  nu2=rep(NA,length(tt));nu3=rep(NA,length(tt))
  nu4=rep(NA,length(tt));nu5=rep(NA,length(tt))
  sigma=matrix(0,6,6)
  for (k in 1:length(tt)){
  	temp=(X-Xbar[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu[k]=sum(na.omit(temp))
    temp = (X*phi[k,1]-Xbar1[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu1[k]=sum(na.omit(temp))
    temp = (X*phi[k,2]-Xbar2[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu2[k]=sum(na.omit(temp))
    temp = (X*phi[k,3]-Xbar3[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu3[k]=sum(na.omit(temp))
    temp = (X*phi[k,4]-Xbar4[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu4[k]=sum(na.omit(temp))
    temp = (X*phi[k,5]-Xbar5[,k])*(Time>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k])
    nu5[k]=sum(na.omit(temp))
    tmp1=rbind((X-Xbar[,k])*(Time>=tt[k]),(X*phi[k,1]-Xbar1[,k])*(Time>=tt[k]),(X*phi[k,2]-Xbar2[,k])*(Time>=tt[k]),
               (X*phi[k,3]-Xbar3[,k])*(Time>=tt[k]),(X*phi[k,4]-Xbar4[,k])*(Time>=tt[k]),(X*phi[k,5]-Xbar5[,k])*(Time>=tt[k]))
    tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*phi[k,1]*dN_[,k]-dGbar1[,k]),(X*phi[k,2]*dN_[,k]-dGbar2[,k]),
               (X*phi[k,3]*dN_[,k]-dGbar3[,k]),(X*phi[k,4]*dN_[,k]-dGbar4[,k]),(X*phi[k,5]*dN_[,k]-dGbar5[,k]))
    tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  est.n=ginv(sigma)%*%c(sum(nu),sum(nu1),sum(nu2),sum(nu3),sum(nu4),sum(nu5))
  
   
  delta.s=(as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1)/tt)
  t.grid <- seq(0, end_time, by = gap_time)
  delta.s.hat <- approx(x = tt, y = delta.s,
                      xout = t.grid, rule = 2)$y
  
  out=list("est" = eta_hat, "est.t" = as.vector(delta.s.hat),"t_grid" = t.grid)
  
	if(var) {	
  		v=matrix(rexp(n*n.resample),nrow=n)
  		temp_list <- lapply(seq_len(ncol(v)), function(j) {
  resam_nonlinear(v[, j], X = X, Time = Time, Delta = Delta, obsT = obsT, Y = Y, gap_time=gap_time)
})
	temp <- do.call(cbind, temp_list)
  		est.se=sd(temp[1,])
  		est.up = quantile(temp[1,], 0.975,na.rm=TRUE)
  		est.low = quantile(temp[1,], 0.025,na.rm=TRUE)
  		
  		delta.s.se=apply(temp,1,sd,na.rm=TRUE)[-1]
  		delta.s.up=apply(temp,1,quantile,0.975,na.rm=TRUE)[-1]
  		delta.s.low=apply(temp,1,quantile,0.025,na.rm=TRUE)[-1]
		out=list("est" = as.vector(eta_hat), "est_t" = as.vector(delta.s.hat), "t_grid" = t.grid, "SE_est"=est.se, "SE_est_t" = as.vector(delta.s.se), "CI_lower_est" = est.low, "CI_upper_est" = est.up, "CI_lower_est_t"= as.vector(delta.s.low),"CI_upper_est_t"= as.vector(delta.s.up))
	}
	return(out)
}
