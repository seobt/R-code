################################################################################
### DSE: Non-paraller version
#################################################################################
DSE = function(x,Y,m=5,min_lambda=0.05,max_lambda=0.85){
  n=length(Y) # Sample size
  nbeta=dim(x)[2] # number of covariates
  ###########################################
  ####### fit lse
  fitlse = lm(Y ~ x)
  
  sigma = sig = summary(fitlse)$sigma
  stepsize = 0.05 * sig
  ## Set lambda range
  lambda = rev(seq(min_lambda * sig, max_lambda*sig, by = stepsize))
  lambdaveclen = length(lambda)
  con = gradient = max.grad = rep(0, lambdaveclen)
  bic = numpara = ll = scalepara = mse = rep(Inf, lambdaveclen)
  out = vector("list", lambdaveclen)
  ## Set initial values
  mix = dden(0, 1)
  betas=fitlse$coef  

  ######## Construct pseudo DS data ####################
  #m = 5
  if (m==1) {
    xa=x
    repY=Y
    pseudo2=0
  }
  else {
  id = rep(1:n, m)
  sp = 1 / (m + 1)
  pseudo = qnorm(seq(sp, 1 - sp, sp))
  pseudo2 = pseudo / sqrt(sum(pseudo^2) / m)
  repY = rep(Y, m)
  xa = apply(as.matrix(x), 2, rep, m)}

  ########################################################
  # Find initial value of betas and mix that converges
  test=max(lambda) 
  for (jj in 1:20) {
    u = pseudo2 * test
    xxa=cbind(repY+rep(u, each = n),xa)
    qqq=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = test*1.5)
    if (qqq$con==0) break else test=test*1.1
  }
  mix=qqq$mix
  betas=qqq$beta
  ####################################################################
  ####################################################################
  zxc=NULL
  for (ii in 1:lambdaveclen) {

    u = pseudo2 * lambda[ii]
    xxa=cbind(repY+rep(u, each = n),xa)
    out=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda[ii])
    mix              = out$mix 
    betas            = out$beta
    colnames(betas)=NULL
    con          = out$convergence
    gradient     = max(abs(out$grad))
    max.grad     = out$max.gradient
    numpara      = 2 * length(mix$pt) - 1 + length(betas)
    if (con==3) ll=-Inf else ll=out$ll
    bic= -2 * ll/m + log(n) * numpara
    pt=rep(0,n);pr=rep(0,n);
    pt[1:length(mix$pt)]=mix$pt
    pr[1:length(mix$pr)]=mix$pr
    aaa=c(con,gradient,max.grad,numpara,ll,bic,out$beta,pt,pr)
    zxc=rbind(zxc,aaa)

  }
  endpoint=7+length(betas)
  con=zxc[1:lambdaveclen,1];gradient=zxc[1:lambdaveclen,2];max.grad=zxc[1:lambdaveclen,3];numpara=zxc[1:lambdaveclen,4];ll=zxc[1:lambdaveclen,5];bic=zxc[1:lambdaveclen,6]
  outbeta=zxc[1:lambdaveclen,c(7:(endpoint-1))]
  pt=zxc[,c(endpoint:(endpoint+n-1))];pr=zxc[,c((endpoint+n):(endpoint+2*n-1))]
  ######## Check convergence and rerun with different initial values 
  infind=which(bic==Inf)
  for (nf in c(infind)){

    betas=outbeta[nf-1,]
    nonzeroweight=which(pr[nf-1,]!=0)
    mix=dden(pt[nf-1,nonzeroweight],pr[nf-1,nonzeroweight])
    u = pseudo2 * lambda[nf]
    xxa=cbind(repY+rep(u, each = n),xa)
    out=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda[nf])
    if (out$convergence==3){
      betas=outlse
      mix = dden(0, 1)
      out=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda[nf])
    }
    pt[nf,1:length(out$mix$pt)]=out$mix$pt
    pr[nf,1:length(out$mix$pr)]=out$mix$pr
    numpara[nf]=2*length(out$mix$pt)-1+length(out$betas)
    ll[nf]=out$ll
    con[nf]=out$convergence
    print(out$convergence)
    outbeta[nf,]=out$beta
    bic[nf]=-2 * out$ll/m + log(n) *  (2 * length(out$mix$pt) - 1 + length(out$beta))
  }
  #### Only use convergent result #############################################
  converged=which(con==0)
  lambda_ic=list(lambda=lambda[converged], bic = bic[converged], numpara = numpara[converged], ll = ll[converged])
  indexb=which.min(lambda_ic$bic)
  outsyb=outbeta[indexb,]
  if (indexb==length(lambda_ic$bic) | indexb==1) print('Problem with lambda range')
  support=pt[indexb,]
  weight=pr[indexb,]
  support=support[-which(weight==0)];weight=weight[-which(weight==0)]
  list(beta=outsyb,bic=lambda_ic$bic[indexb],lambda=lambda_ic$lambda[indexb],loglik=lambda_ic$ll[indexb],mix=dden(support,weight))
} 
################################################################################
### DSE: parallel version
#################################################################################
DSE_parallel = function(x,Y,m=5,min_lambda=0.05,max_lambda=0.85){
  n=length(Y) # Sample size
  nbeta=dim(x)[2] # number of covariates
  ###########################################
  ####### fit lse, huber, tukey
  fitlse = lm(Y ~ x) # Use LSE as initial values
  #fitlse=rlm(Y~x, psi=psi.bisquare,method ="M",maxit = 500)

  
  sigma = sig = summary(fitlse)$sigma
  stepsize = 0.05 * sig
  ## Set lambda range
  lambda = rev(seq(min_lambda * sig, max_lambda*sig, by = stepsize))
  lambdaveclen = length(lambda)
  con = gradient = max.grad = rep(0, lambdaveclen)
  bic = numpara = ll = scalepara = mse = rep(Inf, lambdaveclen)
  out = vector("list", lambdaveclen)
  ## Set initial values
  mix = dden(0, 1)
  betas=fitlse$coef  
  ######## Construct pseudo DS data ####################
  #m = 5
  if (m==1) {
    xa=x
    repY=Y
    pseudo2=0
  }
  else {
  id = rep(1:n, m)
  sp = 1 / (m + 1)
  pseudo = qnorm(seq(sp, 1 - sp, sp))
  pseudo2 = pseudo / sqrt(sum(pseudo^2) / m)
  repY = rep(Y, m)
  xa = apply(as.matrix(x), 2, rep, m)
  }
  ########################################################
  # Find initial value of betas and mix that converges
  test=max(lambda) 
  for (jj in 1:50) {
    u = pseudo2 * test
    xxa=cbind(repY+rep(u, each = n),xa)
    qqq=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = test*1.5)
    if (qqq$con==0) break else test=test*1.1
  }
  mix=qqq$mix
  betas=qqq$beta

  ####################################################################
  ####################################################################
  zxc=foreach(ii=1:lambdaveclen, .combine='rbind',.packages='nnls') %dopar% {
    #for (ii in 1:lambdaveclen) {
    source("./symmetric.R")
    #source("./simstuff.R")
    u = pseudo2 * lambda[ii]
    xxa=cbind(repY+rep(u, each = n),xa)
    out=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda[ii])
    mix              = out$mix 
    betas            = out$beta
    con          = out$convergence
    gradient     = max(abs(out$grad))
    max.grad     = out$max.gradient
    numpara      = 2 * length(mix$pt) - 1 + length(betas)
    if (con==3) ll=-Inf else ll=out$ll
    bic= -2 * ll/m + log(n) * numpara
    pt=rep(0,n);pr=rep(0,n);
    pt[1:length(mix$pt)]=mix$pt
    pr[1:length(mix$pr)]=mix$pr
    c(con,gradient,max.grad,numpara,ll,bic,out$beta,pt,pr)
    #aaa=c(con,gradient,max.grad,numpara,ll,bic,out$beta,pt,pr)
    #zxc=rbind(zxc,aaa)
  }
  rownames(zxc)=NULL
  endpoint=7+length(betas)
  con=zxc[1:lambdaveclen,1];gradient=zxc[1:lambdaveclen,2];max.grad=zxc[1:lambdaveclen,3];numpara=zxc[1:lambdaveclen,4];ll=zxc[1:lambdaveclen,5];bic=zxc[1:lambdaveclen,6]
  outbeta=zxc[1:lambdaveclen,c(7:(endpoint-1))]
  pt=zxc[,c(endpoint:(endpoint+n-1))];pr=zxc[,c((endpoint+n):(endpoint+2*n-1))]
  ######## Check convergence and rerun with different initial values 
  infind=which(bic==Inf)
  for (nf in c(infind)){
    betas=outbeta[nf-1,]
    nonzeroweight=which(pr[nf-1,]!=0)
    mix=dden(pt[nf-1,nonzeroweight],pr[nf-1,nonzeroweight])
    u = pseudo2 * lambda[nf]
    xxa=cbind(repY+rep(u, each = n),xa)
    out=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda[nf])
    if (out$convergence==3){
      betas=outlse
      mix = dden(0, 1)
      out=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda[nf])
    }
    pt[nf,1:length(out$mix$pt)]=out$mix$pt
    pr[nf,1:length(out$mix$pr)]=out$mix$pr
    numpara[nf]=2*length(out$mix$pt)-1+length(out$betas)
    ll[nf]=out$ll
    con[nf]=out$convergence
    print(out$convergence)
    outbeta[nf,]=out$beta
    bic[nf]=-2 * out$ll/m + log(n) *  (2 * length(out$mix$pt) - 1 + length(out$beta))
  }
  #### Only use convergent result #############################################
  converged=which(con==0)
  lambda_ic=list(lambda=lambda[converged], bic = bic[converged], numpara = numpara[converged], ll = ll[converged])
  indexb=which.min(lambda_ic$bic)
  outsyb=outbeta[indexb,]
  if (indexb==length(lambda_ic$bic) | indexb==1) print('Problem with lambda range')
  support=pt[indexb,]
  weight=pr[indexb,]
  support=support[-which(weight==0)];weight=weight[-which(weight==0)]
  list(beta=outsyb,bic=lambda_ic$bic[indexb],lambda=lambda_ic$lambda[indexb],loglik=lambda_ic$ll[indexb],mix=dden(support,weight))
} 
#####################################################
#####################################################
################################################################################
### DSE: parallel version with fixed bandwidth
#################################################################################
DSE_h= function(x,Y,m=5,lambda,mix,betas=fitlse$coef){
  n=length(Y) # Sample size
  nbeta=dim(x)[2] # number of covariates
  ###########################################
  ####### fit lse, huber, tukey
  fitlse = lm(Y ~ x)
  
  sigma = sig = summary(fitlse)$sigma
  if (m==1) {
    xa=x
    repY=Y
    pseudo2=0
  } else{ 
    id = rep(1:n, m)
    sp = 1 / (m + 1)
    pseudo = qnorm(seq(sp, 1 - sp, sp))
    pseudo2 = pseudo / sqrt(sum(pseudo^2) / m)
    repY = rep(Y, m)
    xa = apply(as.matrix(x), 2, rep, m)
  }
  ########################################################
  u = pseudo2 * lambda
  xxa=cbind(repY+rep(u, each = n),xa)
  return(cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda))
} 


############################################################################
############################################################################
################################################################################
### DSE: parallel version with fixed bandwidth
#################################################################################
DSE_h2= function(x,Y,m=5,lambda,mix,betas=fitlse$coef){
  n=length(Y) # Sample size
  nbeta=dim(x)[2] # number of covariates
  ###########################################
  ####### fit lse, huber, tukey
  fitlse = lm(Y ~ x)
  
  sigma = sig = summary(fitlse)$sigma
  ######## Construct pseudo DS data ####################
  #m = 5
  if (m==1) {
    xa=x
    repY=Y
    pseudo2=0
  } else{ 
    id = rep(1:n, m)
    sp = 1 / (m + 1)
    pseudo = qnorm(seq(sp, 1 - sp, sp))
    pseudo2 = pseudo / sqrt(sum(pseudo^2) / m)
    repY = rep(Y, m)
    xa = apply(as.matrix(x), 2, rep, m)
  }
  ########################################################
  u = pseudo2 * lambda
  xxa=cbind(repY+rep(u, each = n),xa)
  #ini_result=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = 2)
  result=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda)
  if (result$con>0) {
    test=lambda 
    for (jj in 1:20) {
      u = pseudo2 * test
      xxa2=cbind(repY+rep(u, each = n),xa)
      qqq2=cnmms(xxa2, init = list(mix = mix, beta = betas), lambda = test*1.5)
      if (qqq2$con==0) break else test=test*1.1
    }
    mix=qqq2$mix
    betas=qqq2$beta
    result=cnmms(xxa, init = list(mix = mix, beta = betas), lambda = lambda)
  }

  return(result)
  
  #mix=qqq$mix
  #betas=qqq$beta
  ####################################################################
  #list(beta=qqq$beta,bic=lambda_ic$bic[indexb],lambda=lambda_ic$lambda[indexb],loglik=lambda_ic$ll[indexb],mix=dden(support,weight))
} 










