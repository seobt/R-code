library(lsei)
################################################################################
initial = function(x, lambda, beta, mix, kmax, start) UseMethod("initial")
logd = function(x, lambda, beta, pt, which) UseMethod("logd")
valid = function(x, beta, mix) UseMethod("valid")

cnmms = function(x, lambda, init = NULL, plot = F, start = 1, which = c(1, 1, 1), kmax = 10, maxit = 1000, tol = 1e-10, grid = 200, llt = NULL) {
  x = sym(x)
  if(is.null(init) || is.null(init$mix)) init = initial.snpmle(x, lambda, init, kmax, start) else init = init
  beta = init$beta
  mix = init$mix
  ll1 = -Inf
  convergence = 1
  for(i in 1:maxit) {
    l = logd(x, lambda, beta, mix$pt, which = c(1, 0, 0, 0))$ld
    ma = apply(l, 1, max)
    dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100 
    gridpoints = seq(0, max(abs(range(x, beta))), length = grid)  
    g = maxgrad(x, lambda, beta, dmix, ma, grid = gridpoints, tol = -Inf)
    gradient = max(g$grad)
    mix = dden(c(mix$pt, g$pt), c(mix$pr, rep(0, length(g$pt))))
    lpt = logd(x, lambda, beta, mix$pt, which = c(1, 0, 0, 0))$ld
    dpt = pmin(exp(lpt - ma), 1e100)
    a = cbind(dpt / dmix - 2)
    r = nnls(rbind(a, rep(1, length(mix$pt))), c(rep(0, nrow(a)), 1))
    sol = r$x / sum(r$x)
    r = lsch(mix, beta, dden(mix$pt, sol), beta, x, lambda, which = c(1, 0, 0))
    mix = collapse.snpmle(r$mix, beta, x, lambda)
    r = bfgs(mix, beta, x, lambda, which = which)
    mix = r$mix
    beta = r$beta
    if(r$conv == 3) {
      l = logd(x, lambda, beta, mix$pt, which = c(1, 0, 0, 0))$ld 
      ma = apply(l, 1, max)
      dmix = drop(exp(l - ma) %*% mix$pr) + 1e-100
      grad = grad(x, lambda, mix$pt, beta, dmix, ma, order = 0)$d0
      gradient = max(grad)
      if(length(gradient) == 1 && isTRUE(all.equal(gradient, 0))) r$conv = 0
      convergence = r$conv
      break
    }       
    if(is.null(llt)) {if(r$ll >= ll1 && r$ll <= ll1 + tol) {convergence = 0; break}}
    else {if(r$ll >= llt) {convergence = 0; break}
          else if(r$ll >= ll1 && r$ll <= ll1 + 1e-16) {convergence = 1; break}}
    ll1 = r$ll
  }
  if(isTRUE(all.equal(mix$pt[1], .Machine$double.eps))) mix$pt[1] = 0   
  if(sum((eps.range <- range(x, beta)) > 0) != 1) print("all residuals are either positive or negative")  
  if(plot) plotgrad(x, lambda, mix, beta, mix)  
  if(anyDuplicated(round(mix$pt, 6))) mix = unique(mix, prec = c(min(diff(sort(mix$pt))), 0))
  if(anyDuplicated(round(mix$pt, 6))) {print("sym - may collapse"); print(mix)}
  list(mix = mix, lambda = lambda, scale = sqrt(sum((mix$pt^2 + lambda^2) * mix$pr)), beta = beta, eps.range = eps.range, num.iterations = i, ll = r$ll, grad = r$grad, max.gradient = gradient, convergence = convergence)
}
################################################################################
maxgrad = function(x, lambda, beta, dmix, ma, grid = 100, tol = -Inf, maxit = 100) {
  np = length(grid)
  dg = grad(x, lambda, grid, beta, dmix, ma, order = 1)$d1
  jmax = (1:(np - 1))[dg[1:(np - 1)] > 0 & dg[2:np] < 0]
  if(length(jmax) < 1) return
  pt = (grid[jmax] + grid[jmax + 1]) * .5
  left = grid[jmax]
  right = grid[jmax + 1]
  pt.old = left
  if(length(pt) != 0) {
    for(i in 1:maxit) {
      step = pt - pt.old
      d1.old = grad(x, lambda, pt.old, beta, dmix, ma, order = 1)$d1
      pt.old = pt
      d1 = grad(x, lambda, pt, beta, dmix, ma, order = 1)$d1
      left[d1 > 0] = pt[d1 > 0]
      right[d1 < 0] = pt[d1 < 0]
      pt = pt - (d1 * step)/(d1 - d1.old)
      j = is.na(pt) | pt < left | pt > right
      pt[j] = (left[j] + right[j]) * .5
      if( max(abs(pt - pt.old)) <= 1e-14 * diff(range(grid)) ) break
    }
  }
  else i = 0
  if(dg[np] >= 0) pt = c(pt, grid[np])
  if(dg[1] <= 0) pt = c(grid[1], pt)
  if(length(pt) == 0) stop("no new support point found")
  g = grad(x, lambda, pt, beta, dmix, ma, order = 0)$d0
  names(pt) = names(g) = NULL
  j = g >= tol
  list(pt = pt[j], grad = g[j], num.iterations = i)
}

grad = function(x, lambda, pt, beta, dmix, ma, order = 0) {
  if(is.dden(dmix)) {
    l = logd(x, lambda, beta, dmix$pt, which=c(1, 0, 0))$ld
    ma = apply(l, 1, max)
    dmix = drop(exp(l - ma) %*% dmix$pr) + 1e-100
  }
  g = vector("list", length(order))
  names(g) = paste("d", order, sep="")
  which = c(1, 0, 0)
  if(order == 1) which[3] = 1
  dl = logd(x, lambda, beta, pt, which = which)
  dpt = pmin(exp(dl$ld - ma), 1e100)
  if(order == 0) g$d0 = colSums(dpt / dmix) - length(x)
  if(order == 1) g$d1 = colSums(dpt * dl$dt1 / dmix)
  g
}

plotgrad = function(x, lambda, mix, beta, dmix, len = 1000, xlab = expression(theta), ylab = expression(d(theta * "; " * G)), cex = 1, lower, upper, ...) {
  if( missing(lower) || missing(upper) ) {
    endpt = max(abs(range(x, beta)))
    if( missing(lower) ) lower = - endpt * 1.05
    if( missing(upper) ) upper = endpt * 1.05
  }
  pt = seq(lower, upper, len = len)
  g = grad(x, lambda, pt, beta, dmix)$d0
  plot(pt, g, type = "l", col = "blue", xlab = xlab, ylab = ylab, cex = cex, cex.axis = cex, cex.lab = cex, ... )
  if(is.dden(mix)) {
    j = mix$pr != 0
    points(c(-mix$pt[j], mix$pt[j]), rep(0, 2 * length(mix$pt[j])), pch = 20, col = "red")     
  }
  lines(c(lower, upper), c(0, 0), col = "black")
}
################################################################################
bfgs = function(mix, beta, x, lambda, tol = 1e-16, maxit = 100, which = c(1, 1, 1), D = NULL) {
  k1 = if(which[1]) length(mix$pr) - 1 else 0
  k2 = if(which[2]) length(mix$pt) else 0
  k3 = if(which[3]) length(beta) else 0
  if(k1 == 0) which[1] = 0
  if(sum(which) == 0) stop("No parameter specified to be updated in bfgs()")
  if(is.null(D)) D = diag(-1, nrow = k1 + k2 + k3)
  else if(nrow(D) != k1 + k2 + k3) stop("Provided D has unmatching dimensions")
  dl = dll.snpmle(x, lambda, mix, beta, which = c(1, which))
  ll = dl$ll
  grad = c(if(which[1]) dl$dp[-(k1 + 1)] - dl$dp[k1 + 1] else NULL, if(which[2]) dl$dt else NULL, if(which[3]) dl$db else NULL)
  r = list(mix = mix, beta = beta, ll = ll, grad = grad, convergence = 1)
  prmt = c(if(which[1]) r$mix$pr[-(k1+1)] else NULL, if(which[2]) r$mix$pt else NULL, if(which[3]) r$beta else NULL)  
  for(i in 1:maxit) {
    old.r = r
    old.prmt = prmt
    prmt2 = drop(prmt - D %*% r$grad)
    d1 = prmt2 - prmt
    alpha = 1
    if(which[1]) {
      pr = prmt2[1:k1]
      pr2 = c(pr, 1 - sum(pr))
      if(any(pr2 < 0)) {
        if(any(r$mix$pr == 0)) {
          j = r$mix$pr == 0
          mix = dden(r$mix$pt[!j], r$mix$pr[!j])
          r = bfgs(mix, r$beta, x, lambda, tol, maxit, which)
          return(r)
        }
        step = pr2 - r$mix$pr
        ratio = pmax(-pr2, 0) / abs(step)
        jmax = which.max(ratio)
        alpha = 1 - ratio[jmax]
        pr2 = r$mix$pr + alpha * step
        pr2[jmax] = 0
      }
    }
    else pr2 = r$mix$pr
    prmt2 = prmt + alpha * d1
    if(which[2]) pt2 = prmt2[(k1+1):(k1+k2)] else pt2 = r$mix$pt
    if(which[3]) beta2 = prmt2[(k1+k2+1):(k1+k2+k3)] else beta2 = r$beta
    mix2 = dden(pt2, pr2) 
    r = lsch(r$mix, r$beta, mix2, beta2, x, lambda, which = which, brkt = ifelse(alpha == 1, TRUE, FALSE))
    if(r$conv != 0) break
    if(r$ll >= old.r$ll && r$ll <= old.r$ll + tol) {convergence = 0; break}
    prmt = c(if(which[1]) r$mix$pr[-(k1 + 1)] else NULL, if(which[2]) r$mix$pt else NULL, if(which[3]) r$beta else NULL)
    d = prmt - old.prmt
    g = r$grad - old.r$grad
    dg = sum(d * g)
    if(dg < 0) D = D + (1 + drop(t(g) %*% D %*% g) / dg) * outer(d, d) / dg - (outer(d, g) %*% D + D %*% outer(g, d)) / dg
  }
  r$num.iterations = i
  r
}
################################################################################
lsch = function(mix1, beta1, mix2, beta2, x, lambda, maxit = 100, which = c(1, 1, 1), brkt = FALSE) {
  k = length(mix1$pt)
  convergence = 1
  d1 = c(if(which[1]) mix2$pr[-k] - mix1$pr[-k] else NULL, if(which[2]) mix2$pt - mix1$pt else NULL, if(which[3]) beta2 - beta1 else NULL)
  d1.norm = sqrt(sum(d1 * d1))
  s = d1 / d1.norm
  dl1 = dll.snpmle(x, lambda, mix1, beta1, which = c(1, which))
  ll1 = dl1$ll
  grad1 = c(if(which[1]) dl1$dp[-k] - dl1$dp[k] else NULL, if(which[2]) dl1$dt else NULL, if(which[3]) dl1$db else NULL)
  g1s = sum(grad1 * s)
  g1d1 = sum(grad1 * d1)
  if(d1.norm == 0 || g1s <= 0) return(list(mix = mix1, beta = beta1, grad = grad1, ll = ll1, convergence = 3))
  a = 0
  b = 1
  if(which[1] && any(mix2$pr == 0)) brkt = FALSE
  if(!which[3]) brkt = FALSE
  for(i in 1:maxit) {
    repeat {
      m = dden((1 - b) * mix1$pt + b * mix2$pt, (1 - b) * mix1$pr + b * mix2$pr)
      beta = (1 - b) * beta1 + b * beta2
      if(valid.snpmle(x, beta, m)) break
      brkt = FALSE
      b = a + 0.9 * (b - a)
      if(b <= .Machine$double.xmin) stop("valid failed")
    }
    dl = dll.snpmle(x, lambda, m, beta, which = c(1, which))
    ll = dl$ll
    grad = c(if(which[1]) dl$dp[-k] - dl$dp[k] else NULL, if(which[2]) dl$dt else NULL, if(which[3]) dl$db else NULL)
    gs = sum(grad * s)
    if(brkt && gs > g1s * .5 && ll >= ll1 + g1d1 * b * .33) {a = b; b = 2 * b} else break
  }
  if(i == maxit) brkt = FALSE
  alpha = b
  for(i in 1:maxit) {
    g1d = g1d1 * alpha
    if(ll >= ll1 && ll + g1d <= ll) {convergence = 2; break}
    if(brkt) {
      if(ll >= ll1 + g1d * .33 && abs(gs) <= g1s * .5) {convergence = 0; break}
      if(ll >= ll1 + g1d * .33 && gs > g1s * .5) a = alpha else b = alpha
    }
    else {if(ll >= ll1 + g1d * .33) {convergence = 0; break} else b = alpha}
    alpha = (a + b) * .5
    m = dden((1 - alpha) * mix1$pt + alpha * mix2$pt, (1 - alpha) * mix1$pr + alpha * mix2$pr)
    beta = (1 - alpha) * beta1 + alpha * beta2
    dl = dll.snpmle(x, lambda, m, beta, which = c(1, which))
    ll = dl$ll
    grad = c(if(which[1]) dl$dp[-k] - dl$dp[k] else NULL, if(which[2]) dl$dt else NULL, if(which[3]) dl$db else NULL)
    gs = sum(grad * s)
  }
  mix = dden((1 - alpha) * mix1$pt + alpha * mix2$pt, (1 - alpha) * mix1$pr + alpha * mix2$pr)
  beta = (1 - alpha) * beta1 + alpha * beta2
  list(mix = mix, beta = beta, grad = grad, ll = ll, convergence = convergence, num.iterations = i)
}
################################################################################
logLik.snpmle = function(x, lambda, beta, mix) {
  ld = sweep(logd(x, lambda, beta, mix$pt, which = c(1, 0, 0))$ld, 2, log(mix$pr), "+")
  ma = apply(ld, 1, max)
  pid = exp(ld - ma)
  pis = rowSums(pid)
  ll = sum(log(pis) + ma)
  ll
}

dll.snpmle = function(x, lambda, mix, beta, which = c(1, 0, 0, 0)) {
  r = list()
  dl = logd(x, lambda, beta, mix$pt, which = c(1, which[4:3]))                  
  lpt = dl$ld                                                                   
  ma = apply(lpt, 1, max)
  explptma = exp(lpt - ma)
  dpt = pmin(explptma, 1e100)
  dmix = drop(explptma %*% mix$pr) + 1e-100
  dp = dpt / dmix
  p = sweep(dp, 2, mix$pr, "*")
  if(which[1] == 1) r$ll = sum(log(rowSums(explptma %*% mix$pr)) + ma)
  if(sum(which[2:4]) == 0) return(r)
  if(which[2] == 1) r$dp = colSums(dp)
  if(sum(which[3:4]) == 0) return(r)
  if(which[3] == 1) r$dt = colSums(p * dl$dt1)
  if(which[4] == 0) return(r)
  r$db = colSums(apply(sweep(dl$db1, c(1, 2), p, "*"), c(1, 3), sum))
  r
}

initial.snpmle = function(x, lambda, init = NULL, kmax = NULL, start = 1) initial(x, lambda, init$beta, init$mix, kmax, start)

valid.snpmle = function(x, beta, mix) valid(x, beta, mix)

collapse.snpmle = function(mix, beta, x, lambda) {
  mix = sort.dden(mix)
  j = mix$pr == 0
  if(any(j)) mix = dden(mix$pt[!j], mix$pr[!j])
  logLikc = logLik.snpmle(x, lambda, beta, mix)
  if( any(mix$pr <= 1e-3) ) {
    j = mix$pr > 1e-3
    mixt = mix
    mixt$pt = mixt$pt[j]
    mixt$pr = mixt$pr[j] / sum(mix$pr[j])
    logLikt = logLik.snpmle(x, lambda, beta, mixt)
    if(logLikt + 1e-10 >= logLikc) mix = mixt                                                         
  }
  repeat {
    if(length(mix$pt) == 1) break
    prec = min(diff(mix$pt))
    mixt = unique(mix, prec = c(prec, 0))
    logLikt = logLik.snpmle(x, lambda, beta, mixt)      
    if(logLikt + 1e-10 >= logLikc) mix = mixt else break                                 
  }
  mix
}
################################################################################
sym = function(x) {
  xymat = cbind(x[, 1], rep(1, nrow(x)), x[, -1])
  dimnames(xymat) = list(NULL, c("y", paste("x", 1:(ncol(xymat) - 1), sep = "")))
  class(xymat) = "sym"
  xymat
}

length.sym = function(x) nrow(x)

initial.sym = function(x, lambda, beta = NULL, mix = NULL, kmax = NULL, start = 1) {
  if(is.null(beta)) {
    if(ncol(x) == 2) beta = mean(x[, 1]) else beta = lm(x[, 1, drop = FALSE] ~ x[, -(1:2), drop = FALSE])$coef
    names(beta) = NULL
  } 
  if(is.null(kmax)) kmax = 10
  if(is.null(mix)) {
    xij = x[, -1, drop = FALSE]
    xbeta = drop(xij %*% beta) 
    epsilon = x[, 1, drop = TRUE] - xbeta 
    mixa = dden(seq(0, max(abs(range(epsilon))), length = kmax), mix$pr)
    mixb = sort(dden(unique(abs(quantile(epsilon[epsilon < quantile(epsilon, .5)], probs = seq(0, 1, length = kmax), type = 1)))))
    mixc = sort(dden(unique(abs(quantile(epsilon[epsilon > quantile(epsilon, .5)], probs = seq(0, 1, length = kmax), type = 1)))))   
    lla = logLik.snpmle(x, lambda, beta, mixa)
    llb = logLik.snpmle(x, lambda, beta, mixb)
    llc = logLik.snpmle(x, lambda, beta, mixc)
    if(start == 1) mix = switch(which.max(c(lla, llb, llc)), mixa, mixb, mixc) else mix = switch(which.min(c(lla, llb, llc)), mixa, mixb, mixc)   
  }                        
  list(beta = beta, mix = mix)
}

logd.sym = function(x, lambda, beta, pt, which = c(1, 0, 0)) {
  dl = vector("list", 3)
  names(dl) = c("ld", "db1", "dt1")
  xij = x[, -1, drop = FALSE]
  xbeta = drop(xij %*% beta)
  epsilon = x[, 1, drop = TRUE] - xbeta
  outerf1 = outer(epsilon, pt, "-") 
  outerf2 = outer(epsilon, pt, "+") 
  f1 = dnorm(outerf1, 0, lambda)
  f2 = dnorm(outerf2, 0, lambda) 
  sumf1f2 = f1 + f2
  sumf1f2[sumf1f2 < 1e-100] = 1e-100
  denom = lambda^2 * sumf1f2
  c1 = outerf1 * f1 / denom 
  c2 = outerf2 * f2 / denom
  if(which[1] == 1) dl$ld = log(.5 * sumf1f2)
  if(which[2] == 1) dl$db1 = sweep(array(c1 + c2, dim = c(length(x), length(pt), length(beta))), c(1, 3), xij, "*")
  if(which[3] == 1) dl$dt1 = c1 - c2
  dl
}
         
range.sym = function(x, beta, ...) {
  if(ncol(x) == 2) return(range(x[, 1, drop = TRUE]))
  xij = x[, -1, drop = FALSE]
  xbeta = drop(xij %*% beta)
  epsilon = x[, 1, drop = TRUE] - xbeta     
  range(epsilon)
}
  
valid.sym = function(x, beta, mix) all(mix$pt >= 0) && all(mix$pr >= 0)

depsilon = function(x, lambda, pt) .5 * (dnorm(x, pt, lambda) + dnorm(x, -pt, lambda))
dmixepsilon = function(x, lambda, mix, log = FALSE) {
  out = rep(0, length(x))  
  for (i in 1:length(mix$pt)) out = out + mix$pr[i] * depsilon(x, lambda, mix$pt[i])
  if(log) out = log(out)
  drop(out)
}
################################################################################
################################################################################
dden = function(pt, pr=1) {
  if (is.null(pt)) d = list(pt=NULL, pr=NULL)
  else {
    k = max(length(pt), length(pr), na.rm=TRUE)
    pt = rep(pt, len=k)
    if(is.null(pr)) pr = rep(1, length(pt)) else pr = rep(pr, len=k)
    d = list(pt=pt, pr=pr/sum(pr))
  }
  class(d) = "dden"
  d
}

is.dden = function (d) any(class(d) == "dden")

print.dden = function (d, ...) {
  if (is.null(d)) b = matrix(nrow=0, ncol=2)
  else b = cbind(d$pt, d$pr)
  dimnames(b) = list(NULL, c("pt", "pr"))
  print(b, ...)
}

sort.dden = function(d) {
  if( is.null(d) ) return(d)
  index = order(d$pt)
  d$pt = d$pt[index]
  d$pr = d$pr[index]
  d
}

is.unsorted.dden = function(d) is.unsorted(d$pt)

unique.dden = function(d, prec=0) {
  if( length(d$pt) == 1 ) return(d)
  if( is.unsorted.dden(d) ) d = sort.dden(d)
  prec = rep(prec, len=2)
  if ( all(prec < 0) ) return(d)
  pt2 = pt = d$pt
  pr2 = pr = d$pr
  j  = pr <= prec[2]
  pt = pt[!j]
  pr = pr[!j]
  index = 0
  repeat {
    if( length(pt) == 0 ) break
    j = abs(pt[1] - pt) <=  prec[1]
    index = index + 1
    pt2[index] = weighted.mean(pt[j], pr[j])
    pr2[index] = sum( pr[j] )
    pt = pt[!j]
    pr = pr[!j]
  }
  dden(pt=pt2[1:index], pr=pr2[1:index])
}
########################################
################################################################################
int <- function(f, x, tol=1e-10) {
  n <- length(x)
  if(n == 1) return(0)
  simp <- (sum(range(x[-1] - x[-n]) * c(-1, 1)) < tol)
  if(!simp) o <- (1 / 2) * sum((f[-1] + f[-n]) * (x[-1] - x[-n]))
  else o <- ((x[2] - x[1]) / 3) * (sum(f) + sum(f[2:(n - 1)]) + 2 * sum(f[seq(2, n - 1, 2)]))
  return(o)
}
################################################################################
pr.reg <- function(X, d, U, f0, w, perm) {
  n <- length(X)
  if(missing(f0)) {f0 <- 1 + 0 * U; f0 <- f0 / int(f0, U)}
  if(missing(w)) w <- function(i) 1 / (i + 1)
  N <- ncol(perm)
  f.avg <- 0 * U
  L.avg <- 0
  wt <- 0 * X
  for(j in 1:N) {
    f <- f0
    L <- 0
    for(i in 1:n) {
      num <- d(X[perm[i, j]], U) * f
      den <- int(num, U)
      wt[perm[i, j]] <- wt[perm[i, j]] + int(num / U**2, U) / den
      L <- L + log(den)
      f <- (1 - w(i)) * f + w(i) * num / den
    }
    f.avg <- f.avg + f / N
    L.avg <- L.avg + L / N
  }
  return(list(f=f.avg, L=-L.avg, wt=wt / N))
}
################################################################################
################################################################################
prem.reg <- function(X, Y, perm, get.ci=FALSE) {
  tol <- 1e-05
  ols <- lm(Y ~ X - 1)
  sigma <- summary(ols)$sigma
  UUmax <- max(5 * sigma, 50)
  UU <- seq(1e-05, UUmax, len=201)
  K <- function(z, u) dnorm(z, 0, u)
  if(missing(perm)) {
    n <- length(Y)
    nperm <- 25
    perm <- matrix(0, ncol=nperm, nrow=n)
    perm[,1] <- 1:n
    for(j in 2:nperm) perm[,j] <- sample(n)
  }   
  else nperm <- ncol(perm)
  b.old <- as.numeric(ols$coefficients)
  f <- 1 + 0 * UU
  f <- f / int(f, UU)
  lik <- c()
  repeat {
    Z <- as.numeric(Y - X %*% b.old)
    o <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)
    lik <- c(lik, -o$L)
    ww <- o$wt
    b <- as.numeric(lm(Y ~ X - 1, weights=ww)$coefficients)
    if(sum(abs(b - b.old)) < tol) break else b.old <- b
  }
  Z <- as.numeric(Y - X %*% b)
  o <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)
  lik <- c(lik, -o$L)
  ww <- o$wt
  if(get.ci) {
    if(!exists("hessian")) library(numDeriv)   #requires the R package "numDeriv"
    prml <- function(bb) {
      Z <- as.numeric(Y - X %*% bb)
      L <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)$L
      return(L)
    }
    J <- solve(hessian(prml, b))
    ci <- b + 1.96 * outer(sqrt(diag(J)), c(-1, 1))
  } else ci <- NULL
  return(list(b=b, Yhat=as.numeric(X %*% b), ci=ci, f=o$f, UU=UU, wt=ww, lik=lik))
}


