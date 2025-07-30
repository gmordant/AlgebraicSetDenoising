#--------------------------------------
#    functions 
#--------------------------------------
rm(list = ls()) 
multiexponents <- function(dim, deg) {
  # Generate all combinations of n non-negative integers that sum to d
  combs <- expand.grid(rep(list(0:deg), dim))
  combs <- combs[rowSums(combs) == deg, ]
  combs
}
multiexponentsL<- function(dim,deg)
{
  temp<-c(0,0)
  l=1
  while (l <=deg){
    temp<-rbind(temp, multiexponents(dim,l))
    l=l+1
  }
  return(temp)
}
Coef<-function(d,k)
{
  return (choose(d, 2*k)/(2^k)*factorial(2*k)/factorial(k))
}
# This function currently only works for diag cov Gaus noise and dimension 2
# We assume X of the form n x d 
Vandermonde<- function(X, U, sigma){
  size<- dim(U)[1]
  d<- dim(X)[2]
  Mhat<- matrix(0, nrow=size, ncol=size)
  for (k in 1:size)
    for (l in k:size)
    {
      exps<-as.integer(as.vector(U[k,]+U[l,]))
      Mhat[k,l]<- InnerFunc(X,exps[1],exps[2],sigma)
    }
  toRet<- Mhat+t(Mhat)
  toRet<- toRet- diag(diag(toRet)/2)
  return (toRet)
}
InnerFunc<-function(X, k,l, sigma){
  partSum=0
  for (kappa in 0:floor(k/2))
    for (lambda in 0:floor(l/2))
      partSum= partSum + (-1)^(kappa+lambda) *Coef(k,kappa )* Coef(l,lambda) * sigma^(2*( kappa+ lambda)) * mean ( X[,1]^(k-2*kappa) *X[,2]^(l-2*lambda))
  
  return( partSum)
}

SolveOrder2<-function(q, lB=-2, uB=2, step=0.04, shift=0)
{
  
  #------   first pass --------
  # Initialize vectors for storing results
  y1a <- numeric()
  y2a <- numeric()
  y2b <- numeric()
  
  # Loop to solve the quadratic equation
  seq_i <-  seq(lB, uB, by = step)
  for (i in seq_i) {
    aa <- q[6]
    bb <- q[5] * i + q[3]
    cc <- q[1] + q[2] * i + q[4] * i^2 + shift
    disc <- bb^2 - 4 * aa * cc
    if (disc >= 0) {
      y2a <- c(y2a, (-bb + sqrt(disc)) / (2 * aa))
      y2b <- c(y2b, (-bb - sqrt(disc)) / (2 * aa))
      y1a <- c(y1a, i)
    }
  }
  #------   second pass --------
  
  # Initialize vectors for storing results
  x1a <- numeric()
  x2a <- numeric()
  x2b <- numeric()
  
  # Loop to solve the quadratic equation
  seq_i <-  seq(lB, uB, by = step)
  for (i in seq_i) {
    aa <- q[4]
    bb <- q[5] * i + q[2]
    cc <- q[1] + q[3] * i + q[6] * i^2 + shift
    disc <- bb^2 - 4 * aa * cc
    if (disc >= 0) {
      x2a <- c(x2a, (-bb + sqrt(disc)) / (2 * aa))
      x2b <- c(x2b, (-bb - sqrt(disc)) / (2 * aa))
      x1a <- c(x1a, i)
    }
  }
  
  l1<-c(x2a,x2b)
  l2<-c(x1a,x1a)
  

  s1 <- c(y1a, y1a)
  s2 <- c(y2a, y2b)
  
  s1<- c(l1,s1)
  s2<-c(l2,s2)
  
  
  return(cbind(s1,s2))
}

SemiAlg<-function(q,lambda, lB=-2, uB=2, step=0.04)
{
  Pts<-SolveOrder2(q, lB=-2, uB=2, step=step)
  for (eta in seq(-lambda, lambda, length.out=70) )
  {
    Pts<-rbind(Pts,SolveOrder2(q, lB=-2, uB=2, step=step, shift=eta))
  }
  return(Pts)
}

#--------------------------------------
# Example : Circle 
#--------------------------------------

set.seed(123)

# Number of points
n <- 600
dim<- 2
deg <- 2

# Radius of the circle
radius <- sqrt(2)

sigma<-0.4

# Generate angles evenly spaced around the circle
angles <- seq(0, 2*pi, length.out = n)

# Compute the x and y coordinates
x1 <- radius * cos(angles) 
y1 <- radius * sin(angles) 

x2 <- x1 + sigma*rnorm(n)
y2 <- y1 + sigma*rnorm(n)


# Produce the data matrix, compute the unbiased moment matrix and extract kernel.

X<-cbind(x2, y2)

# Without correction
M<- Vandermonde(X, U=multiexponentsL(dim,deg), sigma=0)
temp<- svd(M)
q <- temp$v[,6]
q <- q / q[1]

S<-SolveOrder2(q)


# With correction
svd_result <- svd(Vandermonde(X, U=multiexponentsL(dim,deg), sigma=sigma))
q1 <- svd_result$v[,6]
q1 <- q1 / q1[1]

S1<-SolveOrder2(q1)

SA<-SemiAlg(q1,lambda=n^{-1/2}*log(n))

idx<- sample(1:nrow(SA),2000)

# Plotting the results

plot(x2, y2, col = "black", pch = 19,, cex=0.4, xlab = "", ylab = "", xlim= c(-2.5,2.5), ylim = c(-2.5,2.5))
#points(SA[idx,1], SA[idx,2], col = rgb(red = 0.5, green = 0, blue = 1, alpha = 0.15), pch = 19, cex = 0.4) 
points(x1, y1, col = "green", pch = 19, cex=0.4)
points(S[,1], S[,2], col = "red", pch = 19, cex = 0.4)  # Modify 's1', 's2' as needed
points(S1[,1], S1[,2], col = "blue", pch = 19, cex = 0.4)  # Modify 's1', 's2' as needed

