## ----setup, include=FALSE---------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(tinytex.verbose = TRUE)


## ---- cache=TRUE------------------------------------------------------
library(demogR)
tw <- read.table("data/taiwan1906.txt", header=TRUE)
A <- leslie.matrix(tw$lx, tw$mx, L=FALSE, peryear=2.5, infant.class=FALSE)
## calculate a bunch of quantities
ea <- eigen.analysis(A)
## in case you wondered if Taiwan (1906) was unusual
## we can calculate the annual growth rate
log(ea$lam)/2.5
eee <- ea$elas
plot(eee,ylab="Elasticity", peryear=2.5)
legend("topright", c("survival","fertility"),lwd=2, col=c("black","grey"))


## ---- cache=TRUE------------------------------------------------------
## extract the subdiagonal
elas <- eee[row(eee)==col(eee)+1]
age <- cumsum(rep(2.5,20))
plot(age,1/elas, type="l", lwd=2, axes=FALSE, frame=TRUE,
     xlab="Age", ylab="Mortality Rate")
axis(1)


## ---- cache=TRUE------------------------------------------------------
# first, calculate eigenvalues/vectors of A
ev <- eigen(A)
# find the dominant eigenvalue and extract
# make sure only to take real part (since dominant eigenvalue is real)
lmax <- which(Re(ev$values)==max(Re(ev$values)))
lambda <- Re(ev$values[lmax])
# matrix of eigenvectors
U <- ev$vectors
# dominant eigenvector
u <- abs(Re(U[,lmax]))
# matrix of left eigenvectors (complex conjugate of the inverse of U)
V <- Conj(solve(U))
# dominant left eigenvector
v <- abs(Re(V[lmax,]))
# outer product of v and u
s <- v%o%u
# non-existent transitions -> 0 
s[A == 0] <- 0
(ssurv <- s[row(s)==col(s)+1])
plot(age,ssurv, type="l", lwd=2, col="orange", 
     xlab="Age", ylab="Fitness Sensitivity")
# weird



## ---- cache=TRUE------------------------------------------------------
e <- s*A/lambda
(esurv <- e[row(e)==col(e)+1])
## a neat trick
sum(e)
plot(age,esurv, type="l", lwd=2, col="cornflowerblue", 
     xlab="Age", ylab="Fitness Elasticity")


## Plot tangent to curve
f <- function(x) {
  1 - exp(-1*(x-1))
}
# derivative of the constraint function
fprime <- function(x) {
  exp(-x)*exp(1)
}

x <- seq(1.5,3.25,length=100)
plot(2.5,0.6, type="n", axes=FALSE, 
     xaxs="i", yaxs="i",
     xlab=expression(a[ij]), 
     ylab="", 
     xlim=c(1,4), ylim=c(0,1))
mtext(side=2, text=expression(lambda), las=1, line=2, cex=1.5)
curve(1-exp(-1*(x-1)), from=1.5, to=3.5, lwd=2, add=TRUE)
box()
lines(x, f(2.5)+fprime(2.5)*(x-2.5), lwd=1, col="red")
segments(2.5,0,2.5,f(2.5),lty=3, col="black")
segments(0,f(2.5),2.5,f(2.5), lty=3, col="black")


