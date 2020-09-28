
# L = \lambda = exp(r), multiplicative rate of increase
# a = AFR
# w = ALR
# R = recruitment fraction (i.e., l(\alpha))
# b = annual fertility
# p = annual survival probability

## the function from Cole (1954)
cole <- function(r,a,w,b) exp(-r) + b*exp(-r*a) - b*exp(-r*(w+1)) - 1

uniroot(cole, lower=0.01, upper=0.1, a=20, b=0.125, w=50)

## Compare great apes
# chimps
uniroot(cole, lower=0.01, upper=0.1, a=14, b=0.1, w=40)
# gorillas
uniroot(cole, lower=0.01, upper=0.2, a=8, b=0.25, w=30)
# orangutans
uniroot(cole, lower=0.01, upper=0.1, a=14, b=0.0625, w=40)
# Hutterites
uniroot(cole, lower=0.01, upper=0.2, a=20, b=0.3, w=50)

## expanded version from Slade et al. (1998)
cole1 <- function(L,a,w,R,b,p)  p*(L^-1) + R*b*(L^-a) - R*b*(p^(w - a +1))*(L^-(w+1)) - 1

uniroot(cole1, lower=1, upper=1.2, a=20, b=0.125, w=45, R=0.8, p=0.98)

# use r instead of \lambda

cole2 <- function(r,a,w,R,b,p){
  p*exp(-r) + R*b*exp(-r*a) - R*b*(p^(w - a +1))*exp(-r*(w+1)) - 1
}

uniroot(cole2, lower=0, upper=0.05, a=20, b=0.125, w=45, R=0.8, p=0.98)

# what's a good value of R? 75% survival to age 5 and then 2% annual mortality
# thereafter is not unreasonable
(R <- 0.75*0.98^15)

# Ache-like
uniroot(cole2, lower=0, upper=0.05, a=20, b=0.2, w=45, R=0.55, p=0.98)
## surprisingly close to Hill & Hurtado's observed r for the Forest Period


aaa <- seq(15,24,by=0.5)
bbb <- seq(0.1,0.25,length=19)

## this basically doesn't work. The low fertilities create problems
rvals2 <- matrix(NA,nr=19,nc=19)
for(i in 1:length(aaa)){
  for(j in 1:length(bbb)){
    rvals2[i,j] <- uniroot(cole2, lower=0.01, upper=0.1, a=aaa[i], b=bbb[j],
                      w=45, R=0.8, p=0.98)$root
  }
}

# Contour plot
contour(x=aaa,y=bbb,z=rvals2,
       lwd=2,
       col="grey",
       xlab="AFR",
       ylab="Birth Rate")


