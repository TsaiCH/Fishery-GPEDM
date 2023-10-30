######################################################################################################
## Single-species Pella-Tomlinson production model
## parms from Thorson et al.2012 CJFAS
######################################################################################################
n <- 1.478 # (among-stock sd: 0.844, stocks' min:0.3 max:2.9)
K <- 1
m <- 0.404  # (among-stock sd: 0.135, stocks' min: 0.15 max: 0.75)
f <- seq(0.3,1.6,length.out=100)
#f <- c(seq(0.3,1.6,length.out=50),seq(1.6,0.3,length.out=50))
B.all <- C.all <- c()
for(j in 1:150){
  B <- K
  s <- b <- catch <- c()
  for(i in 1:length(f)){
    P <- (n^(n/(n-1))/(n-1))*m*(B[i]/K-(B[i]/K)^n)
    B[i+1] <- (B[i]+P-f[i]*B[i]) * exp(rnorm(1,0,0.1))
    s[i] <- P/m
    b[i] <- B[i]/K
    catch[i] <- f[i]*B[i]
  }
  B.all <- cbind(B.all,B)
  C.all <- cbind(C.all,catch)
}
B.all <- tail(B.all,100)
C.all <- tail(C.all,100)
############################################################################################################################################
save(B.all,C.all,f,file="PellaTomlinson_scenario_increasing_example.RData")

############################################################################################################################################
## Ricker-type prey-predator dynamics with prey harvested
## Theoretical MSY : s=[0.3,1.6]; exp(s)*r2*(r1-d-s)/(d^2+r1*r2) ~ 1.362
###################################################################################################################################################
f <- seq(0.3,1.6,length.out=100)
#f <- c(seq(0.3,1.6,length.out=50),seq(1.6,0.3,length.out=50))
c <- 0.1
rx <- 2.2
ry <- 1.8
nsim <- 150
B1.all <- B2.all <- C.all <- c()
for(j in 1:nsim){
  x <- y <- catch <- c()
  x[1] <- y[1] <- 0.5
  for(i in 1:length(f)){
    x[i+1] <- x[i]*exp(rx*(1-x[i])-c*y[i]-f[i])*exp(rnorm(1,0,0.1))
    y[i+1] <- y[i]*exp(ry*(1-y[i])+c*x[i])*exp(rnorm(1,0,0.1))
    catch[i] <- x[i]*exp(f[i])
  }
  B1.all <- cbind(B1.all,x)
  B2.all <- cbind(B2.all,y)
  C.all <- cbind(C.all,catch)
}
prey.all <- tail(B1.all,100)
predator.all <- tail(B2.all,100)
C.all <- tail(C.all,100)
########################################################################################################################################
save(prey.all,predator.all,C.all,f,file="PreyPredator_scenario_increasing_example.RData")

