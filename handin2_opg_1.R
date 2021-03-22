

Bach_price <- function(spot, strike, timetomat, sigma) {
	d <- (spot - strike)/(sigma*sqrt(timetomat))
	(spot - strike)*pnorm(d) + sigma*sqrt(timetomat)*dnorm(d)
}

impl <- function(obs_price, spot, strike, timetomat) {
	diff <- function(sigB, obs_price, spot, strike, timetomat) {
		v <- abs(Bach_price(spot, strike, timetomat, sigB) - obs_price)
	}
	
	uniroot(diff, c(10^(-6),1000),obs_price=obs_price,spot=spot,timetomat=timetomat,strike=strike)$root
}

impl_vec <- Vectorize(impl)

S0 <- 100
sigma<-15
Strikes<-S0+(-10:10)*2

obs_price <- Bach_price(S0, 95, 1, sigma)

fff <- function(x) {
	Bach_price(S0, 100, 1, 10) - Bach_price(S0, 95, 1, x)
}
plot(fff, xlim = c(0, 50))
points(0:50, rep(0, 51), type="l")
impl(obs_price, 100, 95, 1)



impl_vec(obs_price, 100, 95, 1)

-Bach_price(S0, 94, 1, 15) + obs_price

out <- impl_vec(obs_price, S0, Strikes, 1)


plot(95:115, sqrt(out), type="l", ylab="Implied Volatility (%)", xlab="Strike")


plot(75:115, sqrt(sig_lv), type="l")


lines(75:115, rep(15/100, 41), col="red")
source("BlackScholesFormula.R")
BachelierFormula  <- function (spot,timetomat,strike, sigma){
	d<-(spot-strike)/(sigma*sqrt(timetomat))
	BachelierFormula<-(spot-strike)*pnorm(d)+sigma*sqrt(timetomat)*dnorm(d)
}



S0<-100
sigmaBach<-15
timetomat<-1
Strikes<-S0+(-20:20)*2
impvol<-1:length(Strikes)

for(i  in 1:length(impvol)){
	Bachprice<-BachelierFormula(S0,timetomat,Strikes[i],sigmaBach)
	impvol[i]<-BlackScholesImpVol(Bachprice,S0,timetomat,Strikes[i],r=0, q=0, opttype=1)
}

plot(Strikes,impvol,xlab="Strikes",ylab="Implied Vol",main="Implied volatility in the Bachelier model", type="l")

library(numDeriv)

# Partial T
BP_T <- function(capT, spot = 100, strike = 95, sigma = 15) {
	Bach_price(spot, strike, capT, sigma)
}

grad(BP_T, capT=10, strike = 75)


# Partial K 

BP_K <- function(strike, capT = 1, spot = 100, sigma = 15) {
	Bach_price(spot, strike, capT, sigma)
}

grad(BP_K, 95)

BP_dK2 <- function(strike) {
	grad(BP_K, strike)
}


N <- length(Strikes)
sig_lv <- rep(NA, N)

for(j in 1:N) {
	i <- j + 74
	sig_lv[j] <- ( grad(BP_T, 1,strike = i)) / ( 0.5 * ((i)^2) * grad(BP_dK2, i))
}

plot(Strikes, sqrt(sig_lv) , type="l", lty=1, ylab="sigma", main="Local volatility in the Bachelier model")
lines(Strikes, impvol, col="red", type="l", lty=1)
abline(h=0.15, col="black", lty=5)
legend("topright", c("Local volatility per. Dupires", "Implied volatility", "Theoretical local volatility"), col=c("black", "red","black"), lty=c(1, 1,5))


impvol2<-1:length(Strikes)

for(i  in 1:length(impvol)){
	Bachprice<-BachelierFormula(S0,timetomat,Strikes[i],sig_lv[i]*S0)
	impvol2[i]<-BlackScholesImpVol(Bachprice,S0,timetomat,Strikes[i],r=0, q=0, opttype=1)
}

