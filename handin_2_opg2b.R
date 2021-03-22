source("BlackScholesFormula.R")
source("Fourier.R")

timetoexp<-1.0

S0<-100
R<-0.02
V0<-0.15^2
kappa<-2
theta<-0.2^2
sigma<-1
rho<--0.5

Nsim<-10^3

# Which delta?
BS_DELTA <- F
HESTON_TRUE <- F
RISK_MIN_DELTA <- T

# STRIKE
vec_hest_fourier <-	Vectorize(Heston.Fourier, vectorize.args = c("spot", "V"))
K = 105

S<-rep(S0,Nsim); V<-rep(V0,Nsim)

Nsteps <- 250
dt<-timetoexp/Nsteps

if(BS_DELTA) {
	# initoutlay <- BlackScholesFormula(S0, timetoexp, K, R, 0, sqrt(V0), 1,1) # Should we use BS or Heston price?
	initoutlay <- Heston.Fourier(S0,timetoexp,K,R,0,V0,theta,kappa,sigma,rho,greek=1)
	Vpf <- rep(initoutlay, Nsim)
	
	a <- BlackScholesFormula(S0, timetoexp, K, R, 0, sqrt(V0), 1,2)
	b <-  Vpf - S0*a
}

if(HESTON_TRUE) {
	initoutlay <- Heston.Fourier(S0,timetoexp,K,R,0,V0,theta,kappa,sigma,rho,greek=1)
	Vpf <- rep(initoutlay, Nsim)
	
	a <- rep(Heston.Fourier(S0,timetoexp,K,R,0,V0,theta,kappa,sigma,rho,greek=2), Nsim)
	b <- Vpf - S0 * a

}
if(RISK_MIN_DELTA) {
	initoutlay <- Heston.Fourier(S0,timetoexp,K,R,0,V0,theta,kappa,sigma,rho,greek=1)
	Vpf <- rep(initoutlay, Nsim)
	
	a <- rep(Heston.Fourier(S0,timetoexp,K,R,0,V0,theta,kappa,sigma,rho,greek=2), Nsim) + rho * sigma * Heston.Fourier(S0,timetoexp,K,R,0,V0,theta,kappa,sigma,rho,greek=4) /S0
	b <- Vpf - S0 * a
	
}



for (i in 1:Nsteps){
		#dW2_old <- dW2
	cat(i, "\n")
		dW1<-sqrt(dt)*rnorm(Nsim,0,1)
		dW2<-rho*dW1+sqrt(1-rho^2)*sqrt(dt)*rnorm(Nsim,0,1)
		
	# Full truncation scheme
			X<-log(S)
			X<-X+(R-0.5*pmax(V,0))*dt+sqrt(pmax(V,0))*dW1
			S<-exp(X)
			V<-V+kappa*(theta-pmax(V,0))*dt+sigma*sqrt(pmax(V,0))*dW2
			
#		cat(V, " \n")	
		
		#S_bs <- S_bs + R*S_bs * dt + sqrt(pmax(V,0))*S_bs * dW1
		#S_bs <- S_bs * exp((R - 0.5 *sqrt(pmax(V,0))^2)*dt + sqrt(pmax(V,0)) * dW1)
		
		Vpf <- a * S + b*exp(R*dt)
		
		if(BS_DELTA) {
			a <- BlackScholesFormula(S, timetoexp -(i - 1)*dt, K, R, 0, sqrt(pmax(V,0)), 1,2)
			b <- (Vpf - a*S)
		} else if(HESTON_TRUE) {
			# Heston fourier is not vectorized
			for(j in 1:Nsim) {
				a[j] <- Heston.Fourier(S[j],timetoexp - (i - 1)*dt,K,R,0,V[j],theta,kappa,sigma,rho,greek=2)

			}
			#a <- vec_hest_fourier(S, timetoexp -(i-1)*dt, K, R, 0, V, theta, kappa, sigma, rho, greek=2)
				
			b <- Vpf - a*S
		} else if(RISK_MIN_DELTA) {
			for(j in 1:Nsim) {
				a[j] <-Heston.Fourier(S[j],timetoexp -(i-1)*dt,K,R,0,V[j],theta,kappa,sigma,rho,greek=2) + 
					rho * sigma * Heston.Fourier(S[j],timetoexp - (i-1)*dt,K,R,0,V[j],theta,kappa,sigma,rho,greek=4) /S[j]
				
			}
			
			#a <- vec_hest_fourier(S,timetoexp,K,R,0,V,theta,kappa,sigma,rho,greek=2) + rho * sigma * vec_hest_fourier(S,timetoexp,K,R,0,V,theta,kappa,sigma,rho,greek=4) /S
			b <- Vpf - a*S
			}
		
		
		
}  
dW1<-sqrt(dt)*rnorm(Nsim,0,1)
dW2<-rho*dW1+sqrt(1-rho^2)*sqrt(dt)*rnorm(Nsim,0,1)

# Full truncation scheme
X<-log(S)
X<-X+(R-0.5*pmax(V,0))*dt+sqrt(pmax(V,0))*dW1
S<-exp(X)
V<-V+kappa*(theta-pmax(V,0))*dt+sigma*sqrt(pmax(V,0))*dW2


Vpf <- a * S + b*exp(R*dt)

plot(S, Vpf)
points(0:200, pmax(0:200 - K,0), type="l", lwd=3, col="red")

j1000 <- list(S, Vpf)
j250 <- list(S, Vpf)
j10 <- list(S, Vpf)
library(scales)

j1223 <- list(S, Vpf)

par(mfrow=c(3,1))
plot(j1000[[1]], j1000[[2]], xlab="S", ylab="Payoff", main="Risk minimizing Delta, N=1000. Nsims = 1000, Nsteps = 1000")
points(0:200, pmax(0:200 - K,0), type="l", lwd=3, col="grey")

plot(j250[[1]], j250[[2]], xlab="S", ylab="Payoff", main="Risk minimizing Delta, N=1000. Nsims = 1000, Nsteps = 250")
points(0:200, pmax(0:200 - K,0), type="l", lwd=3, col="grey")

plot(j10[[1]], j10[[2]], xlab="S", ylab="Payoff",main="Risk minimizing Delta, N=1000. Nsims = 1000, Nsteps = 10")
points(0:200, pmax(0:200 - K,0), type="l", lwd=3, col="grey")

get_hedge_error <- function(j) {
	j[[2]] - pmax(j[[1]] - K,0)	
}

mean(get_hedge_error(j1000))
mean(get_hedge_error(j250))
mean(get_hedge_error(j10))

mean(get_hedge_error(j1223))



# Risk minimzing hegde error 
# 1000 = -0.04830919
# 250 = -0.05594793
# 10 = -0.8496946

# True delta
# 1000 = 0.06879613
# 250 = 0.1934852 #-0.05020561
# 10 = -0.4849947

# BS Delta
# 1000 = -0.1038097
# 250 = 0.07402711
# 10 = -0.7622119