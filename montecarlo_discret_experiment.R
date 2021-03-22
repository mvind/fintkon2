source("BlackScholesFormula.R")
source("Fourier.R")

timetoexp<-1.0

S0<-100
R<-0.02
V0<-0.15^2
kappa<-2
theta<-0.2^2
sigma<-1.0
rho<--0.5

AbsAtZero<-F
MaxAtZero<-F
FullTruncation<-F
LeifMethod <- T

sim_exp <- function() {
	S_bs <- rep(S0, Nsim)
	S<-rep(S0,Nsim); V<-rep(V0,Nsim)
	RunTime<-system.time(
		for (i in 1:Nsteps){
			#dW2_old <- dW2
			dW1<-sqrt(dt)*rnorm(Nsim,0,1)
			dW2<-rho*dW1+sqrt(1-rho^2)*sqrt(dt)*rnorm(Nsim,0,1)
			
			S_bs <- S_bs + R*S_bs * dt + sqrt(V0)*S_bs * dW1
			
			if (AbsAtZero){ 
				S<-S+R*S*dt+sqrt(V)*S*dW1
				V<-abs(V + kappa*(theta-V)*dt + sigma*sqrt(V)*dW2)
			}
			if (MaxAtZero){ 
				S<-S+R*S*dt+sqrt(V)*S*dW1
				V<-pmax(V + kappa*(theta-V)*dt + sigma*sqrt(V)*dW2,0)
			}
			
			if (FullTruncation){
				X<-log(S)
				X<-X+(R-0.5*pmax(V,0))*dt+sqrt(pmax(V,0))*dW1
				S<-exp(X)
				V<-V+kappa*(theta-pmax(V,0))*dt+sigma*sqrt(pmax(V,0))*dW2
			}
			
			if(LeifMethod) {
				S<-S+R*S*dt+sqrt(V)*S*dW1
				#X<-log(S)
				#X<-X+(R-0.5*V)*dt+sqrt(V)*dW1
				#S<-exp(X)
				
				numerator <- 0.5*sigma^2 * ((kappa)^(-1)) * V * (1 - exp(-2*kappa*dt))
				denumerator <- (exp(-kappa * dt) * V + (1 - exp(-kappa*dt))*theta)^2
				gamma2 <- ((dt)^(-1)) * log(1 + numerator/denumerator)
				
				V <- (exp(-kappa * dt) * V + (1 - exp(-kappa*dt))*theta) * exp(-0.5 * gamma2 * dt + sqrt(gamma2) * (dW2))
			}
		}  
	)[3]
	
	
	StdErrSimCall<-TrueCall<-SimCall<-Strikes<-S0+(-10:10)*2
	ConfBandIVSimCall<-IVSimCall<-IVTrueCall<-SimCall
	
	BSCall<-CV<-CV_beta<-IV_CV<-IV_CV_beta <- TrueCall
	
	for (i in 1:length(SimCall)) {
		## These are the payoff we SHOULD recieve, i.e. the payoff of holding S
		SimCall[i]<-exp(-R*timetoexp)*mean(pmax(S-Strikes[i],0)) # Mean of simulated payoff for each strike for S process!!
		
		StdErrSimCall[i]<-sd(exp(-R*timetoexp)*(pmax(S-Strikes[i],0)))/sqrt(Nsim) # Std. err of simulated payoff for S process 
		
		# Implied vol where obsprice is the true price for l 96,  strike, and time to mat is 1. For regular BS call
		IVSimCall[i]<-BlackScholesImpVol(SimCall[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
		
		# This is implied vol for the obsprice plus the std.err of the true payoff
	 	ConfBandIVSimCall[i]<-BlackScholesImpVol(SimCall[i]+1.96*StdErrSimCall[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
		
		# Calculated heston call price
		TrueCall[i]<-Heston.Fourier(S0,timetoexp,Strikes[i],R,0,V0,theta,kappa,sigma,rho,greek=1)
		
		
		# Implied vol for calculated heston price	
		IVTrueCall[i]<-BlackScholesImpVol(TrueCall[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
		
		BSCall[i] <- BlackScholesFormula(S0, timetoexp, Strikes[i], R, 0, sqrt(V0), 1, 1)
		
		# hat g
		P_bs <- exp(-R*timetoexp)*mean(pmax(S_bs - Strikes[i], 0))
		P_he <- exp(-R*timetoexp)*mean(pmax(S-Strikes[i],0))
		
		CV[i] <- P_he + BSCall[i] - P_bs
		
		
		# Lets calculate beta
		beta <- cov(exp(-R*timetoexp)*pmax(S-Strikes[i],0), exp(-R*timetoexp)*pmax(S_bs-Strikes[i],0)) / var(exp(-R*timetoexp)*pmax(S_bs-Strikes[i],0))
		
		CV_beta[i] <-  P_he + beta * (BSCall[i] - P_bs)
		
		IV_CV[i]<-BlackScholesImpVol(CV[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
		IV_CV_beta[i]<-BlackScholesImpVol(CV_beta[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
		
	}
	IVSimCall - IV_CV_beta
}


out <- list()
for(k in 1:100) {
	out[[k]] <- sim_exp()
}

ab <- matrix(unlist(out), ncol = 21, byrow = T)


plot(Strikes, colSums(ab)/100, type="ol")

leif_montecarlo <- colSums(ab)/100
full_trun_montecarlo <- colSums(ab)/100
abs_at_zero_montecarlo <- colSums(ab)/100 
max_at_zero_montecarlo <- colSums(ab)/100


ymin <- min(leif_montecarlo, full_trun_montecarlo, abs_at_zero_montecarlo, max_at_zero_montecarlo)
ymax <- max(leif_montecarlo, full_trun_montecarlo, abs_at_zero_montecarlo, max_at_zero_montecarlo)

plot(leif_montecarlo, type="ol",ylim = c(ymin, ymax))
lines(full_trun_montecarlo, type="ol", col="red")
lines(abs_at_zero_montecarlo, type="ol", col="green")
lines(max_at_zero_montecarlo, type="ol", col="blue")

write(max_at_zero_montecarlo, file="max_at_zero_montecarlo")

read("abs_at_zero_montecarlo")

write(abs)

leif_mp_CV_IV <- colSums(ab)/100
leif_mp_CV_IV_beta <- colSums(ab)/100


write(leif_mp_CV_IV_beta, file="leif_mp_CV_IV_beta")

plot(Strikes, abs_at_zero_mp, type="l", ylim = c(-0.04, 0.05))
plot(Strikes, abs_at_zero_mp_CV_IV, col="red", type="l")
lines(Strikes, abs_at_zero_mp_CV_IV_beta, col="red", type="l")