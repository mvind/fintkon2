#setwd("C:/Users/fqx738/Dropbox/FinKont2")

#source("LiveAtLectures_HestonSimulation.R")

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

Params<-paste("S0=", S0, ", sqrt(V0)=",sqrt(V0),", r=",R, ", T=", timetoexp)
ModelParams<-paste("kappa =", kappa, ", theta =", theta, ", rho =", rho, ", sigma =", sigma)


Nsim<-10^4
NstepsPerYear<-1*252

Nsteps<-round(timetoexp*NstepsPerYear)
dt<-timetoexp/Nsteps

S<-rep(S0,Nsim); V<-rep(V0,Nsim)

AbsAtZero<-FALSE
if (AbsAtZero) title<-"Heston: Euler-S + |V| @0"
MaxAtZero<-FALSE
if (MaxAtZero) title<-"Heston: Euler-S + V^+ @0"
FullTruncation<-T
LeifMethod <- F
# From page 6 in https://www.dropbox.com/s/nw7uzmf8k0imq0t/LeifHestonWP.pdf?dl=0
# "The scheme that appears to produce the smallest discretization bias
if (FullTruncation) title<-"Heston\n Euler-ln(S) + full truncation-V"

ToFile<-FALSE; GraphFile="FullTruncation.png" 
IVspace<-TRUE

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
plot(S, type="l")
plot(V)
ExperimentParams<-paste("#steps/year =",NstepsPerYear, ", #paths =", Nsim, ", runtime (seconds) = ", round(RunTime,2))

plot(S_bs, type="l")

StdErrSimCall<-TrueCall<-SimCall<-Strikes<-S0+(-10:10)*2
ConfBandIVSimCall<-IVSimCall<-IVTrueCall<-SimCall

BSCall<-CV<-CV_beta<-IV_CV<-IV_CV_beta <- TrueCall


if(ToFile) png(GraphFile,width=14,height=14,units='cm',res=300)

for (i in 1:length(SimCall)) {
	## These are the payoff we SHOULD recieve, i.e. the payoff of holding S
	SimCall[i]<-exp(-R*timetoexp)*mean(pmax(S-Strikes[i],0)) # Mean of simulated payoff for each strike for S process!!
	
	StdErrSimCall[i]<-sd(exp(-R*timetoexp)*(pmax(S-Strikes[i],0)))/sqrt(Nsim) # Std. err of simulated payoff for S process 
	
	# Implied vol where obsprice is the true price for l 96,  strike, and time to mat is 1. For regular BS call
	if (IVspace) IVSimCall[i]<-BlackScholesImpVol(SimCall[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
	
	# This is implied vol for the obsprice plus the std.err of the true payoff
	if (IVspace) ConfBandIVSimCall[i]<-BlackScholesImpVol(SimCall[i]+1.96*StdErrSimCall[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
	
	# Calculated heston call price
	TrueCall[i]<-Heston.Fourier(S0,timetoexp,Strikes[i],R,0,V0,theta,kappa,sigma,rho,greek=1)
	

	# Implied vol for calculated heston price	
	if (IVspace) IVTrueCall[i]<-BlackScholesImpVol(TrueCall[i],S0,timetoexp,Strikes[i],R, q=0, opttype=1)
	
	# Remeber this is for each strike 
	# Control variate method
	# Hat prices, i.e. discounted call prices	
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

plot(SimCall, type="ol")
lines(CV, col="blue")
lines(CV_beta, col="red")

plot(Strikes, TrueCall - CV, type="ol")
lines(Strikes, TrueCall - CV_beta, type="ol", col="red")

plot(Strikes, IV_CV - IVSimCall, type="l", ylim=c(-.003, .003))
lines(Strikes, IV_CV_beta - IVSimCall, col="red")
lines(Strikes, IVTrueCall - IVSimCall, col="blue")

if(!IVspace){
	dummy<-c(min(SimCall,TrueCall),max(SimCall,TrueCall))
	plot(Strikes,SimCall,type='b',ylim=dummy,col='blue',ylab="Call price",main=title,xlab="Strike")
	points(Strikes,TrueCall,type='l')
	points(Strikes,SimCall+1.96*StdErrSimCall,type='l',lty=2,col='blue')
	points(Strikes,SimCall-1.96*StdErrSimCall,type='l',lty=2,col='blue')
	
	text(min(Strikes),min(dummy)+0.10*(max(dummy)-min(dummy)),"Black: Closef-form Heston ala Lipton",adj=0,cex=0.5)
	text(min(Strikes),min(dummy)+0.05*(max(dummy)-min(dummy)),"Blue: Sim' est' (o) + 95% conf' band' (---)", col="blue",adj=0,cex=0.5)
	
	text(max(Strikes),max(dummy),Params,adj=1, cex=0.5)
	text(max(Strikes),0.95*max(dummy),ModelParams,adj=1, cex=0.5)
	text(max(Strikes),0.90*max(dummy),ExperimentParams,adj=1, cex=0.5)
	text(max(Strikes),0.85*max(dummy),".../Dropbox/FinKont2/LiveAtLectures_HestonSimulation.R",adj=1, cex=0.5)
}

if(IVspace){
	dummy<-c(min(2*IVSimCall-ConfBandIVSimCall,IVTrueCall),max(IVTrueCall,ConfBandIVSimCall) )
	plot(Strikes,IVSimCall,type='b',col='blue',ylim=dummy, ylab="Implied volatility",xlab="Strike")
	points(Strikes,IVTrueCall,type='l')
	points(Strikes,ConfBandIVSimCall,type='l',lty=2,col='blue')
	points(Strikes,2*IVSimCall-ConfBandIVSimCall,type='l',lty=2,col='blue')
	
	text(min(Strikes),min(dummy)+0.10*(max(dummy)-min(dummy)),"Black: Closef-form Heston ala Lipton",adj=0,cex=0.5)
	text(min(Strikes),min(dummy)+0.05*(max(dummy)-min(dummy)),"Blue: Sim' est' (o) + 95% conf' band' (---)", col="blue",adj=0,cex=0.5)
	
	text(max(Strikes),max(dummy),Params,adj=1, cex=0.5)
	text(max(Strikes),0.95*max(dummy),ModelParams,adj=1, cex=0.5)
	text(max(Strikes),0.90*max(dummy),ExperimentParams,adj=1, cex=0.5)
	text(max(Strikes),0.85*max(dummy),".../Dropbox/FinKont2/LiveAtLectures_HestonSimulation.R",adj=1, cex=0.5)
}

if (ToFile) dev.off()


