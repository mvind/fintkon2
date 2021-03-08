BlackScholesFormula  <- function (spot,timetomat,strike,r, q=0, sigma, opttype=1, greektype=1)
{ 
	
	d1<-(log(spot/strike)+ ((r-q)+0.5*sigma^2)*timetomat)/(sigma*sqrt(timetomat))
	d2<-d1-sigma*sqrt(timetomat)

	if (opttype==5 && greektype==1) result <- strike * pnorm(-d2) - spot * exp(r * timetomat) * pnorm(-d1)
	
	if (opttype==5 && greektype==2) result <- pnorm(d1) - 1
	# if (opttype==1 && greektype==1) result<-spot*exp(-q*timetomat)*pnorm(d1)-strike*exp(-r*timetomat)*pnorm(d2)
	# 
	# if (opttype==2 && greektype==1) result<- -spot*exp(-r*timetomat)*pnorm(-d1)+strike*exp(-q*timetomat)*pnorm(-d2)
	# 
	# if (opttype==1 && greektype==2) result<-exp(-q*timetomat)*pnorm(d1)
	# 
	# if (opttype==2 && greektype==2) result<-exp(-q*timetomat)*(pnorm(d1)-1)
	# 
	# if (greektype==3) result<-exp(-q*timetomat)*dnorm(d1)/(spot*sigma*sqrt(timetomat))
	# 
	# if (greektype==4) result<-exp(-q*timetomat)*spot*dnorm(d1)*sqrt(timetomat)
	# 
	BlackScholesFormula<-result
	
}

Y0 <- 1/100
X0 <- 1/110
S0 <- 50 
	
# r_US <- 0.05
# r_J <- 0.03

r_US <- 0.05
r_J <- 0.03

sig_X <- c(0.1, 0.02)
sig_J <- c(0, 0.25)

capT<-1
strike<-45

Nhedge<-12
Nrep<-1000
dt<-capT/Nhedge

q_drift<- r_J - c(matrix(sig_X, nrow = 1) %*% matrix(sig_J))
q_diff <- sqrt(sum(sig_J^2))

St <- rep(S0, Nrep)
Xt <- rep(X0, Nrep)

initialoutlay<- Y0 * exp(-r_US*capT) * BlackScholesFormula(S0,capT,strike, r = q_drift,q=0,sigma=q_diff,5,1)

Vpf<-rep(initialoutlay,length=Nrep)

# a is in dollars 
a<-Y0 * exp((q_drift - r_US) * capT) * BlackScholesFormula(S0,capT,strike, r = q_drift,q=0,sigma=q_diff,5,2)

c <- -a * S0

# So S0 * X0 is in dollars aswell
b<-Vpf-(a*S0)
b<-Vpf-(a*S0) - c

for(i in 2:Nhedge){
	sig_f_noise <- rep(NA, Nrep)
	# simulate Xt
	sig_x_noise <- rep(NA, Nrep) 
	for(j in 1:Nrep) {
		w <- matrix(c(rnorm(1), rnorm(1)), ncol = 1, nrow = 2)
		sig_f <- matrix(sig_J, nrow = 1, ncol = 2)
		
		sig_f_noise[j] <-  sig_f %*% w
		
		sig_x <- matrix(sig_X, nrow = 1, ncol = 2)
		
		sig_x_noise[j] <- sig_x %*% w
	}
	W <- rnorm(Nrep)
	
	St<-St*exp((q_drift-0.5*q_diff^2)*dt +sqrt(dt)*sig_f_noise)	
	Xt<-Xt*exp((r_US - r_J -0.5*(0.1)^2)*dt +sqrt(dt)*sig_x_noise)	
	
	Vpf<- a * St + b * exp(dt*r_US) 
	Vpf <- Vpf + c * exp(dt*r_J)
	a<- Y0 * exp((q_drift - r_US) * capT) * BlackScholesFormula(St,(capT-(i-1)*dt),strike, r = q_drift,q=0,q_diff,5,2)
	# b<- Vpf - a*St
	c <- -a * St
	b<- Vpf - a*St - c
}

St<-St*exp((q_drift-0.5*q_diff^2)*dt +q_diff*sqrt(dt)*W)	
Xt<-Xt*exp((r_US - r_J -0.5*(0.1)^2)*dt +(0.1)*sqrt(dt)*W)	

Vpf<-a*St+b*exp(dt*r_US) 
Vpf <- Vpf + c * exp(dt*r_J)

plot(St, Vpf)
points(0:100,Y0*pmax(strike - 0:100,0),type='l',lwd=2, col="red")


c12 <- c(St, Vpf)
c100 <- c(St, Vpf)
c1000 <- c(St, Vpf)

hedgeerror<-(Vpf-Y0 * pmax(strike - St,0))
mean(hedgeerror)

plot(c12[1:1000], c12[1001:2000], xlab = "S_j(T)", ylab="Value of hedge portfolio", 
		 main="Discrete hedge experiment of quanto-put \n (N=1000, H = #Equidistant hedge points)")
points(c100[1:1000], c100[1001:2000], col="blue")
points(c1000[1:1000], c1000[1001:2000], col="green")
points(0:100,Y0*pmax(strike - 0:100,0),type='l',lwd=2, col="red")
legend(84, .195, c("H = 12", "H = 100", "H = 1000", "Theoretical Payoff"), col = c("black", "blue", "green", "red"), lwd=1, lty=1)


plot(c12[1001:2000], c12[1001:2000] - Y0 * pmax(strike - c12[1:1000]), type="p")
