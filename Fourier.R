BS.Fourier <- function(spot,timetoexp,strike,r,divyield,sigma,greek=1){
X<-log(spot/strike)+(r-divyield)*timetoexp

integrand<-function(k){
	integrand<-Re(exp( (-1i*k+0.5)*X- 0.5*(k^2+0.25)*sigma^2*timetoexp)/(k^2+0.25))
} 
dummy<-integrate(integrand,lower=-Inf,upper=Inf)$value

BS.Fourier<-exp(-divyield*timetoexp)*spot-strike*exp(-r*timetoexp)*dummy/(2*pi)
}


Heston.Fourier <- function(spot,timetoexp,strike,r,divyield,V,theta,kappa,epsilon,rho,greek=1){
	
	X<-log(spot/strike)+(r-divyield)*timetoexp
	kappahat<-kappa-0.5*rho*epsilon
	xiDummy<-kappahat^2+0.25*epsilon^2
	
	integrand<-function(k){
		xi<-sqrt(k^2*epsilon^2*(1-rho^2)+2i*k*epsilon*rho*kappahat+xiDummy)
		Psi.P<--(1i*k*rho*epsilon+kappahat)+xi
		Psi.M<-(1i*k*rho*epsilon+kappahat)+xi
		alpha<--kappa*theta*(Psi.P*timetoexp + 2*log((Psi.M + Psi.P*exp(-xi*timetoexp))/(2*xi)))/epsilon^2
		beta<--(1-exp(-xi*timetoexp))/(Psi.M + Psi.P*exp(-xi*timetoexp))
		numerator<-exp((-1i*k+0.5)*X+alpha+(k^2+0.25)*beta*V)
		
		if(greek==1) dummy<-Re(numerator/(k^2+0.25))
		if(greek==2) dummy<-Re((0.5-1i*k)*numerator/(spot*(k^2+0.25)))
		if(greek==3) dummy<--Re(numerator/spot^2)
		if(greek==4) dummy<-Re(numerator*beta)
		
		integrand<-dummy
	} 
	
	
	dummy<-integrate(integrand,lower=-100,upper=100,stop.on.error = FALSE)$value
	
	if (greek==1) dummy<-exp(-divyield*timetoexp)*spot-strike*exp(-r*timetoexp)*dummy/(2*pi)
	
	if(greek==2) dummy<-exp(-divyield*timetoexp)-strike*exp(-r*timetoexp)*dummy/(2*pi)
	
	if(greek==3) dummy<--strike*exp(-r*timetoexp)*dummy/(2*pi)
	
	if(greek==4) dummy<--strike*exp(-r*timetoexp)*dummy/(2*pi)
	
	Heston.Fourier<-dummy
}

Andreasen.Fourier <- function(spot,timetoexp,strike,Z,lambda,beta,epsilon){
	
	X<-log(spot/strike)
	
	integrand<-function(k){
		neweps<-lambda*epsilon
		xi<-sqrt(k^2*neweps^2+beta^2+0.25*neweps^2)
		Psi.P<--beta+xi
		Psi.M<-beta+xi
		A<--beta*(Psi.P*timetoexp + 2*log((Psi.M + Psi.P*exp(-xi*timetoexp))/(2*xi)))/(epsilon^2)
		B<-(1-exp(-xi*timetoexp))/(Psi.M + Psi.P*exp(-xi*timetoexp))
		integrand<-Re(exp( (-1i*k+0.5)*X+A-(k^2+0.25)*B*lambda^2*Z)/(k^2+0.25))
	} 
	
	dummy<-integrate(integrand,lower=-Inf,upper=Inf)$value
	
	Andreasen.Fourier<-spot-strike*dummy/(2*pi)
}