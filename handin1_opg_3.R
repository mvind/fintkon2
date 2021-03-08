


## Lets first estimate r* 

K = 4.5


# Vasicek parameters
# HANDIN NOTATION
#      kappa   theta^Q       sigma
# dr = 0.1 * (0.05 - r) dt + 0.015 r dW
# r(0) = r0

# Convert to book notation
a <- 0.1
b <- 0.1 * 0.05
sig <- 0.015
r0 <- 0.02

# prop 24.3

v_p <- function(t, capT, r, a, b, sigma) {
	# B function
	B <- (1/a) * (1 - exp(-a*(capT - t)))
	
	# A function
	A_1 <- ((B - capT + t) * (a*b - 0.5 * sigma^2))	/ a^2
	A_2 <- ((sigma^2) * (B^2)) / (4*a)
	A <- A_1 - A_2
	
	exp(A - B*r)
}


# Build optim function
op <- function(r_star) {
	sum <- v_p(1, 2, r = r_star, a = a, b = b, sigma = 0.015) +
		v_p(1, 3, r = r_star, a = a, b = b, sigma = 0.015) +
		v_p(1, 4, r = r_star, a = a, b = b, sigma = 0.015) +
		v_p(1, 5, r = r_star, a = a, b = b, sigma = 0.015) + 
		v_p(1, 6, r = r_star, a = a, b = b, sigma = 0.015)
	
	K - sum
}

# v_p(0, 6, 0.06, a = a, b = b, sigma = 0.015)
sol <- uniroot(op, c(0, 1))

r_star <- sol$root

# We need to simulate r(t) 
## simulate short rate paths
n <- 1   # MC simulation trials
T <- 6    # total time
m <- 6   # subintervals
dt <- T/m  # difference in time each subinterval

r <- matrix(0,m+1,n)  # matrix to hold short rate paths
r[1,] <- r0

k = 0.1 
theta = 0.05

for(j in 1:n){
	for(i in 2:(m+1)){
		dr <- k*(theta-r[i-1,j])*dt + sig*sqrt(dt)*rnorm(1,0,1)
		r[i,j] <- r[i-1,j] + dr
	}
} 

plot(r, type="l")
## NO WE DONT WE HAVE r0 !!!

# Calculate K_i 
# Use v_p with r = r_star 

# Calculate bond options prop 24.9 (note there is refered to the same formula in chap 26)
# USING 26.13
# T_1 = 1 = capT 
# T_2 = 2


call_price <- function(T1, T2) {
	sig_2  <- ((sig^2) / (2*a^3)) * (1 - exp(-2*a*T1)) * 
		(1 - exp(-a*(T2 - T1)))^2
	
	p_0_2 <- v_p(0, T2, r = r0, a = a, b = b, sigma = sig)
	p_0_1 <- v_p(0, T1, r = r0, a = a, b = b, sigma = sig)
	
	K_2 <- v_p(1, T2, r = r_star, a = a, b = b, sigma = sig)
	d2 <-( log( p_0_2 / (K_2 * p_0_1) ) - 0.5 * sig_2) / sqrt(sig_2)
	d1 <- d2 + sqrt(sig_2)
	
	# time 0 price 
	p_0_2 * pnorm(d1) - K_2 * p_0_1 * pnorm(d2)
	
	
}

sum(call_price(T1 = 1, T2 = 2:6))
