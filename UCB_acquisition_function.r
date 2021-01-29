#For computing the UCB
#From Brochu 2010, "A tutorial on Bayesian Optimization of Expensive Cost Functions, with Application to Active  User Modelling and Hierarchical Reinforcement Learning"
#and from Srinivas "Gaussian Process Optimization in Bandit Setting: No Regret and experimental Design"

tau_t_calc = function(t,d,delta){
	tau_t = 2*log(
		t^(d/2+2)*pi^2/(3*delta)
	)
	return(tau_t)
}

kappa_calc = function(v=1,t,d,delta){ # can change v, the original refernece set it to 1
	tau_t = tau_t_calc(t,d,delta)
	kappa = sqrt(v*tau_t)
	return(kappa)}
 
#kappa is the multiplier for the variance