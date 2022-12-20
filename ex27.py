import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson as poisson
from scipy.optimize import minimize_scalar as minimize_scalar
import scipy.stats

N_intervals = np.array([1042,860,307,78,15,3,0,0,0,1])
N_interactions = np.arange(0,10,1)
N=np.sum(N_intervals)
dof = 11-1-1 #n=11 measurements, 1 parameter(time interval)


def mean(N_interactions, N_intervals):
    tot_int = 0
    sum = 0
    for i in N_interactions:
        sum += N_interactions[i]*N_intervals[i]
        tot_int += N_intervals[i]
    return sum/tot_int

mu = mean(N_interactions,N_intervals)
print(mu)

def chi2(mean,N_interactions,N_intervals):
    chi2 = 0
    for i in N_interactions:
        N_exp = poisson.pmf(N_interactions[i],mean)*N
        chi2 += (N_intervals[i]-N_exp)**2/N_exp
    return chi2

def chi2_2(mean,N

chi = chi2(mu,N_interactions,N_intervals)
print(chi)

mu_chi=minimize_scalar(chi2,args=(N_interactions,N_intervals))
p_chi = 1 - scipy.stats.chi2.cdf(mu_chi.fun,dof,scale=N)
print(mu_chi)
print(p_chi)

plt.plot(N_interactions,poisson.pmf(N_interactions,mu_chi.x), label="chi 2 fit")
plt.plot(N_interactions,poisson.pmf(N_interactions,mu), label="estimator fit")
plt.plot(N_interactions, N_intervals/N, label="measured")
plt.legend()
plt.show()

#now for 2 poissons


##


