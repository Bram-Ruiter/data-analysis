import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt


N = 10
CI = 0.95 #confidence interval


def k_bounds(N,p,CI): #returns the bounds for k for given p value
    binomial = binom(N,p)
    #determine lower bound for k
    p_lower = 0
    k_min = 0
    while (p_lower + binomial.pmf(k_min) < (1-CI)/2): #lower tail of probability should always be less than (1-CF)/2
        p_lower += binomial.pmf(k_min)
        k_min += 1

    #determine upper bound for k
    p_upper = 0
    k_max = N
    while (p_upper + binomial.pmf(k_max) <(1-CI)/2): #upper tail of probability should always be less than (1-CF)/2
        p_upper += binomial.pmf(k_max)
        k_max -= 1

    return p_lower,p_upper,k_min,k_max

prob = np.linspace(0,1,1000) #probability array
k_min_array = np.array([])
k_max_array = np.array([])
CF_array = np.array([])
for p in prob:
    p_lower, p_upper, k_min, k_max = k_bounds(N,p,CI)
    CF_array = np.append(CF_array,1-p_lower-p_upper)
    k_min_array = np.append(k_min_array,k_min)
    k_max_array = np.append(k_max_array,k_max)

fig,ax = plt.subplots()
ax.set_xlabel("prob")
ax.set_ylabel("k")
ax.plot(prob, k_min_array, label="Lower Bound")
ax.plot(prob,k_max_array, label="Upper Bound")
ax.legend()

ax2 = ax.twinx()
ax2.set_xlabel("prob")
ax2.set_ylabel("Confidence Interval")
ax2.plot(prob, CF_array ,label="Confidence Interval", color ='g')
ax2.set_ylim(0.80,1)


plt.show()