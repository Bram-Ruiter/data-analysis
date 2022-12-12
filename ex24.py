import numpy as np
from scipy.stats import poisson
import  matplotlib.pyplot as plt

b = 3
x=4 #take x=4 as an example
CL = 0.90

def R(x,s): #returns the feldman cousin ratio
    poisson1 = poisson.pmf(x,s+b)
    #consider the physical limitations for s:
    if x<b:
        poisson2 = poisson.pmf(x,b)
    else:
        poisson2 = poisson.pmf(x,x)
    return poisson1/poisson2

def n_bounds(x,s,CL): #returns the bounds for #measurements n for a given s value
    #determine lower bound for s
    n_max = 50 #asume that no more than 50 events can be measured with any reasonable probaility
    ratios = np.array([])
    for x in range(n_max):
        ratios = np.append(ratios, R(x,s))
    ratios_sorted_index = np.argsort(ratios)[::-1 ] #get the indices of the array elements in ratios in descending order
    n_array = np.array([])
    poisson_cum = 0
    i=0
    #sum over all values n in descending order based on the R ratio.
    while (poisson_cum<CL):
        n_array = np.append(n_array, ratios_sorted_index[i])
        poisson_cum += poisson.pmf(ratios_sorted_index[i],s+b)
        i+=1
    x_min = np.amin(n_array)
    x_max = np.amax(n_array)
    return poisson_cum,x_min,x_max

s_array = np.linspace(0,10,1000)
x_min_array = np.array([])
x_max_array = np.array([])
CL_array = np.array([])
for s in s_array:
    CF, x_min, x_max = n_bounds(x,s,CL)
    CL_array = np.append(CL_array,CF)
    x_min_array = np.append(x_min_array,x_min)
    x_max_array = np.append(x_max_array,x_max)

fig,ax = plt.subplots()
ax.set_xlabel("n")
ax.set_ylabel("s")
ax.plot(x_min_array, s_array, label="Lower Bound")
ax.plot(x_max_array, s_array, label="Upper Bound")
ax.legend()
ax.fill_betweenx(s_array, x_min_array,x_max_array, alpha=0.2)
ax.set_xlim(0,12)
ax.set_ylim(0,10)


#legacy: confidence level plot for when plotting with s on the x-axis
#ax2 = ax.twinx()
#ax2.set_xlabel("prob")
#ax2.set_ylabel("Confidence Level")
#ax2.plot(s_array, CL_array ,label="Confidence Level", color = 'g')
#ax2.set_ylim(0.80,1)

plt.show()


def n_bounds_old(x,s,CL): #returns the bounds for #measurements n for a given s value
    #determine lower bound for s
    p_lower = 0
    x_min = 0
    while (p_lower + R(x_min,s) < (1-CL)/2): #lower tail of probability should always be less than (1-CF)/2
        p_lower +=  R(x_min,s)
        x_min += 1

    #determine upper bound for s
    p_upper = 0
    x_max = 50  #take 50 since the chance to measure more than that is just negligible
    while (p_upper +  R(x_max,s) <(1-CL)/2): #upper tail of probability should always be less than (1-CF)/2
        p_upper +=  R(x_max,s)
        x_max -= 1

    return p_lower,p_upper,x_min,x_max