##import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os
from scipy.optimize import minimize as minimize
from scipy.optimize import minimize_scalar as minimize_scalar

##generate data
Afb=0.015

def pdf(x, A=Afb): #the pdf which follows from the formula for the differential cross section
    return 3/8*(1+x**2)+A*x


def pdf_background(x, A=Afb): #pdf corrected for 20% uniform background
    return pdf(x,A)*0.8+(1/2)*0.2


def generate_data(N, A=Afb): #Generates data on the domain -1,1 following the pdf
    x=np.random.uniform(0,1,N)
    y = np.array([]) #empty array to store results
    for i in range(N):
        p = [1,4*A,3,4-4*A-8*x[i]] #define the polynomial
        for root in np.roots(p): #solve the polynomial. There are only real roots in our domain, so select for this.
            if abs(root.imag)<1e-5: #if the root is not complex, store the real part of the result.
                y = np.append(y,root.real)
    return y


##see if the generated data indeed follows the pdf
#plot some generated data
N_random=10000 #number of randomly generated points
y = generate_data(N_random)
nbins=50
count, bins, ignored = plt.hist(y, bins=nbins, histtype='bar', label='Sampled distribution')
binwidth = (bins[nbins]-bins[0])/nbins

#plot the pdf
xmin=-1
xmax=1
u = np.linspace(xmin,xmax,1000)
v=pdf(u)*(N_random*binwidth) #make the plot and account for the histogram not being normalized
plt.plot(u,v, label ='pdf')
plt.legend()
plt.xlabel('cos($\\theta$)')
plt.ylabel("Count")
plt.title("Generated values")
print("\n A sample plot of generated data and the corresponding pdf. \n")
plt.show()



## estimation using method of moments
def estimator0(data): #estimate A for a given data set assuming no background
    N =np.size(data)
    M1 = np.sum(data)/N #M1 is the sum of all xi divided by N
    A = 3/2*M1
    var_M1 = 1/N*(2/5-4/9*A**2)
    var_A = (3/2)**2*var_M1
    z_score = A/(var_A**0.5)
    p_value = scipy.stats.norm.sf(abs(z_score)) #p_value from z-score
    return A, var_A, p_value

#determine the p_value from the z-score using that A is assymptotically normally distributed. We can then determine the p_value using that A is normally distributed around A=0 with a variance of var_A. The p-value can then be calculated by determining the area under the normal distribution "to the right" of the measured A value. This is a one tailed test. N>>30 so we don't need to use the t distribution.

def estimator_background(data): #estimate A for a given data set assuming 20% uniform background
    N =np.size(data)
    M1 = np.sum(data)/N #M1 is the sum of all xi divided by N
    A = 15/8*M1
    var_M1 = 1/N*(29/75-(8/15*A)**2)
    var_A = (15/8)**2*var_M1
    z_score = A/(var_A**0.5)
    p_value = scipy.stats.norm.sf(abs(z_score)) #p_value from z-score
    return A, var_A, p_value


## loop over N_values to determine N_min

def determine_N():
    p_crit = 5e-3 #critical p-value
    N=100 #startpoint for number of data points
    Nstep = 10 #iterate over N in steps of Nsteps, can be chosen for desired precision depending on computational power available.
    p_value = 1 #a start value
    while p_value>p_crit: #iterate over N until critical p_value is reached
        N += Nstep #increase #data points
        data = generate_data(N) #generate data
        A, var_A, p_value = estimator0(data) #determine p_value
    return N, A, var_A, p_value

N, A, var_A, p_value = determine_N()
print("For a p-value of 5e-3 at least N =", N, " data points are needed.")
print("Using this \n A =", A, " and \n var_A =", var_A, "are determined. \n")
#Executing determine_N a few times shows it varies from a few hundred to 1000+. This is due to the fact that the numbers are generated randomly so sometimes more favorable nummers are generated resulting in lower variances and hence lower p-values. For iterations which result in relatively low N, the generater could have been lucky in quickly generating a favorable array with low variance.


## Estimation for the data
os.chdir('/home/bram/Documents/data-analysis') #path to file
#os.chdir("C:\\Users\\Bram Ruiter\\Downloads")
data = np.genfromtxt('asymmetry.dat') #load the data

A1, var_A1, p_value1 = estimator_background(data)
print("Measured from data: \n A =", A1, "\n var_A =", var_A1, "\n")
sigma_A1 = var_A1**0.5


## Determine the confindence interval
#the confidence interval is determined as a central confidence interval for a renormalized normal distribution between [-1,1]

prob_interval = scipy.stats.norm.cdf(1,A1,sigma_A1) - scipy.stats.norm.cdf(-1,A1,sigma_A1) #the probability that A lies within [-1,1].
renormalization = 1/prob_interval #renormalization constant so that the probability that A lies within -[1,1] becomes 1.

CL = 0.683 #desired confidence level
prob = 0
Amin=A1
Amax=A1
Astep = (var_A1)**0.5*0.001 #step size of A: using 0.01sigma as a step size

while prob<CL: #iterate over A in steps of Astep untill the desired confidence level is reached
    Amin -= Astep
    Amax += Astep
    prob = renormalization*(scipy.stats.norm.cdf(Amax,A1,sigma_A1) - scipy.stats.norm.cdf(Amin,A1,sigma_A1)) #renormalized probability between Amin and Amax


print("For the desired confidence level, A lies within the interval \n[", Amin, ",", Amax, "]")

