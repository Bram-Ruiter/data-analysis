##import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os
from scipy.optimize import minimize as minimize
from scipy.optimize import minimize_scalar as minimize_scalar

##generate data
N=10000
Afb=0.015
#data = np.genfromtxt('asymmetry.dat', delimiter='/n')


def pdf(x, A=Afb): #the pdf which follows from the formula for the differential cross section
    return 3/8*(1+x**2)+A*x

def cdf(x, A=Afb): #the pdf integrated to a point x
    return 3/8*(x+1/3*x**3)+1/2*A*x**2 - 3/8*(-1+1/3*-1)-1/2*A*1

def pdf_background(x, A=Afb): #pdf corrected for 20% background
    return pdf(x,A)*0.8+(1/2)*0.2

def cdf_background(x, A=Afb): #the cdf of background pdf
    return cdf(x,A)*0.8 +1/2*(x-1)*0.2

def generate_data(N): #Generates data on the domain -1,1 following the pdf
    x=np.random.uniform(0,1,N)
    A=Afb
    y = np.array([])
    for i in range(N):
        p = [1,4*A,3,4-4*A-8*x[i]]
        for root in np.roots(p): #there are only real roots in our domain, so select for this.
            if abs(root.imag)<1e-5: #if the root is not complex, append the real part to the generated array
                y = np.append(y,root.real)
    return y


##see if the generated data indeed follows the pdf
y = generate_data(N)
nbins=50
count, bins, ignored = plt.hist(y, bins=nbins, histtype='bar', label='Sampled distribution')
binwidth = (bins[nbins]-bins[0])/nbins

#plot the pdf
xmin=-1
xmax=1
u = np.linspace(xmin,xmax,1000)
v=pdf(u)*(N*binwidth) #make the plot and account for the histogram not being normalized
plt.plot(u,v, label ='pdf')
plt.legend()
plt.xlabel('cos($\\theta$)')
plt.ylabel("Count")
plt.title("Sampling distribution")
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

def estimator_background(data): #estimate A for a given data set assuming background
    N =np.size(data)
    M1 = np.sum(data)/N #M1 is the sum of all xi divided by N
    A = 15/8*M1
    var_M1 = 1/N*(29/75-(8/15*A)**2)
    var_A = (15/8)**2*var_M1
    z_score = A/(var_A**0.5)
    p_value = scipy.stats.norm.sf(abs(z_score))
    return A, var_A, p_value


def p_value(A, var_A):
    z_score = A/(var_A**0.5)
    p_value = scipy.stats.norm.sf(abs(z_score))
    return p_value


## estimation using chi squared
def binned(x, nbins, start=-1, stop=1): #bins the data in interval between start and stop
    N = np.size(x)
    bins = np.linspace(start,stop,nbins+1) #array containing the bin edges
    binned = np.zeros(nbins)
    for i in range(N): #sum over all measurements
        for j in range(nbins):
            if x[i]>bins[j] and x[i]<bins[j+1]: #determine in which bin the measurement lies
                binned[j] += 1 #add 1 if the measurement is in the bin
                break
    return binned, bins


def chi2(A, data, nbins): #chi squared function which needs to be minimized
    N = np.size(data) #number of measurements
    binned_data, bins = binned(data, nbins)
    chi2 = 0
    for i in range(nbins): #sum over all bins
        x_exp = (cdf(bins[i+1],A)-cdf(bins[i],A))*N #expected number of measurements in a bin
        chi2 += (binned_data[i]-x_exp)**2/x_exp #chi2 component for the bin
    return chi2

def minimize_chi2(data, nbins=50):
    mu_chi = minimize(chi2, x0=0.015, args=(data, nbins)) #minimize chi2 with respect to A
    dof = nbins-1 #N-1 degrees of freedom since 1 unknown parameter
    p_chi = 1-scipy.stats.chi2.cdf(mu_chi.fun,dof)
    return mu_chi, p_chi

def chi2_background(A, data, nbins): #chi squared function which needs to be minimized
    N = np.size(data) #number of measurements
    binned_data, bins = binned(data, nbins)
    chi2 = 0
    for i in range(nbins): #sum over all bins
        x_exp = (cdf_background(bins[i+1],A)-cdf_background(bins[i],A))*N #expected number of measurements in a bin
        chi2 += (binned_data[i]-x_exp)**2/x_exp #chi2 component for the bin
    return chi2

def minimize_chi2_background(data, nbins=50):
    mu_chi = minimize(chi2_background, x0=0.015, args=(data, nbins)) #minimize chi2_background with respect to A
    dof = nbins-1 #N-1 degrees of freedom since 1 unknown parameter
    p_chi = 1-scipy.stats.chi2.cdf(mu_chi.fun,dof)
    return mu_chi, p_chi


## loop over N_values to determine N_min

def determine_N():
    p_crit = 5e-3 #critical p-value
    N=0 #startpoint for number of data points
    Nstep = 10 #iterate over N in steps of Nsteps, can be chosen for desired precision depending on computational power available.
    p_value = 1
    while p_value>p_crit: #iterate over N until critical p_value is reached
        N += Nstep
        data = generate_data(N)
        A, var_A, p_value = estimator0(data)

    return N, A, var_A, p_value

N, A, var_A, p_value = determine_N()
print("For a p-value of 5e-3 at least N =", N, " data points are needed.")
print("Using this A =", A, " and sigma_A =", var_A**0.5, "are determined.")
#Executing determine_N a few times shows it varies around N=1000. This is due to the fact that the numbers are generated randomly so sometimes more favorable nummers are generated resulting in lower variances and hence lower p-values.


## Estimation for the data
#os.chdir('/home/bram/Documents/data-analysis')
os.chdir("C:\\Users\\Bram Ruiter\\Downloads")
data = np.genfromtxt('asymmetry.dat') #load the data

A1, var_A1, p_value1 = estimator_background(data)
print("Measured from data: \n A =", A1, "\n var_A =", var_A1)
sigma_A1 = var_A1**0.5


## Determine the confindence interval
def gaussian(x, mu,var):
    return 1/(var*2*np.pi)**0.5*np.exp(-0.5*((x-mu)**2/var))


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
    prob = renormalization*(scipy.stats.norm.cdf(Amax,A1,sigma_A1) - scipy.stats.norm.cdf(Amin,A1,sigma_A1))


print("For the desired confidence level, A lies within the interval \n[", Amin, ",", Amax, "]")

