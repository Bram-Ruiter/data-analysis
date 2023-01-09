##import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

##generate data
N=50000
Afb=0.015
#data = np.genfromtxt('asymmetry.dat', delimiter='/n')


def pdf(x):
    A=Afb
    return 3/8*(1+x**2)+A*x

def inverse(x):
    return 1/3*(7*np.pi/(-10*(np.pi)**3+54*(np.pi)**2*x+5.1962*(-(np.pi)**4*(9*(np.pi)**2+40*(np.pi)*x-108*x**2))**0.5)**(1/3) + (-10*(np.pi)**3+54*(np.pi)**2*x+5.1962*(-(np.pi)**4*(9*(np.pi)**2+40*(np.pi)*x-108*x**2))**0.5)**(1/3)/(np.pi) -4)


def generate_data(N): #Generates data on the domain -1,1 following the pdf as follows from the differential cross section
    x=np.random.uniform(0,1,N)
    A=Afb
    y = np.array([])
    for i in range(N):
        p = [1,4*A,3,4-4*A-8*x[i]]
        for root in np.roots(p): #there are only real roots in our domain, so select for this.
            if abs(root.imag)<1e-5: #if the root is not complex, append the real part to the generated array
                y = np.append(y,root.real)
    return y


#some rules of code to see if the generated data indeed follows the pdf
y = generate_data(N)
nbins=50
count, bins, ignored = plt.hist(y, bins=nbins, histtype='bar', label='Sampled distribution')
binwidth = (bins[nbins]-bins[0])/nbins

#plot the pdf
xmin=-1
xmax=1
u = np.linspace(xmin,xmax,1000)
v=pdf(u)*(N*binwidth)
plt.plot(u,v)

plt.show()



## estimator part
def estimator0(data): #estimate A for a given data set assuming no background
    N_dat =np.size(data)
    M1 = np.sum(data)/N_dat #M1 is the sum of all xi divided by N
    A = 3/2*M1
    var_M1 = 1/N*(2/5-4/9*A**2)
    var_A = (3/2)**2*var_M1
    return A, var_A



def estimator_background(data): #estimate A for a given data set assuming background
    N_dat =np.size(data)
    M1 = np.sum(data)/N_dat #M1 is the sum of all xi divided by N
    A = 15/8*M1-3/8
    var_M1 = 1/N*(8/25-(8/15*A+1/5)**2)
    var_A = (15/8)**2*var_M1
    return A, var_A


## p value
def p_value(A, var_A): #determine the p_value from the z-score using that A is assymptotically normally distributed.
    z_score = A/(var_A**0.5)
    p_value = scipy.stats.norm.sf(abs(z_score))
    return p_value

## loop over p_values

def determine_N():
    p_crit = 5*10**(-3)
    N=100 #start at 100 data points and iterate in steps of 100
    p_value = 1
    while p_value<p_crit:
        data = generate_data(N)
        A, var_A = estimator0(data)
        p_value = p_value(A, var_A)
        N += 100
    return N, A, var_A, p_value



## estimator with background




## Confidence interval








##
N=10000
a=0*np.pi
b=2*np.pi
x=np.append(np.random.uniform(a,b,N),-np.random.uniform(a,b,N))
y=np.arcsin(x/(b-a))
#y=1/(b-a)*(1-x**2)**(-0.5)

nbins=50
count, bins, ignored = plt.hist(y, bins=nbins, histtype='bar', label='Sampled distribution')
binwidth = (bins[nbins]-bins[0])/nbins


xmin=-np.pi
xmax=np.pi
u=np.linspace(xmin,xmax,1000) #x values

v = np.cos(u)*N*binwidth #y values following exponential, need to get it to the same scale as sampled distrubution so *N *binwidth=amount of samples in a bin
plt.plot(u,v, label='Cosine distribution')
plt.legend()

plt.show()