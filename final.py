import numpy as np
import matplotlib.pyplot as plt

##generate data
N=1000
a=0
b=2*np.pi
Afb=0.015
data = np.genfromtxt('asymmetry.dat', delimiter='/n')

x=np.random.uniform(a,b,N)

def dsigma(x):
    return 3/8*(1+(np.cos(x))**2)+Afb*np.cos(x)



y = dsigma(x)
nbins=50
#count, bins, ignored = plt.hist(y, bins=nbins, histtype='bar', label='Sampled distribution')
binwidth = (bins[nbins]-bins[0])/nbins

xmin=0
#xmax=np.max(y)
xmax=b
u = np.linspace(xmin,xmax,1000)
v=dsigma(u)*(N*binwidth)
plt.plot(u,v)

plt.show()



## estimator part
def estimator0(data): #estimate A for a given data set assuming no background
    N_dat =np.size(data)
    M1 = np.sum(data)/N_dat #M1 is the sum of all xi divided by N
    A = 3/2*M1
    var_M1 = 1/N*(2/5-4/9*A**2)
    var_A = (3/2)**2*var_M1
    sigma_A = (var_A)**0.5
    return A, sigma_A

def estimator_background(data): #estimate A for a given data set assuming background
    N = np.size(data)
    return




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