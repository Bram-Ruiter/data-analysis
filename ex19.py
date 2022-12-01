import numpy as np
import matplotlib.pyplot as plt

mu=0
sigma=5
N = 100 #sample size
normal_array = np.random.normal(mu,sigma,N)

def gaussian(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)

def mu(x,N): #mean estimator
    mu = 0
    for i in range(N):
        mu += x[i]
    return mu/N

def var(x,mu,N): #variance estimator for unkown mean
    var = 0
    for i in range(N):
        var += (x[i]-mu)**2
    return var/(N-1)

def binned(x, nbins, start=-20, stop=20): #bins the data
    bins = np.linspace(start,stop,nbins)
    binned = np.array([])
    for i in range(N):
        for j in range(nbins):
            if x[i]>bins[j] and x[i]<bins[j+1]:
                binned = np.append(bins[j],binned)
                break
    return binned




#unbinnned
mu_sample = mu(normal_array,N)
sigma_sample = np.sqrt(var(normal_array,mu_sample,N))
print('For unbinned data: \n mu = ', mu_sample, '\n sigma = ', sigma_sample )

#plot normally distributed numbers
plt.figure(1)
count, bins, ignored = plt.hist(normal_array, label = 'normally distributed numbers')
binwidth = (bins[1]-bins[0])

#bin 1
nbins1=100
binned1 = binned(normal_array, nbins1)
label1 = 'binned data with nbins = '+ str(nbins1)
mu_sample1 = mu(binned1,N)
sigma_sample1 = np.sqrt(var(binned1,mu_sample1,N))
print('For binned data with nbins = ',nbins1,': \n mu = ', mu_sample1, '\n sigma = ', sigma_sample1 )

#plot the binned function
plt.hist(binned1, fill = False, label = label1)

#plot a gaussian fit
array = np.linspace(-20,20,1000)
plt.plot(array, N*binwidth*gaussian(array,mu_sample,sigma_sample), label ='gaussian fit')
plt.legend()
plt.show()