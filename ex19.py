import numpy as np
import matplotlib.pyplot as plt

mean=0
std=5
N = 200 #sample size
Nsamples = 200 #number of experiments

#normal_array = np.random.normal(mean,std,N)

def gaussian(x,mu,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)

def mu(x,N): #mean estimator
    mu = 0
    for i in range(N):
        mu += x[i]
    return mu/N

def var_est(x,mu,N): #variance estimator for unkown mean
    var = 0
    for i in range(N):
        var += (x[i]-mu)**2
    return var/(N-1)


def binned(x, nbins, start=-20, stop=20): #bins the data
    bins = np.linspace(start,stop,nbins+1) #array of nbins+1 elements meaning it contains nbins
    binned = np.array([])
    for i in range(N):
        for j in range(nbins):
            if x[i]>bins[j] and x[i]<bins[j+1]:
                binned = np.append((bins[j]+bins[j+1])/2,binned)
                break
    return binned



#declare unbinnned
mu_sample = np.array([])
sigma_sample = np.array([])
mu_sample_tot = 0
sigma_sample_tot = 0
mu_sample1_tot = 0
sigma_sample1_tot = 0


#declare binned
mu_sample1 = np.array([])
sigma_sample1 = np.array([])

#loop over all experiments and calculate mu,sigma
for i in range(Nsamples):
    #for unbinned data
    normal_array = np.random.normal(mean,std,N)
    mu_sample = np.append(mu_sample,mu(normal_array,N))
    sigma_sample = np.append(sigma_sample,np.sqrt(var_est(normal_array,mu_sample[i],N)))
    #for binned data
    nbins1=100
    binned1 = binned(normal_array, nbins1)
    mu_sample1 = np.append(mu_sample1, mu(binned1,N-1))
    sigma_sample1 = np.append(sigma_sample1, np.sqrt(var_est(binned1,mu_sample1[i],N-1)))

for i in range(Nsamples):
    mu_sample_tot += mu_sample[i]
    sigma_sample_tot += sigma_sample[i]
    mu_sample1_tot += mu_sample1[i]
    sigma_sample1_tot += sigma_sample1[i]

#now take averages
#for sigma, devide the variance of the population by the sample size to get the variance belonging to the sampling distribution
mu_sample_avg = mu_sample_tot/Nsamples
sigma_sample_avg = sigma_sample_tot/Nsamples*1/np.sqrt(N)
mu_sample1_avg = mu_sample1_tot/Nsamples
sigma_sample1_avg = sigma_sample1_tot/Nsamples*1/np.sqrt(N)

#print sigma
print('For unbinned data: \n mu = ', mu_sample_avg, '\n sigma = ', sigma_sample_avg )
print('For binned data with nbins = ',nbins1,': \n mu = ', mu_sample1_avg, '\n sigma = ', sigma_sample1_avg )
#print('sigma = ', sigma_sample1_avg)


#plot data for mu
plt.figure(1)
counts,bins,ignored = plt.hist(mu_sample, bins=10, label='unbinned')
binwidth = (bins[1]-bins[0])
array = np.linspace(-5,5,1000)
plt.plot(array, Nsamples*binwidth*gaussian(array,mu_sample_avg,sigma_sample_avg), label ='gaussian fit for unbinned')
counts1,bins1,ignored = plt.hist(mu_sample1, bins=10, fill = False, label='binned')
plt.plot(array, Nsamples*binwidth*gaussian(array,mu_sample1_avg,sigma_sample1_avg), label ='gaussian fit for binned')
plt.xlabel('mu')
plt.ylabel('count')
plt.legend()
plt.show()



