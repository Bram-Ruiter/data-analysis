import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

beta = 1/1 #1\tau
N=100
N_experiments = 100
nbins = 5
start = 0
stop = 5

def binned(x, nbins, start=-20, stop=20): #bins the data
    bins = np.linspace(start,stop,nbins+1) #array of nbins+1 elements meaning it contains nbins
    binned = np.array([])
    for i in range(N):
        for j in range(nbins):
            if x[i]>bins[j] and x[i]<bins[j+1]:
                binned = np.append((bins[j]+bins[j+1])/2,binned)
                break
    return binned

def bin_count(x, nbins, start=-20, stop=20): #bins the data
    bins = np.linspace(start,stop,nbins+1) #array of nbins+1 elements meaning it contains nbins
    bin_count = np.zeros(nbins, dtype=int)
    for i in range(N):
        for j in range(nbins):
            if x[i]>bins[j] and x[i]<bins[j+1]:
                bin_count[j] += 1
                break
    return bin_count

#generate data
def generate_data(N):
    generated = np.random.exponential(beta,N)
    #exponential_array[i] = generated
    binned_generated = binned(generated,nbins,start,stop)
    #binned_array[i] = binned_generated
    binned_count = bin_count(generated,nbins,start,stop)
    return generated, binned_generated, binned_count

bins = np.linspace(start,stop,nbins+1)

def likelihood(tau,x): #returns negative likelihood
    L = 1
    for j in range(nbins):
        prob_interval = np.exp(-bins[j]*tau)-np.exp(-bins[j+1]*tau)
        L *= (prob_interval*N)**(x[j])*np.exp(-prob_interval*N)*1/np.math.factorial(x[j])
    return -L

def chi_squared(tau,x): #returns chi_squared
    chi2 = 0
    for i in range(nbins):
        prob_interval = np.exp(-bins[i]*tau)-np.exp(-bins[i+1]*tau)
        chi2 += (x[i]-prob_interval*N)**2
    return chi2/tau

res_l = np.array([])
res_chi = np.array([])
tau_l = 0
tau_chi = 0

#estimate tau for N_experiments and append to create sampling distribution
for i in range(N_experiments):
    data, binned_data, binned_count = generate_data(N)
    res_l = np.append(res_l, minimize(likelihood, x0=0.90, args=binned_count).x[0])
    res_chi = np.append(res_chi, minimize(chi_squared, x0=0.90, args=binned_count).x[0])
    tau_l += res_l[i]/N_experiments
    tau_chi += res_chi[i]/N_experiments

print('From likelihood, tau = ', tau_l)
print('From Chi squared, tau = ', tau_chi)

plt.figure(1)
plt.hist(res_l, label="Maximum likelihood method")
plt.hist(res_chi, label="Chi-squared method", fill=False)
plt.legend()
plt.show()




