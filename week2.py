##Exercise 11
import numpy as np
import matplotlib.pyplot as plt
N=10000
n=1000
nbins=20

avg_poisson = np.array([])
avg_norm = np.array([])
x=np.linspace(9,11,1000)

def gaus(mu,sigma,x):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)

for i in range(n): #sample n times in total
    #sample N numbers from poisson distribution and average
    random_poisson = np.random.poisson(10,N)
    avgp = np.sum(random_poisson)/N
    avg_poisson = np.append(avg_poisson, avgp)
    #Same for gaussian
    random_norm = np.random.normal(10, np.sqrt(10), size = N)
    avgn = np.sum(random_norm)/N
    avg_norm = np.append(avg_poisson, avgn)

#plot sampling distributions
count, bins, ignored = plt.hist(avg_poisson, nbins, label='poisson average')
binwidth = (bins[1]-bins[0])
plt.hist(avg_norm, nbins, label='gaussian average', fill=False)

#gaussian fit with correct mu and sigma
plt.plot(x,n*binwidth*gaus(10,np.sqrt(10/N),x))

plt.xlim(9.5,10.5)
plt.legend()
plt.show()

#plot shows that both means are almost identically normally distributed

##Exercise 12
import numpy as np
N = 10**7

random2d = np.random.uniform(0,1,[2,N]) #generate 2*N random numbers

within =0
for i in range(N):
    if ((random2d[0,i])**2 + (random2d[1,i])**2) > 1:
        pass
    else:
        within +=1

#In total N pair of numbers generated. within/N gives area of quarter circle = pi/4
pi = within/N*4
print('pi =', pi)
print('difference = ', np.pi-pi)

#Darts are bernouli distributed with chance pi/4:
#Variance for pi/4 is: 1/n*pi/4*(1-pi/4)
#Variance for pi: 1/n*4pi*(1-pi/4) -> sigma = 0.0005 which seems to agree with observation

##Exercise 13
import numpy as np
import os
angles = np.genfromtxt('/home/bram/Documents/data_analysis/muon.dat', delimiter='/n')
N=np.size(angles)

# E(t)=0 so  need to use second moment, t=theta
# E(t^2)= pi^2/3 - 2*xi = M2 (t2-bar)
# xi = (M2-pi^2/3)/2
#var(M2) = (M4-M2^2)/N
#M4 = pi^4/5-4(pi^2-6)*xi

x = np.linspace(-np.pi,np.pi,100)

def g(cosx, xi):
    return (1+xi*cosx)/2

#E(xi) = 1/12*t**2(2*cosx*t+3)
# xi=(2*g-1)/cosx

M2_hat = 0
for i in range(N):
    M2_hat += 1/N*angles[i]**2

xi = (M2_hat-(np.pi**(2))/3)/2
print('xi =', xi)

M2 = (np.pi**2)/3-2*xi
M4 = (np.pi**4)/5-4*((np.pi**2)-6)*xi
var = (M4-M2**2)/N
sigma = np.sqrt(var)
print('sigma =', sigma)
''
count, bins, ignored=plt.hist(angles, fill=False)
binwidth = (bins[1]-bins[0])
plt.plot(x,g(np.cos(x),xi)*N*binwidth)
plt.show()

