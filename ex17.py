## Exercise 17
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

angles = np.genfromtxt('muon.dat', delimiter='/n')
data = np.cos(angles)
N=np.size(angles)

# xi_hat = 3M1_hat
# var(xi_hat) = 1/N*(3-xi**2)


y = np.linspace(-np.pi,np.pi,100)

def g(x, xi): #probability function
    return (1+xi*x)/2


def l(xi, data): #negative log likelyhood function, which is to be maximized with respect to xi
    for i in range(N):
        if i == 0:
            l = -1*np.log(g(data[i],xi))
        else:
            l -= np.log(g(data[i],xi))
    return l

res = minimize(l, x0=0, args=data)
xi = res.x[0]
print('xi = ', xi)

def var(xi, data): #see notebook for formula
    for j in range(N):
        if j==0:
            v = -data[j]**2/(1+data[j]*xi)**2
        else:
            v += -data[j]**2/(1+data[j]*xi)**2
    var = -1/v
    return var

print('sigma = ', np.sqrt(var(xi,data)))

plt.figure(1) #see if result makes sense
y = np.linspace(0.3,0.4,1000)
plt.plot(y,l(y,data))
plt.show()



