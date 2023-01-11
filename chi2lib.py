## estimation using chi squared
def cdf(x, A=Afb): #the pdf integrated to a point x
    return 3/8*(x+1/3*x**3)+1/2*A*x**2 - 3/8*(-1+1/3*-1)-1/2*A*1

def cdf_background(x, A=Afb): #the cdf of background pdf
    return cdf(x,A)*0.8 +1/2*(x-1)*0.2

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


## some other functions
def p_value(A, var_A):
    z_score = A/(var_A**0.5)
    p_value = scipy.stats.norm.sf(abs(z_score))
    return p_value

def gaussian(x, mu,var):
    return 1/(var*2*np.pi)**0.5*np.exp(-0.5*((x-mu)**2/var))