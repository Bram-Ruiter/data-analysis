##exercise 5
import numpy as np
import matplotlib.pyplot as plt


#plot y=-ln(x) with x random distributed nubers

N = 2000
a=0
b=1
x = np.random.uniform(a,b,N)
y=-1*np.log(x)

nbins=50
count, bins, ignored = plt.hist(y, bins=nbins, histtype='bar', label='Sampled distribution')
binwidth = (bins[nbins]-bins[0])/nbins
xmin=0
xmax=np.max(y)
u=np.linspace(xmin,xmax,1000) #x values
v = np.exp(-1*u)*N*binwidth #y values following exponential, need to get it to the same scale as sampled distrubution so *N *binwidth=amount of samples in a bin
plt.plot(u,v, label='Exponential distribution')
plt.legend()

plt.show()


## Exercise 6
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial

xmin,xmax = 1,30
x = np.linspace(xmin,xmax,100)
N=4000
nbins=10


def p(x, mu):
    return np.exp(-1*mu)*(mu**(x))*1/factorial(x)

#mu=5
plt.figure(1)
count, bins, ignored = plt.hist(np.random.poisson(5,N), label='mu=5') #random plot
binwidth = (bins[nbins]-bins[0])/nbins
plt.plot(x,p(x,5)*N*binwidth, label='mu = 5') #fit

#mu=10
plt.figure(2)
count, bins, ignored = plt.hist(np.random.poisson(10,N), label='mu=10')
binwidth = (bins[nbins]-bins[0])/nbins
plt.plot(x,p(x,10)*N*binwidth, label='mu=10')

#mu=15
plt.figure(3)
count, bins, ignored = plt.hist(np.random.poisson(15,N), label='mu=15') #random plot
binwidth = (bins[nbins]-bins[0])/nbins
plt.hist(np.random.poisson(10,N) + np.random.poisson(5,N), fill=False, label='combined plot') #sum of random variables plot
#plt.plot(x,(p(x,5)+p(x,10))*N*binwidth, label='combined plot')



#


#plt.plot(x,p(x,5)+p(x,10), label='combined plot')
plt.legend()
plt.show()



































##example
#! /usr/bin/env python

# pull in the code from numpy
import numpy as np
from scipy.stats import norm

# specify the mean and width of the (normal, see below) distribution
mu, sig = 10, 2
# specify the number of samples to be taken from this distribution
nsamples = 1000

# draw 10000 samples from a normal distribution (do e.g. 'help(np.random)' in an interactive python to get help)
s = np.random.normal(mu, sig, nsamples)
xmin, xmax, ymin, ymax = plt.axis()

# pull in the plotting code from matplotlib
import matplotlib.pyplot as plt

# XKCD-like appearance
plt.xkcd()

# after the above import, doing a 'help(plt.hist)' in an interactive python session can be used to see more detail
# note that in older matplotlib versions, the 'normed = False' argument was needed to prevent the plot from being normalised to 1.
nbins = 50
#count, bins, ignored = plt.hist(s, bins=nbins, normed=False, histtype='bar')
count, bins, ignored = plt.hist(s, bins=nbins, histtype='bar')

# the adjustment here is useful primarily for the bottom margin: the x-axis label tends to be cut off if the font size is large
plt.subplots_adjust(wspace = 0.25, hspace = 0.25, top = 0.90, bottom = 0.16, left = 0.14, right = 0.96)

# since we didn't print the distribution normalised to unity, we need to compute the bin width
binwidth = (bins[nbins]-bins[0])/nbins

# also compute the expected curve (use the whole plot range)
xmin, xmax, ymin, ymax = plt.axis()
xg = np.linspace(xmin, xmax, 200)
plt.plot(xg, binwidth * nsamples * norm(loc = mu, scale = sig).pdf(xg), linewidth = 2, color = 'r')

# add some labels
plt.ylabel('number of entries / {:5.2f}'.format(binwidth), size = 'large')
plt.xlabel('a random variable', size = 'large')

# show the results
plt.show()

##ex2
#! /usr/bin/env python

# pull in the code from numpy and scipy
import numpy as np
from scipy.stats import poisson, binom, norm, chi2
xmin, xmax, ymin, ymax = plt.axis()

# pull in the plotting code from matplotlib
import matplotlib.pyplot as plt

#plt.xkcd()

# subsivide the canvas into four plots and immediately "unpack" the results
fig, ax = plt.subplots(2,2, figsize = (10.0, 6.0)) # figure size in inches
fig.suptitle('Example distributions')
axbinom, axpois, axnorm, axchi2 = ax[0][0], ax[0][1], ax[1][0], ax[1][1]

# set the 'wspace' and 'hspace' elements to create a bit more space for axis labels
plt.subplots_adjust(wspace = 0.25, hspace = 0.25, top = 0.92, bottom = 0.08, left = 0.08, right = 0.96)

# in older versions of matplotlib, the 'steps' linestyle offers an alternative to the use of the 'steps' function
# (as is appropriate to display discete probability functions)

# binomial
xb  = np.arange(0,21)
p = 0.3
line_bn0,   = axbinom.step(xb, binom(20,p).pmf(xb), 'b', lw = 2, label = "N=20, p=0.3")
line_bn1,   = axbinom.step(xb, binom(10,p).pmf(xb), 'r', lw = 2, label = "N=10, p=0.3")
#line_bn0,   = axbinom.plot(xb, binom(20,p).pmf(xb), 'b', lw = 2, ls = 'steps', label = "N=20, p=0.3")
#line_bn1,   = axbinom.plot(xb, binom(10,p).pmf(xb), 'r', lw = 2, ls = 'steps', label = "N=10, p=0.3")
axbinom.legend(title = 'binomial distributions', handles = [ line_bn1, line_bn0 ])
axbinom.set_xlabel(r'$n$')
axbinom.set_ylabel('probability')

# Poisson. Adjust the range so that at least 0.999 of the "widest" distribution fits.
xp = np.arange(0,poisson(6).ppf(0.999))
line_ps0,   = axpois.step(xp, poisson(6).pmf(xp), 'b', lw = 2, label = r"$\mu=6$")
line_ps1,   = axpois.step(xp, poisson(3).pmf(xp), 'r', lw = 2, label = r"$\mu=3$")
#line_ps0,   = axpois.plot(xp, poisson(6).pmf(xp), 'b', lw = 2, ls = 'steps', label = r"$\mu=6$")
#line_ps1,   = axpois.plot(xp, poisson(3).pmf(xp), 'r', lw = 2, ls = 'steps', label = r"$\mu=3$")
axpois.legend(handles = [ line_ps1, line_ps0 ])
axpois.legend(title = 'Poisson distributions', handles = [ line_ps1, line_ps0 ])
axpois.set_xlabel(r'$n$')
axpois.set_ylabel('probability')

# normal. Again adjust the range
n = norm(scale = 3)
xn = np.linspace(n.isf(0.999), n.ppf(0.999), 200)
line_nr0,   = axnorm.plot(xn, n.pdf(xn), 'b', lw = 2, label = r"$\sigma = 3$")
line_nr1,   = axnorm.plot(xn, norm(scale = 1).pdf(xn), 'r', lw = 2, label = r"$\sigma = 1$")
axnorm.legend(title = 'normal distributions', handles = [ line_nr1, line_nr0 ])
axnorm.set_xlabel(r'$x$')
axnorm.set_ylabel('probability densit') # this placement is usually felt to be preferable

# chi2. Again adjust the range
xc = np.linspace(0, chi2(20).ppf(0.999), 200)
line_ch0,   = axchi2.plot(xc, chi2(20).pdf(xc), 'b', lw = 2, label = "N=20")
line_ch1,   = axchi2.plot(xc, chi2(10).pdf(xc), 'r', lw = 2, label = "N=10")
axchi2.legend(title = r'$\chi^{2}$ distributions', handles = [ line_ch1, line_ch0 ])
axchi2.set_xlabel(r'$x$')
axchi2.set_ylabel('probability density')#, loc = 'top')

# show everything
plt.show()