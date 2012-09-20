# Bayesian and frequentist interpretations of a Bernoulli process

# Cyrille Rossant's blog post: Introduction to Bayesian thinking
# URL - http://cyrille.rossant.net/introduction-to-bayesian-thinking/

# This Python script simulates 10000 independent Bernoulli processes 
# with 1000 trials each. You need Numpy, Scipy and Matplotlib to 
# execute it. It saves 1000 PNG files in a png folder. A video is 
# made out of those images using VirtualDub.

from numpy import *
from numpy.random import *
from scipy import *
from scipy.stats import *
from pylab import *
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size': 22})
 
# posterior distribution (particular case of beta distribution)
posterior = lambda n,m,p: (n+1)*binom(n, p).pmf(m)
 
# length of the Bernoulli processes
nmax = 1000
 
# number of independent sample paths, for the frequentist simulation
samples = 10000
 
# Bernoulli parameter
p = .35
 
ps = linspace(0., 1., 10000)
 
# stochastic experiment
X = rand(nmax, samples)<=p
 
# sample paths of Bernoulli processes are contained in the columns of Y
Y = cumsum(X, axis=0)
 
# create a figure
fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
 
bins = linspace(0., 1., 50)
for n in xrange(0, nmax):
    print "%d/%d" % (n, nmax)
    ax.cla()
    ax2.cla()
 
    # Bayesian interpretation
    subplot(121)
    post = posterior(n*ones(ps.shape), Y[n,0], ps)
    ax.plot(ps, post)
    if X[n,0]==False:
        c = 'r'
    else:
        c = 'g'
    axvline(X[n,0], color=c, lw=10)
    axvline(p, ls='--', color='k')
    ylim(0, 30)
    xticks([0,p,1])
    yticks([])
    title("Bayesian")
 
    # Frequentist interpretation
    subplot(122)
    hist(Y[n,:]/float(n), bins=bins)
    axvline(p, ls='--', color='k')
    ylim(0, 30*float(samples)/len(bins))
    xticks([0,p,1])
    yticks([])
    title("Frequentist")
    fig.savefig("png/fig%03d.png" % n)