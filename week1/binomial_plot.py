import scipy, scipy.stats
x = scipy.linspace(-10,10,100)
pmf = scipy.stats.binom.pmf(x,10,0.5)
import pylab
pylab.plot(x,pmf)
pylab.show()