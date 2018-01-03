
from pylab import *

res = np.loadtxt('../RadiusResults.txt')
res = array(res)

plot(res[0],res[1],label='no sampling')
#plot(res[2],res[3],label='with sampling')
#legend(loc='upper left')

xlabel('time')
ylabel('R')
title('Maximal radius of the covariance ellipsoid $x^{T}\Sigma^{-1} x = const$')

savefig('../radiusPlot.png')

show()