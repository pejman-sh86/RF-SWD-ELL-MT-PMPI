
import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt('HON_SWD_50points.dat')[:,0]

logw = np.r_[  -1.*np.ones(4),    np.linspace(np.log(.1)/np.log(50.),np.log(.01)/np.log(50.),3),   np.linspace(np.log(.01)/np.log(20.),np.log(.001)/np.log(20.),43)  ]

w = np.r_[   10.**logw[0:4],        50.**logw[4:7],         20.**logw[7:]  ]

plt.plot(p[0:4],w[0:4],        p[4:7],w[4:7],         p[7:],w[7:])
plt.show()
