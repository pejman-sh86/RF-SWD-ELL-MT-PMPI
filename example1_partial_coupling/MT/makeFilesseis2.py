
import numpy as np

I_RF = 1
I_SWD = 0
I_ELL = 0

###########  SWD  ###################
if I_SWD == 1:
    #log10T = np.arange(-4.0, 4.01, 0.01)
    #log10T = np.arange(-8.0, 2.1, 0.1)
    #log10T = np.arange(-4.0,5.0)
    #log10T = np.linspace(25., 250., 21)
    log10F = np.linspace(np.log10(1./100.), np.log10(1./8.),23)
    log10T = -log10F[-1::-1]

    #SWDdata = np.array([[10.0**log10fre, 9.909e3]]) # [period in seconds, veocity in km/s]
    SWDdata = np.zeros((log10T.size,2))
    SWDdata[:,0] = 10.0**log10T
    SWDdata[:,1] = 9.909e3
    np.savetxt('example_SWD.dat', SWDdata, fmt='%.10e')

############  ELL  #################
if I_ELL == 1:
    #log10T = np.arange(-4.0, 4.01, 0.01)
    #log10T = np.arange(-8.0, 2.1, 0.1)
    #log10T = np.arange(-4.0,5.0)
    log10F = np.linspace(np.log10(1./24.), np.log10(1./8.),30)
    log10T = -log10F
    #log10T = np.linspace(np.log10(0.1), np.log10(200.),50)

    #ELLdata = np.array([[10.0**log10fre, 9.909e3]]) # [period in seconds, ellipticity]
    ELLdata = np.zeros((log10T.size,2))
    ELLdata[:,0] = 10.0**log10T[-1::-1]
    ELLdata[:,1] = 9.909e3
    np.savetxt('example_ELL.dat', ELLdata, fmt='%.10e')

############  RF   #################
from obspy.signal.invsim import cosine_taper
if I_RF == 1:
    NTIME = 300
    percentage = 0.05
    time = np.arange(0., float(NTIME))
    taper = cosine_taper(npts=NTIME, p=percentage)

    RFdata = np.zeros((2,time.size))
    RFdata[0,:] = 1.
    #RFdata[1,:] = 1.
    RFdata[1,:] = taper
    np.savetxt('HON_RF.txt', RFdata, fmt = '%.10e')

