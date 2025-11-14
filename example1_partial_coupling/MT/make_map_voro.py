import numpy as np 

filebase = 'HON'

I_RV = 0 
I_SWD = 0
I_ELL = 0 
I_MT = 1

I_VPVS = 1 
NLMX = 30
NRF1 = 1
NMODE = 1
NMODE_ELL = 1

interfaces = np.array([0., 1., 8., 35, 170.])
interfaces = np.array([0., 1., 170.])
#interfaces = np.array([0., 1.40, 8.20, 37.00, 170.])
#interfaces = np.array([0., 1.13, 8.17, 35.12, 76.20, 175.29])
#interfaces = np.array([0., 1.33, 8.17, 14, 18, 28., 35.12, 76.20, 175.29])
Vs = np.array([2., 3.2, 3.8, 4.7, 4.7])
#Vs = np.array([2.1, 3,2593, 3.7945, 4.7244, 4.6517])
#Vs = np.array([1.76, 3.28, 3.81, 4.6666, -100., 4.6607])
#Vs = np.array([2.10, 3.28, 3.90, 3.78, 2.96, 3.6, 4.6666, -100., 4.6607])
VpVs = np.array([1.8, 1.7, -100., -100., 2.])
#VpVs = np.array([1.81, 1.72, 1.70, 1.72, 1.90, 1.89])
#VpVs = np.array([1.81, 1.72, 1.70, 1.72, 1.90, 1.89])
#VpVs = np.array([1.81, 1.72, 1.70, -100., -100., -100., 1.72, 1.90, 1.89])
vel_ref_file = 'HON_vel_ref.txt'
NVELREF = 5
log10_sigma = np.array([-1., -3., -100., -100., -2.])
log10_sigma = np.array([-1., -3., -2.])
#log10_sigma = np.array([-.75, -2.97, -100., -100., -1.68])
#log10_sigma = np.array([-.67, -2.88, -100., -100., -100., -1.64])
#log10_sigma = np.array([-.67, -2.88, -100., -100., -100., -100., -100., -100., -1.64])
sdparR = 1.e-0
sdparSWD = 1.e-0
sdparELL = 1.e-0
sdparMT = 1.e-0
###########################################################
def  getref(z):
    vel_ref = np.loadtxt(vel_ref_file,skiprows=1)
    if z >= vel_ref[NVELREF-1,0]:
        vref = vel_ref[NVELREF-1,1]
        vpvsref = vel_ref[NVELREF-1,2]
    elif z == 0.:
        vref = vel_ref[0,1]
        vpvsref = vel_ref[0,2]
    else:
        iint = 0
        for ipar in range(NVELREF):
            if (z-vel_ref[ipar,0]) <= 0.:
                break
            iint +=  1
        iint -= 1 
        grad = (vel_ref[iint+1,1]-vel_ref[iint,1])/(vel_ref[iint+1,0]-vel_ref[iint,0])
        dz = (z-vel_ref[iint,0])
        vref = vel_ref[iint,1] + dz*grad
        grad = (vel_ref[iint+1,2]-vel_ref[iint,2])/(vel_ref[iint+1,0]-vel_ref[iint,0])
        vpvsref = vel_ref[iint,2] + dz*grad
    return vref, vpvsref

vref, vpvsref  = np.zeros(len(Vs)), np.zeros(len(VpVs))
dVs, dVpVs = np.zeros(len(Vs)), np.zeros(len(VpVs))
for layeridx in range(len(interfaces)):
    vref[layeridx], vpvsref[layeridx] = getref(interfaces[layeridx])
    if Vs[layeridx]>-99.: 
        dVs[layeridx] = Vs[layeridx] - vref[layeridx]
    else:
        dVs[layeridx] = -100.
    if VpVs[layeridx]>-99.:
        dVpVs[layeridx] = VpVs[layeridx] - vpvsref[layeridx]
    else:
        dVpVs[layeridx] = -100.

##########################################
NPL = 1
if (I_RV==-1 or I_SWD==1 or I_ELL==1):
    NPL += 1
    if (I_VPVS==1): NPL += 1
if(I_MT==1): NPL += 1

nlayers = len(interfaces)
NFPMX = NLMX*NPL
ncount = NFPMX + 1 + 3*NRF1 + 3*NRF1 + 2*NMODE + 2*NMODE_ELL + 1

map_voro = np.zeros((1,ncount))

istart = 0
map_voro[0, istart] = float(nlayers)
istart += 1
map_voro[0, istart:NPL*nlayers+1:NPL] = interfaces

if (I_RV==-1 or I_SWD==1 or I_ELL==1): 
    istart += 1
    map_voro[0, istart:NPL*nlayers+1:NPL] = dVs
    if (I_VPVS==1):
        istart += 1 
        map_voro[0, istart:NPL*nlayers+1:NPL] = dVpVs

    if(I_RV==-1): map_voro[0, NFPMX+1:NFPMX+1+NRF1] = sdparR
    if(I_SWD==1): map_voro[0, NFPMX+1+3*NRF1:NFPMX+1+3*NRF1+NMODE] = sdparSWD
    if(I_ELL==1): map_voro[0, NFPMX+1+3*NRF1+NMODE:NFPMX+1+3*NRF1+NMODE+NMODE_ELL] = sdparELL

if (I_MT==1):
    istart += 1
    map_voro[0, istart:NPL*nlayers+1:NPL] = log10_sigma
    map_voro[0, NFPMX+1+3*NRF1+NMODE+NMODE_ELL] = sdparMT

np.savetxt(filebase+'_map_voro_true.dat', map_voro, fmt='%.10e')
