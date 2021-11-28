import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad
import math
from scipy.special import expit

q=1.6*10**(-19) #unit: C
m0=9.11*10**(-31) #unit:  kg
h = 6.626*10**(-34) #unit: J s
h_bar = 1.055*10**(-34) #unit: J s
m_si = 0.19*m0 #For Si, effective mass is 0.19*m0.
temperature = 300# unit:K
bolzcons = 1.38*10**(-23) #unit: JK-1
di_cons_0 = 8.854*10**(-14) #unit: F/cm
di_cons_sio2 = 3.9*di_cons_0
di_cons_si = 11.7*di_cons_0
mobility_Si = 1000 #unit: cm2/(Vs)
t_ox = 1.1 # unit: nm
width = 100 # channel width; unit: nm
L_c = 40 # channel length; unit: nm
V_t = 0.05 # unit: V
NA = 2.7*10**18 #unit: cm-3
NI = 10**10 #unit: cm-3
#square law model
def big_F(psi_s):
    return np.sqrt((np.exp(-q*psi_s/bolzcons/temperature)+q*psi_s/bolzcons/temperature-1)+((NI/NA)**2)*(np.exp(q*psi_s/bolzcons/temperature)+q*psi_s/bolzcons/temperature-1))

def charge_surface(psi_s):
    return np.sqrt(2*di_cons_si*bolzcons*temperature*NA)*big_F(psi_s)

y=[]

psi_s_range=np.linspace(-1,3)
for i in psi_s_range:
    y.append(charge_surface(i))

plt.plot(psi_s_range,y)

plt.xlabel('Vds(V)')
plt.ylabel('Ids(uA)')

textstr = '\n'.join((
    'Square Law model',
    'Oxide thickness(nm):'+str(t_ox),
    'Channel length(nm):'+ str(L_c),
    'Channel width(nm):'+str(width),
    'Threshold voltage(V):'+str(V_t))
)
plt.yscale('log')
plt.text(0, 2500,textstr)

plt.show()