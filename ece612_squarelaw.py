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
#square law model
def charge_channel(vy,vg):
    C_ox = di_cons_sio2/(t_ox*10**-7)
    return -C_ox*(vg-V_t-vy)

def current(vds,vg):
    if vg-V_t <= vds:
        integral = quad(charge_channel, 0, vg-V_t, args=(vg,))[0]
    else:
        integral= quad(charge_channel,0,vds,args=(vg,))[0]
    return -width/L_c*mobility_Si*integral


Vg = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0] #unit: V
vds = np.linspace(0,1)

y=[]

for v in Vg:
    y.append([])
    for i in vds:
        y[-1].append(current(i,v)*10**(6))

for i in y:
    plt.plot(vds,i)
plt.xlabel('Vds(V)')
plt.ylabel('Ids(uA)')

textstr = '\n'.join((
    'Square Law model',
    'Oxide thickness(nm):'+str(t_ox),
    'Channel length(nm):'+ str(L_c),
    'Channel width(nm):'+str(width),
    'Threshold voltage(V):'+str(V_t),
    'Gate range(V):'+str(min(Vg))+'to'+str(max(Vg))
))

plt.text(0, 2500,textstr)

plt.show()