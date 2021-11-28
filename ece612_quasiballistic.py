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
width = 40 # channel width; unit: nm
L_c = 100 # channel length; unit: nm
# V_t = 0.05 # unit: V
V_fb = -1.3 #unit: V
NA = 2.7*10**18 #unit: cm-3
NI = 10**10 #unit: cm-3
psi_b = bolzcons*temperature/q*np.log(NA/NI)
C_ox = di_cons_sio2/(t_ox*10**-7)

#quasi ballistic model
T=0.15
R=1-T
V_t = V_fb+2*psi_b+np.sqrt(2*q*di_cons_si*NA*(2*psi_b))/C_ox
velocity_therm=np.sqrt(2*bolzcons*temperature/(np.pi*m_si))*100

def bracketin(vds):
    return (1-R-(1-R)*np.exp(-q*vds/(bolzcons*temperature)))/(1+R+(1-R)*np.exp(-q*vds/(bolzcons*temperature)))

def current(vds,vg):
    return width*10**(-7)*C_ox*(vg-V_t)*velocity_therm*bracketin(vds)

vds = np.linspace(0,1,100)
Vg = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0] #unit: V



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
    'Quasi Ballistic model',
    'Oxide thickness(nm):'+str(t_ox),
    'Channel length(nm):'+ str(L_c),
    'Channel width(nm):'+str(width),
    'Flat band voltage(V):'+str(V_fb),
    'Gate range(V):'+str(min(Vg))+'to'+str(max(Vg)),
    'Doping Na:'+str(NA),
    'T:'+str(T)
))

plt.text(0.2, 0.002,textstr)

plt.show()