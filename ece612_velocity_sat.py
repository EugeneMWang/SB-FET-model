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
# V_t = 0.05 # unit: V
V_fb = -1.3 #unit: V
NA = 2.7*10**18 #unit: cm-3
NI = 10**10 #unit: cm-3
velocity_sat = 10**7 #unit: cm/s


psi_b = bolzcons*temperature/q*np.log(NA/NI)
C_ox = di_cons_sio2/(t_ox*10**-7)

#velocity saturation model
V_t = V_fb+2*psi_b+np.sqrt(2*q*di_cons_si*NA*(2*psi_b))/C_ox
m=1+np.sqrt(di_cons_si*q*NA/(4*psi_b))/C_ox
def charge_channel(v,vg):
    C_ox = di_cons_sio2/(t_ox*10**-7)
    return -C_ox*(vg-V_t-m*v)

def current_unsat(vds,vg):
    Fv=1/(1+mobility_Si*vds/(velocity_sat*(L_c*10**(-7))))
    integral= quad(charge_channel,0,vds,args=(vg,))[0]
    return -Fv*(width/L_c)*mobility_Si*integral

vds = np.linspace(0,1,100)
def Voltage_sat_find(vg):
    list=[]
    for i in vds:
        list.append(current_unsat(i,vg))
    return list.index(max(list)),max(list)
print(Voltage_sat_find(1.0))

def current(vd,vg):
    if vd<=vds[Voltage_sat_find(vg)[0]]:
        return current_unsat(vd,vg)
    else:
        return Voltage_sat_find(vg)[1]



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
    'Velocity saturation model(long channel)',
    'Oxide thickness(nm):'+str(t_ox),
    'Channel length(nm):'+ str(L_c),
    'Channel width(nm):'+str(width),
    'Flat band voltage(V):'+str(V_fb),
    'Gate range(V):'+str(min(Vg))+'to'+str(max(Vg)),
    'Doping Na:'+str(NA)
))

plt.text(0.2, 0.002,textstr)

plt.show()