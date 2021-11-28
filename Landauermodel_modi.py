#Author: MINGYI WANG @ Purdue University
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad
# import sympy as sy
# x,E_variable,Ec_ch = sy.symbols("x, E_variable,Ec_ch")



q=1.6*10**(-19)
m0=9.11*10**(-31) #unit:  kg
h = 6.626*10**(-34) #unit: J s
h_bar = 1.055*10**(-34) #unit: J s
mh = 0.4*m0 #For Te, effective mass for holes is 0.76*m0, effective mass for electrons is 0.3*m0.
me = 0.6*m0
xm = 0 # Assuming the contact is 0 point
lambda_x = 38*10**(-9) #unit: m
Ev0= -0.3*q #unit: J (we assume SBH is 0.xeV)
Ec0 = 0.3*q
temperature = 300# unit:K
gc = 2 #degeneracy
gv = 2
bolzcons = 1.38*10**(-23) #unit: JK-1
vds = 0.1 #unit: V

##########################################################electrons###########################################################

def integrand_T (x,E_variable, Ec_ch): # this is basically the kai
    E_barrier = ((Ec0-Ec_ch)/lambda_x)*(lambda_x-x)-(E_variable-Ec_ch)
    return (1/h_bar)*np.sqrt(2*me*(E_barrier))

def INTE_function_T(E_variable,Ec_ch):
    if E_variable < Ec_ch:
        return 0
    elif E_variable >= Ec0:
        return 1
    else:
        end_tunneling = lambda_x * (Ec0-E_variable)/(Ec0-Ec_ch)
        # return np.exp(-sy.integrate(integrand_T(x,E_variable,Ec_ch),(x,xm,end_tunneling)))
        return np.exp(-quad(integrand_T,xm,end_tunneling,args = (E_variable,Ec_ch))[0])
        # return np.exp(-fixed_quad((integrand_T,xm,end_tunneling,n=5)[0]))

def Tunneling_all (E_variable,Ec_ch):
    Tunneling_source = INTE_function_T(E_variable,Ec_ch)
    Tunneling_drain = INTE_function_T(E_variable,Ec_ch)
    Reflac_source = 1-Tunneling_source
    Reflac_drain = 1-Tunneling_drain
    Tunneling = Tunneling_drain*Tunneling_source/(1-Reflac_source*Reflac_drain)
    return Tunneling

def number_mode(E_variable,Ec_ch):
    if E_variable > Ec_ch:
        return (gc/(h_bar*np.pi))*np.sqrt(2*me*(E_variable-Ec_ch))
    else:
        return 0

def fermi_distri(E_variable):
    return 1/((np.exp(E_variable/(bolzcons*temperature)))+1)

def fermi_distri_1(E_variable):
    return np.exp(-E_variable/(bolzcons*temperature))/(np.exp(-E_variable/bolzcons*temperature))+1

def fermi_reduced(E_variable):
    x = fermi_distri(E_variable)-fermi_distri(E_variable-vds * q)
    return abs(x)

def integrand_E(E_variable,Ec_ch):
    return INTE_function_T(E_variable,Ec_ch)*number_mode(E_variable,Ec_ch)\
        *fermi_reduced(E_variable)

def INTE_function (Ec_ch):
    return quad(integrand_E,Ec_ch,5*q,args=Ec_ch,epsabs=1.49e-25)[0]
    # return sy.integrate(integrand_E(E_variable,Ec_ch),(E_variable,Ec_ch,2*q))
def current_electron(Ec_ch):
    return INTE_function(Ec_ch)*2*q/h

############################################################################################################################

############################################################Holes###########################################################
def integrand_T_v (x,E_variable, Ev_ch): # this is basically the kai
    E_barrier = ((Ev_ch-Ev0)/lambda_x)*(lambda_x-x)-(Ev_ch-E_variable)
    return (1/h_bar)*np.sqrt(2*mh*(E_barrier))

def INTE_function_T_v(E_variable,Ev_ch):
    if E_variable > Ev_ch:
        return 0
    elif E_variable <= Ev0:
        return 1
    else:
        end_tunneling = lambda_x * (E_variable-Ev0)/(Ev_ch-Ev0)
        return np.exp(-quad (integrand_T_v,xm,end_tunneling, args = (E_variable,Ev_ch))[0])

def Tunneling_all_v (E_variable,Ev_ch):
    Tunneling_source = INTE_function_T_v(E_variable,Ev_ch)
    Tunneling_drain = INTE_function_T_v(E_variable,Ev_ch)
    Reflac_source = 1-Tunneling_source
    Reflac_drain = 1-Tunneling_drain
    Tunneling = Tunneling_drain*Tunneling_source/(1-Reflac_source*Reflac_drain)
    return Tunneling

def number_mode_v(E_variable,Ev_ch):
    if E_variable < Ev_ch:
        return (gv/(h_bar*np.pi))*np.sqrt(2*mh*(Ev_ch-E_variable))
    else:
        return 0

def fermi_distri_1_v(E_variable):
    return np.exp(-E_variable/(bolzcons*temperature))/(np.exp(-E_variable/bolzcons*temperature))+1

def fermi_distri_v(E_variable):
    return 1-1/((np.exp(E_variable/(bolzcons*temperature)))+1)

def fermi_reduced_v(E_variable):
    x = fermi_distri_v(E_variable)-fermi_distri_v(E_variable-vds * q)
    return abs(x)

def integrand_E_v(E_variable,Ev_ch):
    return INTE_function_T_v(E_variable,Ev_ch)*number_mode_v(E_variable,Ev_ch)\
        *fermi_reduced_v(E_variable)

def INTE_function_v (Ev_ch):
    return quad(integrand_E_v,-5*q,Ev_ch,args=Ev_ch,epsabs=1.49e-18)[0]

def current_hole(Ev_ch):
    return INTE_function_v(Ev_ch)*2*q/h



#test code
channelshift = np.linspace(-0.4*q,0.4*q,endpoint=False)
Ec_range=[]
Ev_range=[]
for i in channelshift:
    Ec_range.append(Ec0+i)
    Ev_range.append(Ev0+i+vds*q)

y1=[]
y2=[]
yz=[]
for i in Ec_range:
    y1.append((current_electron(i)/10**6))

for i in Ev_range:
    y2.append(current_hole(i)/10**6)

l=0
for i in y1:
   yz.append(y2[l]+i)
   l+=1

x=[]
for i in channelshift:
    x.append(-i/q*10)
plt.yscale('log')

plt.plot(x,yz,'o',mfc='none',label='$I_{hole}$+$I_{electron}$')
print(x)
print(yz)
plt.plot(x,y1,label='$I_{electron}$')
plt.plot(x,y2,label='$I_{hole}$')
plt.ylim(5*10**-9,10**-5)
plt.xlabel('$V_{g}$(V)')
plt.ylabel('$I_{ds}(A)$')
plt.legend()
plt.show()