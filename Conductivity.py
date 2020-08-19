import matplotlib.pyplot as plt
import numpy as np
from fdint import *

#############constants##########################
kB = 1.38064852 * 10**-23               #Boltzman Constant SI unit
e = 1.6*10**-19                         #Charge on electron SI unit
h = 6.626*10**-34                       #Planks Constant SI unit
pi = np.pi                              #pi
m_e = 9.10938356 * 10**-31              #mass of an electron in SI unit

mp = 4.7                               #Seebeck mass of valence band
mn =10.5                                #Seebeck mass of Conduction band

wmp = 0.03
wmn = 0.07

T1 = 10
T2 = 800

Ee = 0.7

Efs = -1                             #Starting Fermi-level
Efe = 1                              #Ending Fermi-level
npts = 1500                              #Number of points between starting and ending Ef

tpts = 100 

data=np.loadtxt('resistivity.txt')

T = [T for T in np.linspace(T1,T2,tpts)]                                  #Temperature for hall-measurement in K
T = np.array(T)

Sig = []
eF = []

Efs = -0.3                             #Starting Fermi-level
Efe = 0.1                              #Ending Fermi-level

Ef = [Ef for Ef in np.linspace(Efs,Efe,npts)]
Ef = np.array(Ef)

#########################################p-type BAND####################################################
Eg = 0.5                                #Band-gap in eV

#wmp = 0.03             
#wmn = 0.07

for i in range(tpts):

	wmp = 0.03*(300/T[i])**1.5             
	wmn = 0.07*(300/T[i])**1.5 

	p = [(4*pi*(((2*mp*m_e*kB*T[i])/(h**2))**1.5)*fdk(0.5,(-e*Ef)/(kB*T[i])))/(10**6) for Ef in np.linspace(Efs,Efe,npts)]          #hole carrier concentration in cm-3
	n = [(4*pi*(((2*mn*m_e*kB*T[i])/(h**2))**1.5)*fdk(0.5,(e*(Ef-Eg))/(kB*T[i])))/(10**6) for Ef in np.linspace(Efs,Efe,npts)]       #hole carrier concentration in cm-3

	p=np.array(p)
	n=np.array(n)

	c = p-n-4*10**20
	c = np.absolute(c)
	j=np.argmin(c)

	sigP = (8/3)*np.pi*e*((2*m_e*kB*T[i])/(h**2))**1.5*fdk(0,(-e*Ef[j])/(kB*T[i]))*(wmp) 
	sigN = (8/3)*np.pi*e*((2*m_e*kB*T[i])/(h**2))**1.5*fdk(0,(e*(Ef[j]-Eg))/(kB*T[i]))*(wmn) 

	sig = sigN + sigP
#	print(sig) 
	Sig.append(sig)
	eF.append(Ef[j])

########################################Two band model###################################################
Sig = np.array(Sig)
res2 = 1/Sig
res2=res2*10**8

eF = np.array(eF)

plt.rcParams['figure.figsize'] = 8, 6

#plt.plot(T,eF,s=5)

plt.plot(T,res2)

plt.scatter(data[:,0],data[:,1],s=10)

plt.xlabel('Temperature (K)', fontsize= 16)
plt.ylabel('Resistivity ($\mu \Omega cm$)', fontsize =16)

plt.xticks(fontsize= 16)
plt.yticks(fontsize= 16)

#plt.xscale('log')
plt.yscale('log')

#plt.xlim(0,300)
#plt.ylim(100,10000)

plt.show()
