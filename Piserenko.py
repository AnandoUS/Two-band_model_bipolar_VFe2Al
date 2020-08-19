import matplotlib.pyplot as plt
import numpy as np
from fdint import *

#############constants##########################
kB = 1.38064852 * 10**-23               		#Boltzman Constant SI unit
e = 1.6*10**-19                         		#Charge on electron SI unit
h = 6.626*10**-34                       		#Planks Constant SI unit
pi = np.pi                              		#pi
m_e = 9.10938356 * 10**-31              		#mass of an electron in SI unit
mp = 4.7                               			#Seebeck mass of valence band
mn =12.75                                		#Seebeck mass of Conduction band
T = 300                                 		#Temperature for hall-measurement in K

Ee = 0.7

Efs = -0.7                             			#Starting Fermi-level
Efe = 1                              			#Ending Fermi-level
npts = 1500                              		#Number of points between starting and ending Ef

v = 5.76**3/(16*10**4)

#########################################p-type BAND####################################################
wmp = 0.03

p = [(-4*pi*(((2*mp*m_e*kB*T)/(h**2))**1.5)*fdk(0.5,(-e*Ef)/(kB*T)))/(10**6) for Ef in np.linspace(Efs,Efe,npts)]            #hole carrier concentration in cm-3
sigP = [(8/3)*np.pi*e*((2*m_e*kB*T)/(h**2))**1.5*fdk(0,(-e*Ef)/(kB*T))*(wmp) for Ef in np.linspace(Efs,Efe,npts)]
SP = [10**6*(kB/e)*((2*fdk(1,(-e*Ef)/(kB*T))/fdk(0,(-e*Ef)/(kB*T)))-((-e*Ef)/(kB*T))) for Ef in np.linspace(Efs,Efe,npts)]

p=np.array(p)
sigP = np.array(sigP)
SP = np.array(SP)

#########################################n-type BAND####################################################
Eg = 0.02                                #Band-gap in eV
wmn = 0.07

n = [(4*pi*(((2*mn*m_e*kB*T)/(h**2))**1.5)*fdk(0.5,(e*(Ef-Eg))/(kB*T)))/(10**6) for Ef in np.linspace(Efs,Efe,npts)]         		 #hole carrier concentration in cm-3
sigN = [(8/3)*np.pi*e*((2*m_e*kB*T)/(h**2))**1.5*fdk(0,(e*(Ef-Eg))/(kB*T))*(wmn) for Ef in np.linspace(Efs,Efe,npts)]
SN = [-10**6*(kB/e)*((2*fdk(1,(e*(Ef-Eg))/(kB*T))/fdk(0,(e*(Ef-Eg))/(kB*T)))-((e*(Ef-Eg))/(kB*T))) for Ef in np.linspace(Efs,Efe,npts)]

n=np.array(n)
sigN = np.array(sigN)
SN = np.array(SN)

########################################Two band model###################################################
S = (sigN*SN+sigP*SP)/(sigN+sigP)
N = (p+n)

plt.rcParams['figure.figsize'] = 9.1, 6.5

#################DATA##########################
Co=np.loadtxt('Co_2009.txt', skiprows=1)
Ge=np.loadtxt('Ge_2006.txt', skiprows=1)
Si=np.loadtxt('Si_2001.txt', skiprows=1)
Si1=np.loadtxt('Si_2008.txt', skiprows=1)
Si2=np.loadtxt('Si_2009.txt', skiprows=1)
Sn=np.loadtxt('Sn_2009.txt', skiprows=1)
Mo=np.loadtxt('Mo_2002.txt', skiprows=1)
W=np.loadtxt('W.txt', skiprows=1)
Pt=np.loadtxt('Pt.txt', skiprows=1)
Ti=np.loadtxt('Ti.txt', skiprows=1)
Zr=np.loadtxt('Zr.txt', skiprows=1)


plt.scatter((Co[:,0]-6)/v, Co[:,1],s=50)				#Co		BLUE
plt.scatter((Ge[:,0]-6)/v, Ge[:,1],s=50)				#Ge		ORANGE
plt.scatter((Si[:,0]-6)/v, Si[:,1],s=50)				#Si_2001	GREEN
plt.scatter((Si1[:,0]-6)/v, Si1[:,1],s=50)				#Si_2008	RED
plt.scatter((Si2[:,0]-6)/v, Si2[:,1],s=50)				#Si_2009	PURPLE
plt.scatter((Sn[:,0]-6)/v, Sn[:,1],s=50)				#Sn_2009	BROWN
plt.scatter((Mo[:,0]-6)/v, Mo[:,1],s=50)				#Mo_2002	PINK
plt.scatter((W[:,0]-6)/v, W[:,1],s=50)					#W		GREY
plt.scatter((Pt[:,0]-6)/v, Pt[:,1],s=50)				#Pt		MUDDY_YELLOW
plt.scatter((Ti[:,0]-6)/v, Ti[:,1],c='cadetblue',marker='p',s=67.5)	#Ti		MUDDY_YELLOW
plt.scatter((Zr[:,0]-6)/v, Zr[:,1],c='lightcoral',marker='p',s=67.5)	#Zr

plt.plot(N/10**20,S,linewidth=2)

#plt.plot([N[0],N[npts-1]], [25,25])

plt.plot([N[0],N[npts-1]], [0,0],linewidth=0.3,c='k')
plt.plot([0,0], [-200,200],linewidth=0.3,c='k')

plt.xlabel('$n$ (X 10$^{20}$ cm$^{-3}$)', fontsize= 20)
plt.ylabel('Thermopower |S| ($\mu$V/K)', fontsize =20)

plt.xticks(fontsize= 18)
plt.yticks(fontsize= 18)

plt.xlim(-100,250)
plt.ylim(-200,200)

plt.show()
