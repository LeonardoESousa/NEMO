import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import os
from scipy.linalg import expm
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import matplotlib.pylab as pylab
#'figure.figsize': (8, 4.8),
#params = {'xtick.top': False,
#			'xtick.direction': 'out',
#			'ytick.direction': 'out',
#			'legend.fontsize': 13,
#			'figure.titlesize': 15,
#			'axes.titlesize': 22,
#			'font.family': 'Times New Roman',
#         'axes.labelsize': 18,
#         'xtick.labelsize':15,
#         'ytick.labelsize':15,
#		 "errorbar.capsize":2.5}
params = {'xtick.top': False,
			'xtick.direction': 'out',
			'ytick.direction': 'out',
			'legend.fontsize': 18,
			'figure.titlesize': 18,
			'axes.titlesize': 25,
			'font.family': 'Times New Roman',
         'axes.labelsize': 21,
         'xtick.labelsize':18,
         'ytick.labelsize':18,
		 "errorbar.capsize":2.5}
pylab.rcParams.update(params)

tipos = ['S1','T1', 'T2']

def painel():
    socsS1 = np.loadtxt('S1\SOCs.txt')
    socsT1 = np.loadtxt('T1\SOCs.txt')
    socsS1 = socsS1[:,0]*0.12398#/1000 results in meV
    socsT1 = socsT1[:,0]*0.12398#/1000 results in meV
    
    opts = np.loadtxt('S1\opt\SOCs.txt').flatten()*0.12398
    optt = np.loadtxt('T1\opt\SOCs.txt').flatten()*0.12398
    
    opts = opts[0]
    optt = optt[0]
    
    
    exs = np.loadtxt('S1\ENGs.txt')
    ext = np.loadtxt('T1\ENGs.txt')
    n = int(np.shape(exs)[1]/2)
    
    deltaS = exs[:,0] - exs[:,n]
    deltaT = ext[:,0] - ext[:,n]
    
    ops = np.loadtxt('S1\opt\ENGs.txt')
    opt = np.loadtxt('T1\opt\ENGs.txt')
    ops = ops[0] - ops[n]
    opt = opt[0] - opt[n]
    
    
    fig, ax = plt.subplots(1,4, figsize=(20,5))
    A = ax[0].hist(socsS1, bins=int(np.sqrt(len(socsS1))) ,density = False, label="S1 Geoms")
    ax[0].plot(np.zeros(100)+opts,np.linspace(0,max(A[0]),100),':',lw=4,color='red')
    B = ax[1].hist(socsT1, bins=int(np.sqrt(len(socsS1))) ,density = False, label="T1 Geoms")
    ax[1].plot(np.zeros(100)+optt,np.linspace(0,max(B[0]),100),':',lw=4,color='red')
    
    A2 = ax[2].hist(deltaS, bins=int(np.sqrt(len(deltaS))) ,density = False, label="S1 Geoms")
    ax[2].plot(np.zeros(100)+ops,np.linspace(0,max(A2[0]),100),':',lw=4,color='red')
    B2 = ax[3].hist(deltaT, bins=int(np.sqrt(len(deltaT))) ,density = False, label="T1 Geoms")
    ax[3].plot(np.zeros(100)+opt,np.linspace(0,max(B2[0]),100),':',lw=4,color='red')
    
    
    ax[0].legend(loc="best")
    ax[1].legend(loc="best")
    ax[2].legend(loc="best")
    ax[3].legend(loc="best")
    ax[0].set_xlabel("SOC (meV)")
    ax[1].set_xlabel("SOC (meV)")
    ax[0].set_ylim([0,1.1*max(A[0])])
    ax[1].set_ylim([0,1.1*max(B[0])])
    ax[0].set_xlim([0,1.2*max(A[1])])
    ax[1].set_xlim([0,1.2*max(B[1])])
    
    
    ax[0].set_yticks(np.arange(0,1.1*max(A[0] ),np.ceil(1.1*max(A[0] )/5)))
    ax[1].set_yticks(np.arange(0,1.1*max(B[0] ),np.ceil(1.1*max(B[0] )/5)))
    ax[2].set_yticks(np.arange(0,1.1*max(A2[0]),np.ceil(1.1*max(A2[0])/5)))
    ax[3].set_yticks(np.arange(0,1.1*max(B2[0]),np.ceil(1.1*max(B2[0])/5)))
    
    ax[0].set_xticks(np.round(np.linspace(0,1.1*max(A[1] ),4),1))#(1.1*max(A[1] )/4),1)))
    ax[1].set_xticks(np.round(np.linspace(0,1.1*max(B[1] ),4),1))#(1.1*max(B[1] )/4),1)))
    ax[2].set_xticks(np.round(np.linspace(0,1.2*max(A2[1]),4),1))#(1.1*max(A2[1])/4),1)))
    ax[3].set_xticks(np.round(np.linspace(0,1.2*max(B2[1]),4),1))#(1.1*max(B2[1])/4),1)))
    
    print(np.round(np.linspace(0,1.2*max(A2[1]),4),1),np.round(np.linspace(0,1.2*max(B2[1]),4),1))
    
    ax[2].set_xlabel("$\Delta E_{ST}$ (eV)")
    ax[3].set_xlabel("$\Delta E_{ST}$ (eV)")
    ax[2].set_ylim([0,1.1*np.round(max(A2[0]),1)])
    ax[3].set_ylim([0,1.1*np.round(max(B2[0]),1)])
    ax[2].set_xlim([0,1.2*np.round(max(A2[1]),1)])
    ax[3].set_xlim([0,1.2*np.round(max(B2[1]),1)])
    
    ax[0].set_title('e')
    ax[1].set_title('f')
    ax[2].set_title('g')
    ax[3].set_title('h')
    
    fig.tight_layout()
    #plt.show()
    plt.savefig("painel.pdf")

 
def gauss(x,v,s):
    y =  (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-v)/s)**2)
    return y

def isc(tipo,n_triplet):
    T = 300 #K
    hbar = 6.582*(10**(-16)) #Reduced Planck's constant
    kb = 8.617*(10**(-5)) # Boltzmann constant    
    socs_complete = np.loadtxt(tipo+'\SOCs.txt')*0.12398/1000
    exs = np.loadtxt(tipo+'\ENGs.txt')
    n = np.shape(exs)[1]
    n = int(n/2)
    total = 0
    rates = []
    for j in range(np.shape(socs_complete)[1]):
        try:
            lambdas = np.loadtxt(tipo+'\lambda.txt')[j]    #(tipo+'\opt\lambdas.txt')[0]
        except:
            try:
                lambdas = np.loadtxt(tipo+'\lambda.txt')[0]
            except:
                lambdas = np.loadtxt(tipo+'\lambda.txt')    
        if tipo == 'S1':
            delta = exs[:,n+j] - exs[:,0]  #Tn (final) - S1 (initial)    
        elif tipo == 'T1' or tipo == 'T2':
            delta = exs[:,0] - exs[:,n+j] #S1 (final) - Tn (initial)

        delta = delta.flatten()
        const = (2*np.pi/hbar)
        rate = 0
        socs = socs_complete[:,j]
        sigma = np.sqrt(2*lambdas*kb*T + (kb*T)**2)
        for i in range(np.shape(socs)[0]):
            rate += (socs[i]**2)*gauss(delta[i]+lambdas,0,sigma)
        rate /= len(socs)
        rate *= const
        rates.append(rate)
        total += rate
        print(np.round(np.mean(delta),2),'ISC rate:',tipo, j+1,format(rate, "5.2e"),'s^-1')
        #print('ISC lifetime:', tipo,format(rate**-1, "5.2e"),'s')
    print('Total',format(total, "5.2e"))
    print('---------------------------------------------------------')
    return rates[n_triplet-1]  #[0] 

def avgmlj(tipo):
    Lambdas = np.loadtxt(tipo+'\opt\lambdas.txt')
    T = 300 #K
    hbar = 6.582*(10**(-16)) #Reduced Planck's constant
    kb = 8.617*(10**(-5)) # Boltzmann constant    
    socs = np.loadtxt(tipo+'\SOCs.txt')*0.12398/1000
    exs = np.loadtxt(tipo+'\ENGs.txt')
    n = np.shape(exs)[1]
    n = int(n/2)
    s1 = np.tile(exs[:,0],(n,1)).transpose()   
    t1 = np.tile(exs[:,n],(n,1)).transpose()
    if tipo == 'S1':
        delta = exs[:,n] - exs[:,0]#exs[:,0] - exs[:,n]#exs[:,n] - exs[:,0]
        Lambda = Lambdas[0]
    elif tipo == 'T1':
        delta = exs[:,0] - exs[:,n]
        Lambda = Lambdas[1]
    
    socs = socs[:,0]
    
    const = (socs**2)*(2*np.pi/hbar)*(1/(np.sqrt(4*np.pi*kb*T*Lambda)))
    lqm = 0.05 #eV
    wqm = 0.14878 # eV
    sqm = lqm/wqm
    rate = 0
    for i in range(np.shape(socs)[0]):
        rates = 0
        for n in range(0,50):
            exponent = (-1/(4*kb*T*Lambda))*(Lambda+ n*wqm + delta[i])**2
            rates += const[i]*(np.exp(-1*sqm)*(sqm**n)/np.math.factorial(n))*np.exp(exponent)
        rate += rates
    rate /= len(socs)
    #print(np.mean(SH),max(SH),min(SH))
    print('AVG rate:', tipo,format(rate, "5.2e"),'s^-1')
    #print('AVG lifetime:', tipo,format(rate**-1, "5.2e"),'s')
    #plt.hist(SH,bins=int(np.sqrt(len(SH))))
    #plt.show()
    return rate
 
def mlj(tipo):
    T = 300 #K
    Lambdas = np.loadtxt(tipo+'\opt\lambdas.txt')
    hbar = 6.582*(10**(-16)) #Reduced Planck's constant
    kb = 8.617*(10**(-5)) # Boltzmann constant    
    socs = np.loadtxt(tipo+'\opt\SOCs.txt')
    socs = socs.flatten()[0]*0.12398/1000
    exs = np.loadtxt(tipo+'\opt\ENGs.txt')
    n = len(exs)
    n = int(n/2)
    if tipo == 'S1':
        delta = exs[n] - exs[0]#exs[0] - exs[n]#exs[n] - exs[0]
        Lambda = Lambdas[0]
    elif tipo == 'T1':
        delta = exs[0] - exs[n]
        Lambda = Lambdas[1]
        
    #print('SOC',socs)
    const = (socs**2)*(2*np.pi/hbar)*(1/(np.sqrt(4*np.pi*kb*T*Lambda)))
    lqm = 0.05 #eV
    wqm = 0.14878 # eV
    sqm = lqm/wqm
    rate = 0
    for n in range(0,50):
        exponent = (-1/(4*kb*T*Lambda))*(Lambda+ n*wqm + delta)**2
        rate += (np.exp(-1*sqm)*(sqm**n)/np.math.factorial(n))*np.exp(exponent)
    rate *= const
    #print(tipo,socs,delta,Lambda)
    print('MLJ rate:',tipo,format(rate, "5.2e"),'s^-1')
    #return rate 

def abs_spectra():
    T = 300 #K
    nr = 1.4969
    epsilon0 = 8.854187817*10**(-12) #F/m
    hbar = 6.582119514*10**(-16) #eV s
    hbar2 = 1.054571800*10**(-34) #J s
    mass = 9.10938356*10**(-31) # kg
    c = 299792458 #m/s
    e = 1.60217662*10**(-19) #C
    pi = np.pi
    kb = 8.6173303*(10**(-5)) #eV/K
    amu = 1.660539040*10**(-27) #kg
    
    const = (np.pi*(e**2)*hbar)/(2*nr*mass*c*epsilon0)*10**(20)
    exs = np.loadtxt('S0/ENGs.txt')
    exs = exs[:,:int(np.shape(exs)[1]/2)]
    F = np.loadtxt('S0/OSCs.txt')
    
    x = np.linspace(min(exs.flatten())-3*(kb*T), max(exs.flatten())+ 3*(kb*T), 1000)
    y = np.zeros(len(x))
    for n in range(np.shape(F)[1]):
        deltaE = exs[:,n]
        f = F[:,n]
        
        
        for i in range(0,len(deltaE)):
            y += const*f[i]*gauss(x,deltaE[i],kb*T) 
    y = y/len(deltaE)
    print('Max Abs',x[list(y).index(max(y))])
    plt.plot(x,y)
    plt.show()
    #x = (10**9)*4.135667696*10**(-15)*c/x
    with open('abs.lx', 'w') as f:
        f.write('#Energy(eV)    Cross-Section(A^2)\n')
        for i in range(len(x)):
            f.write('    '.join([str(x[i]),str(y[i])])+'\n')
    
    #return rate, x, y

def emi_spectra():
    T = 300 #K
    nr = 1.4969
    pot = 2
    epsilon0 = 8.854187817*10**(-12) #F/m
    hbar = 6.582119514*10**(-16) #eV s
    hbar2 = 1.054571800*10**(-34) #J s
    mass = 9.10938356*10**(-31) # kg
    c = 299792458 #m/s
    e = 1.60217662*10**(-19) #C
    pi = np.pi
    kb = 8.6173303*(10**(-5)) #eV/K
    amu = 1.660539040*10**(-27) #kg
    
    const = ((e**2)/(2*np.pi*hbar*mass*(c**3)*epsilon0))
    exs = np.loadtxt('S1/ENGs.txt')
    deltaE = exs[:,0]                            ##- 0.4
    f = np.loadtxt('S1/OSCs.txt')[:,0]
    
    x = np.linspace(min(deltaE)-3*(kb*T), max(deltaE)+ 3*(kb*T), 1000)
    y = np.zeros(len(x))
    for i in range(0,len(deltaE)):
        y += (nr**pot)*const*(deltaE[i]**2)*f[i]*gauss(x,deltaE[i],kb*T) 
    y = y/len(deltaE)
    print('Max Fluor',x[list(y).index(max(y))])
    IntEmi = np.trapz(y,x)
    rate = (IntEmi/hbar)
    print('Emission rate:',format(rate, "5.2e"),'s^-1')
    #print('Emission Lifetime:',format(rate**-1, "5.2e"),'s')
    #x = (10**9)*4.135667696*10**(-15)*c/x
    with open('fluor.lx', 'w') as f:
        f.write('#Energy(eV)    DiffRate\n')
        for i in range(len(x)):
            f.write('    '.join([str(x[i]),str(y[i])])+'\n')
    return rate, x, y

def phosph(n_triplet):
    T = 300 #K
    nr = 1.4969
    pot = 2
    epsilon0 = 8.854187817*10**(-12) #F/m
    hbar = 6.582119514*10**(-16) #eV s
    hbar2 = 1.054571800*10**(-34) #J s
    mass = 9.10938356*10**(-31) # kg
    c = 299792458 #m/s
    e = 1.60217662*10**(-19) #C
    pi = np.pi
    kb = 8.6173303*(10**(-5)) #eV/K
    amu = 1.660539040*10**(-27) #kg
    
    const = (1/3)*(1/(3*np.pi*(hbar**3)*(c**3)*epsilon0))
    exs = np.loadtxt('T'+str(n_triplet)+'/ENGs.txt')
    n = np.shape(exs)[1]
    n = int(n/2)
    deltaE = exs[:,n+n_triplet-1]
    f = np.loadtxt('T'+str(n_triplet)+'/RTP.txt')[:,n_triplet-1]
    x = np.linspace(max(min(deltaE)-3*(kb*T),0), max(deltaE)+ 3*(kb*T), 1000)
    y = np.zeros(len(x))
    count = 0
    for i in range(0,len(deltaE)):
        y += (nr**pot)*const*(abs(deltaE[i])**3)*f[i]*gauss(x,deltaE[i],kb*T) 
        count += 1
    y = y/count
    y *= 6.2415093433E+18
    print('Max Phosph',x[list(y).index(max(y))])
    
    IntEmi = np.trapz(y,x)
    rate = (IntEmi/hbar)
    print('RTP rate:',format(rate, "5.2e"),'s^-1')
    #print('RTP Lifetime:',format(rate**-1, "5.2e"),'s')
    #plt.hist(shift,bins=int(np.sqrt(len(deltaE+shift))))
    #x= (10**9)*4.135667696*10**(-15)*c/x
    plt.plot(x,y)
    plt.show()
    with open('phosph_'+str(n_triplet)+'_.lx', 'w') as f:
        f.write('#Energy(eV)    DiffRate\n')
        for i in range(len(x)):
            f.write('    '.join([str(x[i]),str(y[i])])+'\n')
    return rate, x, y


def kinetics(n_triplet):
    kisc =  isc('S1',n_triplet)   
    krisc = isc('T'+str(n_triplet),n_triplet)  
    print('KISC:',format(kisc, "5.2e"))
    print('KRISC:',format(krisc, "5.2e"))
    kem, x, y = emi_spectra() 
    kp, x, y = phosph(n_triplet)
    #print('KF:',format(kem, "5.2e"))
    #print('KP:',format(kp, "5.2e"))
    
    
    knr =  0   #np.loadtxt('knr.txt')
    M = np.zeros((5,5))
    M[0,0] = -(kem + kisc +knr)
    M[0,1] = krisc
    M[1,0] = kisc
    M[1,1] = -(krisc+kp)
    M[2,0] = kem
    M[3,1] = kp
    M[4,0] = knr


    P = np.zeros((5,1))
    P[0,0] = 100
    R = P
    deltat = 10**(-9)
    x = [0] #deltat*np.linspace(0,10000,10000)
    #for i in range(1,20000):   #,len(x)):
    while np.sum(P[2:,0]) < 99.999:#9999:
        P = np.matmul(expm(M*deltat),P)
        R = np.hstack((R,P))
        x.append(x[-1]+deltat)
        print(P[0,0],P[1,0],np.sum(P[2:,0]),end="\r", flush=True)
        if x[-1] > 10**-7:
            deltat = 10**-8
    x = np.array(x)
    print('\n\nsoma',np.sum(R[:,-1]))
    np.savetxt('kinetics.txt', R, delimiter='\t')
    np.savetxt('x.txt', x, delimiter='\t')
    #return x, R

def func(x,a1,b1,a2,b2):
    return a1*np.exp(-x*b1) + a2*np.exp(-x*b2)

def plot(n_triplet):
    x = np.loadtxt('x.txt')
    R = np.loadtxt('kinetics.txt')
    fig, ax = plt.subplots(2,1,figsize=(5,7))
    y1 = R[0,:]
    y2 = R[1,:]
    ax[0].plot(x,y1,label='S1',lw=2)
    ax[0].plot(x,y2,label='T'+str(n_triplet),lw=2)
 
    y3 = R[2,-1]-R[2,:]
    y4 = R[3,-1]-R[3,:]
    
    ax[1].plot(x,y3,label=r'$S_1 \to S_0$',lw=2,color='blue')
    ax[1].plot(x,y4,label=r'$T_'+str(n_triplet)+r' \to S_0$',lw=2,color='red')
    
    
    popt, pcov = curve_fit(func, x, y3,p0=[10,1E9,10,1E6])
    print(format(1/popt[1], "5.2e"),format(1/popt[3], "5.2e"))
    print(popt)
    
    #6.87894669e+01 ,6.98972848e+05 ,2.30961228e+01 ,7.15766116e+07 CZIPN
    #8.52047687e+01 ,7.02254646e+07 ,1.47884150e+01 ,1.19988239e+06 ORTH
    
    
    #ax[1].plot(x,func(x,popt[0],popt[1],popt[2],popt[3]),label=r'Exp',lw=2,color='black')
    
    print(R[2,-1])
    print(R[3,-1])
    ax[0].legend()
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    #ax[0].set_xlim([10**-9,10**1])
    ax[0].set_ylim([10**-2,2.0*10**2])
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('Population (%)')
    #ax[0].set_xticks([10**-9,10**-7,10**-5,10**-3,10**-1])
    #ax[1].set_yscale('log')
    ax[1].legend()
    for j in range(len(y3)):
        if y3[j]/y3[0] < 0.1:
            limite = j
            break
    #ax[1].set_xlim([-0.05*x[limite]*1E6,x[limite]*1E6])
    #ax[1].set_ylim([0,105])
    ax[1].set_xlabel('Time (s)')
    ax[1].set_ylabel('Intensity (a.u.)')
    
    ax[0].set_title('a')
    ax[1].set_title('b')
    plt.tight_layout()
    plt.savefig('kinetics.png')
    #plt.show()
    
def spectra():
    absor  = np.loadtxt('abs.lx')
    fluor  = np.loadtxt('fluor.lx')
    phosph = np.loadtxt('phosph_1_.lx')
    phosph2 = np.loadtxt('phosph_2_.lx')
    print(absor)
    
    fig, ax = plt.subplots(1,4,figsize=(12,6))
    #y *= 1E9
    #y2 *= 1E12
    ax[0].plot(absor[:,0],absor[:,1],   color='blue', lw=2, label='Abs')
    ax[1].plot(fluor[:,0],fluor[:,1],   color='red',  lw=2, label='Fluor')
    ax[2].plot(phosph[:,0],phosph[:,1], color='green',  lw=2, label='Phosph $T_1$')
    ax[3].plot(phosph2[:,0],phosph2[:,1], color='green',  lw=2, label='Phosph $T_2$')


    #ax[0].set_xlabel('Wavelength (nm)')
    #ax[1].set_xlabel('Wavelength (nm)')
    ax[0].set_xlabel('Energy (eV)')
    ax[1].set_xlabel('Energy (eV)')
    ax[2].set_xlabel('Energy (eV)')
    ax[3].set_xlabel('Energy (eV)')
    ax[0].set_ylabel(r'Cross Section ($\AA^2$)')
    ax[1].set_ylabel(r'Diff. rate')
    ax[2].set_ylabel(r'Diff. rate')
    ax[3].set_ylabel(r'Diff. rate')

    ax[0].legend(loc='best', fontsize=12)
    ax[1].legend(loc='best', fontsize=12)
    ax[2].legend(loc='best', fontsize=12)
    ax[3].legend(loc='best', fontsize=12)
    
    #ax[0].set_xticks(np.round(np.linspace(2.4,5,5),1))
    #ax[1].set_xticks(np.round(np.linspace(2.0,3.4,5),1))
    #ax[2].set_xticks(np.round(np.linspace(1.3,3.0,5),1))
    #ax[3].set_xticks(np.round(np.linspace(1.3,3.0,5),1))
    #ax[0].set_xlim([350,750])
    #ax[1].set_xlim([350,750])
    #ax[0].set_yticks([0,5,10,15,20,25])
    ax[0].set_ylim([0,1.2*max(absor[:,1])])
    ax[1].set_ylim([0,1.2*max(fluor[:,1])])
    ax[2].set_ylim([0,1.2*max(phosph[:,1])])
    ax[3].set_ylim([0,1.2*max(phosph2[:,1])])
    #ax[1].set_yticks([0,1,2,3,4])
    
    ax[0].set_title('a')
    ax[1].set_title('b')
    ax[2].set_title('c')
    ax[3].set_title('d')
    plt.tight_layout()
    plt.savefig('spectra.png')
    #plt.show()
   
def overlap():
    absor  = np.loadtxt('abs.lx')
    fluor  = np.loadtxt('fluor.lx')
    phosph = np.loadtxt('phosph_1_.lx')
    phosph2 = np.loadtxt('phosph_2_.lx')
    print(absor)
    
    fig, ax = plt.subplots(1,1,figsize=(12,6))
    #y *= 1E9
    #y2 *= 1E12
    ax.plot(absor[:,0],absor[:,1]/max(absor[:,1]),   color='blue', lw=2, label='Absorption')
    ax.plot(fluor[:,0],fluor[:,1]/max(fluor[:,1]),   color='red',  lw=2, label='Fluorescence')
    ax.plot(phosph[:,0],phosph[:,1]/max(phosph[:,1]), color='green',  lw=2, label='Phosphorescence $T_1$')
    ax.plot(phosph2[:,0],phosph2[:,1]/max(phosph2[:,1]), color='orange',  lw=2, label='Phosphorescence $T_2$')

    ##ax[0].set_xlabel('Wavelength (nm)')
    ##ax[1].set_xlabel('Wavelength (nm)')
    #ax[0].set_xlabel('Energy (eV)')
    #ax[1].set_xlabel('Energy (eV)')
    #ax[2].set_xlabel('Energy (eV)')
    #ax[0].set_ylabel(r'Cross Section ($\AA^2$)')
    #ax[1].set_ylabel(r'Diff. rate')
    #ax[2].set_ylabel(r'Diff. rate')
    #
    ax.legend(loc='best')
    #ax[1].legend(loc='best')
    #ax[2].legend(loc='best')
    #
    #ax[0].set_xticks(np.round(np.linspace(2.4,5,5),1))
    #ax[1].set_xticks(np.round(np.linspace(2.0,3.4,5),1))
    #ax[2].set_xticks(np.round(np.linspace(1.3,3.0,5),1))
    ##ax[0].set_xlim([350,750])
    ##ax[1].set_xlim([350,750])
    #ax[0].set_yticks([0,5,10,15,20,25])
    #ax[0].set_ylim([0,1.2*max(absor[:,1])])
    #ax[1].set_ylim([0,1.2*max(fluor[:,1])])
    #ax[2].set_ylim([0,1.2*max(phosph[:,1])])
    ##ax[1].set_yticks([0,1,2,3,4])
    #
    #ax[0].set_title('a')
    #ax[1].set_title('b')
    #ax[2].set_title('c')
    plt.tight_layout()
    plt.savefig('overlap.png')
    #plt.show()


#for tipo in tipos:
#     isc(tipo,1)
#     mlj(tipo)
#     avgmlj(tipo)
#print('------------------------')   
#abs_spectra()
#emi_spectra()   
#phosph(2)
#SOC()    
#gaps()    

#kinetics(2)    #
#plot(2)        #x,R)
spectra()    
#painel()    
    

overlap()    
    
    
    
    
    