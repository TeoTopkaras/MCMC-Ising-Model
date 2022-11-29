import numpy as np
import matplotlib.pyplot as plt
from numpy import random
#-------------------------------------------------------------
#Full Ising model simulation
#-------------------------------------------------------------
def ising_mc_simu(n_particles,ext_field,coupl_const,temp_min,temp_max,temp_step,n_sys,n_training,n_eq_lim):
    '''Main part of the Monte Carlo simulation of the Ising Model. User defined input parameters and simulation constants
        are specified in the main_file. Returns lists of parameters for the correlation time, average magnetisation and energy
        per spin, with corresponding uncertainties, magnetisation susceptibility per spin and specific heat per spin. 
    '''
    temp,corr_t,magn_spin,sigma_magn_spin,energy_spin,sigma_energy_spin,magn_sus,spec_heat=init_arrays(temp_min,temp_max,temp_step,n_sys)
    for i in range(0,len(temp)):
    	m_exp,s_m_exp=training_sample(n_particles,ext_field,coupl_const,temp[i],n_training,n_eq_lim)
    	print('Expected Magnetisation per spin at Temperature T='+str(temp[i])+' : '+str(m_exp)+'+/- '+str(s_m_exp))
    	for j in range(0,n_sys):
    		sys=sys_eq(n_particles,ext_field,coupl_const,temp[i],n_eq_lim,m_exp,s_m_exp)
    		corr_t[i][j]=corr_time(sys,n_particles,ext_field,coupl_const,temp[i])
    		magn_spin[i][j],sigma_magn_spin[i][j],energy_spin[i][j],sigma_energy_spin[i][j],magn_sus[i][j],spec_heat[i][j]=block_measurement(sys,corr_t[i][j],n_particles,ext_field,coupl_const,temp[i])
        
    return(temp,corr_t,magn_spin,sigma_magn_spin,energy_spin,sigma_energy_spin,magn_sus,spec_heat)
#-------------------------------------------------------------
#Initialisation functions
#-------------------------------------------------------------
def init_arrays(temp_min,temp_max,temp_step,n_sys):
    '''Initialisation function of the necessary arrays
    Parameters
    ---------------
    temp_min : float
            Lower temperature limit
    temp_max: float
            Upper temperature limit
    temp_step: float
            Temperature increment 
    n_sys: int
            Number of systems simulated per temperature
    
    '''
    temp= np.round(np.arange(temp_min, temp_max, temp_step),1)

    corr_t=[[0. for i in range(n_sys)] for j in range(len(temp))]
    magn_spin=[[0. for i in range(n_sys)] for j in range(len(temp))]
    sigma_magn_spin=[[0. for i in range(n_sys)] for j in range(len(temp))]
    energy_spin=[[0. for i in range(n_sys)] for j in range(len(temp))]
    sigma_energy_spin=[[0. for i in range(n_sys)] for j in range(len(temp))]
    magn_sus=[[0. for i in range(n_sys)] for j in range(len(temp))]
    spec_heat=[[0. for i in range(n_sys)] for j in range(len(temp))]
    
    return(temp,corr_t,magn_spin,sigma_magn_spin,energy_spin,sigma_energy_spin,magn_sus,spec_heat)
#-------------------------------------------------------------
#Algorithms
#-------------------------------------------------------------
def ising_mcmc(sys,n_particles,ext_field,coupl_const,temp,sampling_ratio):
    '''Performs the evolution of a system following a Monte Carlo Markov Chain
    Parameters
    ---------------
    sys : array of size (N,N)
        Configuration of the system's spin orientation
    n_particles : int
        Number of particles
    ext_field : float
        Value of the external magnetic field
    coupl_const : float
        Value of the coupling constant of the system
    temp : float
        Temperature of the system 
    sampling_ratio : int
        Number of iterations (i.e. spin flips) performed 
    
    Returns
    ---------------
    sys : array of size (N,N)
        New configuration of the system's spin orientation
    
    '''
    for t in range(0,sampling_ratio):
        i=np.random.randint(n_particles)
        j=np.random.randint(n_particles)
        
        Sn=sys[(i+1)%n_particles][j]+sys[(i-1)%n_particles][j]+sys[i][(j+1)%n_particles]+sys[i][(j-1)%n_particles]
        dE=2*sys[i][j]*Sn*coupl_const+2*sys[i][j]*ext_field
        
        p_ratio=np.exp(-dE/temp)
        p_test=random.random()
        
        if ((p_test<p_ratio)):
            sys[i][j]=-sys[i][j]
    return(sys)

def training_sample(n_particles,ext_field,coupl_const,temp,n_training,n_eq_lim):
    ''' Determines the expected magnetisation per spin at equilibrium for a given temperature 
    Parameters
    ---------------
    n_particles : int
        Number of particles
    ext_field : float
        Value of the external magnetic field
    coupl_const : float
        Value of the coupling constant of the system
    temp : float
        Temperature of the system 
    n_training : int
        Number of systems to be evolved for the determination of the expected average megnitisation per spin
    n_eq_lim : int
        Number of iterations after which the system is expected to be reaching equilibrium
    
    Returns
    ---------------
    np.mean(m_samp) : float
        Expected magnetisation per spin at equilibrium for a given temperature
        
    np.std(m_samp) : float
        Standard deviation of the expected magnetisation per spin at equilibrium for a given temperature
        
    '''
    
    magn_samp=[]
    for n in range(0,n_training):
        magn=[]
        sys=[[np.random.choice([-1,1], p=[0.2,0.8]) for i in range(n_particles)] for i in range(n_particles)]
        magn.append(np.sum(sys)/n_particles**2)
        for t in range(0,500):
            sys=ising_mcmc(sys,n_particles,ext_field,coupl_const,temp,n_particles**2)
            magn.append(np.sum(sys)/n_particles**2.)
            if (len(magn)>=n_eq_lim) :
                sigma=((np.mean(magn[len(magn)-11:len(magn)-1])-magn[len(magn)-1])/np.std(magn[len(magn)-11:len(magn)-1]))
                if np.abs(sigma)<1 :
                    magn_samp.append(magn[len(magn)-1])
                    break
    return(np.mean(magn_samp),np.std(magn_samp))

def sys_eq(n_particles,ext_field,coupl_const,temp,n_eq_lim,magn_exp,sigma_magn_exp):
    ''' Returns a system configuration representative of a system in equilibrium
    Parameters
    ---------------
    n_particles : int
        Number of particles
    ext_field : float
        Value of the external magnetic field
    coupl_const : float
        Value of the coupling constant of the system
    temp : float
        Temperature of the system 
    n_eq_lim : int
        Number of iterations after which the system is expected to be reaching equilibrium
    magn_exp : float
        Expected magnetisation per spin at equilibrium for a given temperature 
    
    sigma_magn_exp : float
        Standard deviation of the expected magnetisation per spin at equilibrium for a given temperature
        
    
    Returns
    ---------------
    sys :array of size (N,N)
        Configuration of the system's spin orientation
        
    '''
    magn=[]
    sys=[[np.random.choice([-1,1], p=[0.2,0.8]) for i in range(n_particles)] for i in range(n_particles)]
    magn.append(np.sum(sys)/n_particles**2.)
    for t in range(0,2000):
        sys=ising_mcmc(sys,n_particles,ext_field,coupl_const,temp,n_particles**2)
        magn.append(np.sum(sys)/n_particles**2.)
        if (len(magn)>=n_eq_lim):
            sigma_exp=((magn_exp-magn[len(magn)-1])/sigma_magn_exp)
            if (np.abs(sigma_exp)<1) :
                break
    return(sys)

def corr_time(sys,n_particles,ext_field,coupl_const,temp):
    ''' Estimates the system's correlation time
    Parameters
    ---------------
    sys : array of size (N,N)
        Configuration of the system's spin orientation. The system provided is assumed to be in equilibrium
    n_particles : int
        Number of particles
    ext_field : float
        Value of the external magnetic field
    coupl_const : float
        Value of the coupling constant of the system
    temp : float
        Temperature of the system 
    
    Returns
    ---------------
    corr_t : ndarray
        Correlation time of the system
        
    '''
    magn=[]
    magn.append(np.sum(sys)/n_particles**2.)
    for t in range(0,50):
        sys=ising_mcmc(sys,n_particles,ext_field,coupl_const,temp,n_particles**2)
        magn.append(np.sum(sys)/n_particles**2.)
    x_func=[]
    for i in range(0, len(magn)):
        S_tti=[magn[j]*magn[i+j] for j in range(len(magn)-i)]
        S_t=[magn[j] for j in range(len(magn)-i)]
        S_ti=[magn[i+j] for j in range(len(magn)-i)]
        x_func.append((1/(len(magn)-i)*np.sum(S_tti))-((1/(len(magn)-i))**2*np.sum(S_t)*np.sum(S_ti)))
    
    corr_t=0.
    for i in range(0, len(x_func)):
        if x_func[i]>0 :
            corr_t+=(x_func[i]/x_func[0])
        else:
            break
    return(corr_t)

def block_measurement(sys,corr_t,n_particles,ext_field,coupl_const,temp):
    '''Estimates the average magnetisation and energy per spin, the magnetic susceptibility and specific het per spin
            following the blocking method
    Parameters
    ---------------
    sys : array of size (N,N)
        Configuration of the system's spin orientation. The system provided is assumed to be in equilibrium
    corr_t : float
        Correlation time 
    n_particles : int
        Number of particles
    ext_field : float
        Value of the external magnetic field
    coupl_const : float
        Value of the coupling constant of the system
    temp : float
        Temperature of the system 
    
    Returns
    ---------------
    magn_spin : float
        Magnetisation per spin
    sigma_magn_spin : float
        Error of the magnetisation per spin  
    energy_spin : float
        Energy per spin of the system
    sigma_energy_spin : float
        Error of the energy per spin
    magn_sus : float
        Magnetic susceptibility per spin
    spec_heat : float
     Specific heat per spin of the system
    '''
        
    e_tot=[]
    m_tot=[]
    for i in range(0,16):
    	sys=ising_mcmc(sys,n_particles,ext_field,coupl_const,temp,int(corr_t*n_particles**2.))  
    	sn=[[sys[i][j]*(sys[(i+1)%n_particles][j]+sys[(i-1)%n_particles][j]+sys[i][(j+1)%n_particles]+sys[i][(j-1)%n_particles]) for i in range(n_particles)] for j in range(n_particles)]
    	e_tot.append(-(coupl_const*0.5*np.sum(sn)+ext_field*np.sum(sys)))
    	m_tot.append(np.sum(sys))
    e_tot2=[e_tot**2 for e_tot in e_tot]
    m_tot2=[m_tot**2 for m_tot in m_tot]
    
    magn_spin=np.mean(m_tot)/n_particles**2
    sigma_magn_spin=np.sqrt((np.mean(m_tot2)-np.mean(m_tot)**2.)/8.)/n_particles**2. 
    energy_spin=np.mean(e_tot)/n_particles**2
    sigma_energy_spin=np.sqrt((np.mean(e_tot2)-np.mean(e_tot)**2.)/8.)/n_particles**2.
    magn_sus=((np.mean(m_tot2)-np.mean(m_tot)**2.))/(temp*n_particles**2.)
    spec_heat=((np.mean(e_tot2)-np.mean(e_tot)**2.))/(temp*n_particles)**2.
    
    return(magn_spin,sigma_magn_spin,energy_spin,sigma_energy_spin,magn_sus,spec_heat)

def bootstrap_mean(samp,n_bootstrap):
    '''performs bootstrapping for the calculation of the average and standard deviation of the average of an array

    Parameters
    ---------------
    samp : array of size (N)
        list of the values the average has to be calculated for
    n_bootstrap: int
	the number of bootstrap iterations
    
    Returns
    ---------------
    np.mean(samp_m) : float
        Mean value of the array
    np.std(samp_m) : scalar
        Standard deviation of the mean value estimate
    '''
    samp_m=[0. for i in range(n_bootstrap)]
    for i in range(0,n_bootstrap):
        samp_it=np.random.choice(samp, size=len(samp), replace=True)
        samp_m[i]=samp_it.mean()   
    return(np.mean(samp_m),np.std(samp_m))
#-------------------------------------------------------------
#Plotting functions
#-------------------------------------------------------------
def plot_corr_time(temp,corr_t,n_bootstrap):
    '''Performs statistics and plots the correlation time as a function of temperature
    Parameters
    ---------------
    temp : array of size(N)
        List of temperatures explored 
    
    corr_t : array of size(N,S)
        List of correlation times of a system per temperature
    n_bootstrap: int
	the number of bootstrap iterations
   
    '''
    t_m=[0. for i in range(len(temp))]
    s_tm=[0. for i in range(len(temp))]
    for i in range(0,len(temp)):
        t_list=corr_t[i]
        t_m[i],s_tm[i]=bootstrap_mean(t_list,n_bootstrap)
    
    plt.errorbar(temp,t_m, yerr=s_tm, marker='s',fmt='o',color='black')
            
    plt.xlabel('T/k$_B$')
    plt.ylabel('τ/N$^2$')
    plt.show()

def plot_magn_spin(temp,magn_spin,sigma_magn_spin):
    '''Performs statistics and plots the average magnetisation per spin as a function of temperature
    Parameters
    ---------------
    temp : array of size(N)
        List of temperatures explored 
    
    magn_spin : array of size(N,S)
        List of magnetisation per spin of a system per temperature 
        
    sigma_magn_spin : array of size(N,S)
        List of errors of the magnetisation per spin of a system per temperature 
   
    '''

    mgn=[0. for i in range(len(temp))]
    s_mgn=[0. for i in range(len(temp))]
    for i in range(0,len(temp)):
        mgn[i]=np.mean(magn_spin[i])
        sigma2=[sigma_magn_spin**2. for sigma_magn_spin in sigma_magn_spin[i]]
        s_mgn[i]=np.sqrt(np.sum(sigma2))/len(magn_spin[i])
    
    plt.errorbar(temp,mgn, yerr=s_mgn,marker='s',fmt='o',color='black')


    plt.xlabel('T/k$_B$')
    plt.ylabel('<|m|>')
    plt.show()

def plot_energy_spin(temp,energy_spin,sigma_energy_spin):
    '''Performs statistics and plots the average energy per spin as a function of temperature
    Parameters
    ---------------
    temp : array of size(N)
        List of temperatures explored 
    
    energy_spin : array of size(N,S)
        List of energies per spin of a system per temperature 
        
    sigma_energy_spin : array of size(N,S)
        List of errors of the energies per spin of a system per temperature 
   
    '''
    enr=[0. for i in range(len(temp))]
    s_enr=[0. for i in range(len(temp))]
    for i in range(0,len(temp)):
        enr[i]=np.mean(energy_spin[i])
        s2=[sigma_energy_spin**2. for sigma_energy_spin in sigma_energy_spin[i]]
        s_enr[i]=np.sqrt(np.sum(s2))/len(energy_spin[i])
    
    plt.errorbar(temp,enr, yerr=s_enr,marker='s',fmt='o',color='black')


    plt.xlabel('T/k$_B$')
    plt.ylabel('e')
    plt.show()

def plot_magn_sus(temp,magn_sus,n_bootstrap):
    '''Performs statistics and plots the magnetic susceptibility per spin as a function of temperature
    Parameters
    ---------------
    temp : array of size(N)
        List of temperatures explored 
    
    magn_sus : array of size(N,S)
        List of magnetic susceptibilities per spin of a system per temperature 
    n_bootstrap: int
	the number of bootstrap iterations
    '''
    t_m=[0. for i in range(len(temp))]
    s_tm=[0. for i in range(len(temp))]
    for i in range(0,len(temp)):
        t_list=magn_sus[i]
        t_m[i],s_tm[i]=bootstrap_mean(t_list,n_bootstrap)
    
    plt.errorbar(temp,t_m, yerr=s_tm, marker='s',fmt='o',color='black')
            
    plt.xlabel('T/k$_B$')
    plt.ylabel('χ$_M$')
    plt.show()

def plot_spec_heat(temp,spec_heat,n_bootstrap):
    '''Performs statistics and plots the specific heat per spin as a function of temperature
    Parameters
    ---------------
    temp : array of size(N)
        List of temperatures explored 
    
    spec_heat : array of size(N,S)
        List of specific heat per spin of a system per temperature 
    n_bootstrap: int
	the number of bootstrap iterations
    '''
    t_m=[0. for i in range(len(temp))]
    s_tm=[0. for i in range(len(temp))]
    for i in range(0,len(temp)):
        t_list=spec_heat[i]
        t_m[i],s_tm[i]=bootstrap_mean(t_list,n_bootstrap)
    
    plt.errorbar(temp,t_m, yerr=s_tm, marker='s',fmt='o',color='black')
            
    plt.xlabel('T/k$_B$')
    plt.ylabel('C')
    plt.show()


