# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 14:56:13 2015
Reference:
1.Plant, R. E., and M. Kim. "Mathematical description of a bursting pacemaker neuron by a modification of the Hodgkin-Huxley equations." Biophysical journal 16.3 (1976): 227-244.
2. Canavier, C. C., J. W. Clark, and J. H. Byrne. "Routes to chaos in a model of a bursting neuron." Biophysical journal 57.6 (1990): 1245-1251.
kesheng.xu@cinv.cl
@author: ksxuu
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
import os.path
from matplotlib import colors

from matplotlib import rc, cm
import matplotlib.gridspec as gridspec
#for  MPI
#from mpi4py import MPI  
#comm=MPI.COMM_WORLD
#numproc=comm.size
#rank = comm.Get_rank()

def run_kut4(F,t,y,dt,args1):
        K0 = dt*F(y,t,args1)
        K1 = dt*F(y + K0/2.0,t + dt/2.0, args1)
        K2 = dt*F( y + K1/2.0,t + dt/2.0,args1)
        K3 = dt*F( y + K2,t + dt,args1)
        return (K0 + 2.0*K1 + 2.0*K2 + K3)/6.0


def PLANT(Var,t):
    [v,m,n,h,q,r,s]=Var;
    
    I_Na= (gna*m**3*h+gt)*(v-vna);
    I_K = gk*n**4*(v-vk);
    I_A = ga*q*r*(v-vk);
    I_ks = gp*s*(v-vk);
    I_l = gl*(v-vl);
    
                                                                            
    Imemb=I_Na + I_K + I_A + I_ks + I_l;
    
    alpha_m = 0.1*(-26-1.2*v)/(np.exp((-26-1.2*v)/10)-1);
    beta_m = 4*np.exp((-51-1.21*v)/18);
    alpha_n = 0.01*(-21-1.21*v)/(np.exp((-21-1.21*v)/10)-1);
    beta_n = 0.125*np.exp((-31-v)/80);
    alpha_h = 0.07*np.exp((-51-1.21*v)/20);
    beta_h = 1/(np.exp((-21-1.3*v)/10)+1);
                
    minf = alpha_m/(alpha_m+beta_m);                                                              
    ninf = alpha_n/(alpha_n+beta_n); 
    hinf = alpha_h/(alpha_h+beta_h);
    qinf = 1/(1+np.exp(-0.08*(v+45)));
    rinf = 1/(1+np.exp(0.27*(v+50)));
    sinf = 1/(1+np.exp(-0.7*(v+47)));
    
    tau_m = 12.5/(alpha_m+beta_m);
    tau_n = 12.5/(alpha_n+beta_n);
    tau_h = 12.5/(alpha_h+beta_h);
                                                               
    Det=np.array([-Imemb+Iep+I_stim,
                (minf - m)/tau_m,
                (ninf - n)/tau_n,
                (hinf - h)/tau_h,
                (qinf - q)/tau_q,
                (rinf - r)/tau_r,
                (sinf - s)/tau_s])
    return Det


"""
The main function is starting from here          
"""
# The default value of the plant model
gna = 4.0; gk = 0.3; gl=0.003;
gt=0.008; ga=0.06;  gp=0.015;
vna = 30; vk = -75; vl = -40;
tau_q=10;  tau_r=235; tau_s=8000;
Iep=-0.22; 
I_stim=0.1865  # stimulus that you can change according to the paper


# initial value 
v=-40 #Fijamos arbitrariamente un voltaje inicial
t_trans=100000
t_simu=100000
delta_t=0.01


alpha_m=0.1*(-26-1.2*v)/(np.exp((-26-1.2*v)/10)-1);
beta_m=4*np.exp((-51-1.21*v)/18);

alpha_n=0.01*(-21-1.21*v)/(np.exp((-21-1.21*v)/10)-1);
beta_n=0.125*np.exp((-31-v)/80);

alpha_h=0.07*np.exp((-51-1.21*v)/20);
beta_h=1/(np.exp((-21-1.3*v)/10)+1);

m=alpha_m/(alpha_m+beta_m);
n=alpha_n/(alpha_n+beta_n);
h=alpha_h/(alpha_h+beta_h);
q=1/(1+np.exp(-0.08*(v+45)));
r=1/(1+np.exp(0.27*(v+50)));
s=1/(1+np.exp(-0.7*(v+47)));



#initial conditions
y0=[v,m,n,h,q,r,s]
# integrate to get rid of transient behaviour:
time = np.arange(0,t_trans,delta_t)
time1 = np.arange(0,t_simu,delta_t)
Yt=integrate.odeint(PLANT, y0,time,args = ())
Y0=Yt[-1,:]
Y=integrate.odeint(PLANT, Y0,time1,args = ())



fig=plt.figure(1,figsize=(6,4))
plt.subplot(1,2,1)
plt.plot(time1,Y[:,0])
plt.subplot(1,2,2)
plt.plot(Y[:,6],Y[:,0])

plt.savefig('timeseries_%s.TIF'%I_stim)

plt.tight_layout()

#t_max=10000
#N_steps =int(t_max/delta_t);   # number of steps required
#y = Y[-1,:]
##   # integrate both orbits by a time delta_t:
#for I in range(N_steps):
#    y = y + run_kut4(HyB,t,y,dt,tempF)
#    t=t+delta_t
   





#plt.plot(Y[:,0],Y[:,2])
