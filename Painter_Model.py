#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 15:15:31 2020

@author: amithzafal
"""


import numpy as np
from joblib import Parallel, delayed
from numba import jit

@jit(nopython=True)           #Upadation Matrix: The function updates the populations according to reaction happening 
def updation(num,population,N):
            c0        = population[0].copy()
            ci        = population[1].copy()
        
            l   = num
            if l < N:               # 0 to  -1
               c0[l] = 0
               ci[l] = 1
        
            elif l < 2*N:  # 0 to  1
                ci[l-N] = 0
                c0[l-N] = 1
                #w_check = 2
            elif l ==2*N:               # 0 to  -1
                population[2,0] =  population[2,0] + 1#mRNA birth
            elif l ==2*N+1:               # 0 to  -1
                population[2,0] =  population[2,0] - 1 #mRNA decay    
            elif l ==2*N+2:               # 0 to  -1
                population[2,1] =  population[2,1] + 1 #Translation 
               #population[0,0] = population[0,0] - 1 #mRNA translated
            elif l ==2*N+3:               # 0 to  -1
               population[2,2] = population[2,2] + 1 #GFP Maturation
               population[2,1] = population[2,1] - 1
            elif l ==2*N+4:               # 0 to  -1
               population[2,1] = population[2,1] - 1 #GFP decay      
            elif l ==2*N+5:               # 0 to  -1
               population[2,2] = population[2,2] - 1 #Mature GFP decay     
            
            population[0,:] = c0.copy()
            population[1,:] = ci.copy()  
            
            return population
        
        
        
@jit(nopython=True)            
def sample_discrete(probs):              
             
            q = np.random.rand()
            i = 0
            p_sum = 0
            while p_sum < q:
                p_sum = p_sum + probs[i]
                i = i + 1
            return i - 1
        
@jit(nopython=True)             
def Gillespie_Polycomb(population,time):
        
            props = PolycombPropensity(population,time)
            props_sum = props.sum()
            #Time
            time = np.random.exponential(1.0 / props_sum)
            rxn_probs = props / props_sum
            rxn_probs = rxn_probs.flatten()
            rxn      = sample_discrete(rxn_probs)
            return rxn, time
        
@jit(nopython=True)
def RandAssign1D(N):   
        
            permut = np.array([0,-1])
            nucconfig = np.zeros(N, dtype=np.float64)
            for i in range(0,N):
                nucconfig[i] = np.random.choice(permut)
            return nucconfig

@jit(nopython=True)   
def PolycombPropensity(population,time): #nucconf- c1,c0,ci
            
            c0        = population[0].copy()
            ci        = population[1].copy()
            propi0      = np.zeros(N, dtype=np.float64)
            prop0i      = np.zeros(N, dtype=np.float64)
            props      =  np.zeros((3,N), dtype=np.float64)
            propi0      = Piu.copy()                               #-1 to  0 # I to U
            rs          = np.zeros(N, dtype=np.float64)            #Enhanced Painting
            dels        = np.zeros(N, dtype=np.float64)     
            if time <= ts:
                propi0      = Piu.copy()                           #-1 to  0 # I to U
                if r !=0 or delta !=0:
                    for i in range(0,N):
                        for j in range(0,N):
                            if ci[j] == 1 and j != i:
                                #MMij[i][j] = 1/(abs(i-j)**gamma)
                                rs[i] = rs[i] + Oij[j][i]
                                dels[i] = dels[i] + Mij[j][i]
                        prop0i[i]    = np.float64(Rui[i] + Kme*Eme*r*(rs[i]) + Kme*Eme*delta*(dels[i]) + Kme*Eme*delta*r*(dels[i]))
                else:
                   prop0i    = Rui.copy()
            else:
                propi0      = Piu1.copy()                               #-1 to  0 # I to U
                if r !=0 or delta !=0:
                    for i in range(0,N):
                        for j in range(0,N):
                            if ci[j] == 1 and j != i:
                                #MMij[i][j] = 1/(abs(i-j)**gamma)
                                dels[i] = dels[i] + Mij[j][i]
                        prop0i[i]    = np.float64(Rui1[i] +  Kme*Eme*delta*(dels[i]) + Kme*Eme*delta*r*(dels[i]))
                else:
                    prop0i    = Rui1.copy()  
            mRNA  = population[2,0]
            GFP   = population[2,1]
            GFP_m = population[2,2]
            props[0,:] = prop0i*c0
            props[1,:] = propi0*ci
            #Transcription
            mT         = -1*(np.sum(ci[48:53])/5) #specifies the transcribing region
            aT         = a0*(np.tanh((mT-mc)/d)-np.tanh((-1-mc)/d))/(np.tanh((1-mc)/d)-np.tanh((-1-mc)/d)) + 0.2   
            props[2,0] = aT
            props[2,1] = mRNA*m_decay
            props[2,2] = mRNA*m_gfp
            props[2,3] = GFP*g_mat
            props[2,4] = GFP*g_decay
            props[2,5] = GFP_m*g_decay
            return props
        
@jit(nopython=True)
def Polycomb_Replication(population):     
    
            c0        = population[0].copy()
            ci        = population[1].copy()        
            for j in range(0,N):
                rc = np.random.uniform(0.0,1.0)
                if rc >= 0.5:
                        c0[j] = 1
                        ci[j] = 0
            population[0,:] = c0.copy()
            population[1,:] = ci.copy()
            return population

@jit(nopython=True)       
def Stocha(traject_num):

                initialconf = RandAssign1D(N)
                c0          = np.zeros(N, dtype=np.float64)
                ci          = np.zeros(N, dtype=np.float64)
                for q in range(0,N):
                    if initialconf[q] == -1:
                        ci[q] = np.float64(1)
                    else:
                        c0[q] = np.float64(1)
                outconfF        = np.zeros((snaps+100,N), dtype = np.float64)
                population      = np.zeros((3,N), dtype=np.float64)
                population[0,:] = c0.copy()
                population[1,:] = ci.copy()
                transconfig      = np.zeros((N**3,3), dtype = np.int64)
                t           = 0                              # Start time
                ii           = 0                              #Counter
                repli        = 0                              #Counter Replication
                dTcyc        = Tcyc
                deltaT       = 0         #snaps
                while t < T:
        
        
                    num, dt = Gillespie_Polycomb(population,t)
                    prevpopula    = population.copy()
                    population    = updation(num,population,N)
                    t             = t + dt
                    config        = prevpopula[0]*0 + prevpopula[1]*-1 
                    if t >= dTcyc:                       # DNA Replication
                        t             = dTcyc 
                        repli         = repli + 1
                        dTcyc         = dTcyc + Tcyc
                        population    = Polycomb_Replication(prevpopula)
                        config        = population[0]*0 + population[1]*(-1)
                    while t > deltaT:
                        deltaT              = deltaT + delT
                        outconfF[ii]        = config.copy()
                        transconfig[ii]     = prevpopula[2,0:3]
                        ii                  = ii + 1
                return transconfig[0:snaps]                              #Variable-outconfF saves Epigenetic configurations as a function of time and nucleosome position. 
                                                                          #Variable transconfig saves mRNA, GFP,GFP* levels as a fucntion of time       
########################################################### INPUT Parsed
### STANDARD CODE PARAMETERS
N           = 100               #Total Number of Nucleosomes.
T           = 300               #End Time
Tcyc        = 50                #Replication time
############################################################### MODEL PARAMETERS
### PRE/TRE Profiles
M     = np.zeros(N, dtype=np.float64)       # a function of j; j [0,N] #PRC2 Methylation Protein
#M     = M  + 0.01                              # Background
DM    = np.zeros(N, dtype=np.float64)       # a function of k; k [0,N] #Demethylation Enzyme
DM    = DM + 0                              # Background
#M1    = np.zeros(N, dtype=np.float64)       # a function of k; k [0,N] #Demethylation Enzyme
#M1    = M1 + 0.01 

i_start = 48                     #Nucleosome Position where the enzymes are localized
i_end   = 52
########################################################## Collapsed PARAMETERS
K0          = 0.1              # Nucleosome turnover
Kme         = 2*K0             # Acetylation, Methylation rate
Kdme        = 0                # Demethylation rate
Khdac       = 0                # Deacetylation Rate
Eme         = .6               # Epsilon - spreading efficiency
Edme        = 0                # Demethylation Spreading efficiency
r           = 0                # Allosteric boost
delta       = 0.0              # Reader Writer Mechanism
gamma       = 1                # Contact Probability |i - j| exponent
delT        = 1                # Determines the time interval between recording epigenetic configuration state. (Time resolution) 
snaps       = np.int64(T/delT) # Number of timepoints
limt        = 2400             # Number of Trajectories or cells 
ts          = 20000            # Specifies the time at which painter unbinds    
#TRANSCRIPTION
a0          =  1               # Maximal transcription rate. Actual would be higher considering leaky rate 
d           =  0.05            # Steepness of transcriptional switch  
mc          = -0.4             # Critical magenetization value at which the switch happens
mRNA        = 0
GFP         = 0
GFP_m       = 0
m_decay     = 1/10
g_decay     = 1/30
m_gfp       = 1
g_mat       = 1/10 
for i in range(i_start,i_end+1):
    M[i]  =  1    
Mij  = np.zeros((N, N))
Lij  = np.zeros((N, N))
Piu  = np.zeros(N, dtype=np.float64)
Piu1  = np.zeros(N, dtype=np.float64)
Oij  = np.zeros((N, N))
Oij1  = np.zeros((N, N))
Rui  = np.zeros(N, dtype=np.float64)                              
Rui1  = np.zeros(N, dtype=np.float64)
for k in range(0,N):
    for l in range(0,N):
        if l == k:
            Mij[l][k] = 0
            Lij[l][k] = 0
            Oij[l][k] = 0
            Oij1[l][k] = 0
        else:
            Mij[l][k]  = 1/(abs(k-l)**gamma)
            Lij[l][k]  = DM[l]*Mij[l][k]
            Oij[l][k]  = M[l]*Mij[l][k]
#            Oij1[l][k] = M1[l]*Mij[l][k]
mij_sum = Mij.sum(axis=1)
lij_sum = Lij.sum(axis=0)
oij_sum = Oij.sum(axis=0)
oij_sum1 = Oij1.sum(axis=0)
for i in range(0,N):
                Piu[i] = np.float64(K0 + Kdme*(DM[i] + Edme*(lij_sum[i])))
                Rui[i] = np.float64(Kme*M[i] + Kme*Eme*oij_sum[i]) 
                Piu1[i] = np.float64(K0 + Kdme*(DM[i] + Edme*(lij_sum[i])))
                Rui1[i] = np.float64(0)       
                                  
fin_out = Parallel(n_jobs=-1)(delayed(Stocha)(iii)for iii in range(0,limt))                
np.savetxt('output.txt',fin_out[1])

