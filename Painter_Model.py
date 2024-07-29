#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 15:15:31 2020

@author: amithzafal
"""
# Copyright © 2020 ENS Lyon. All rights reserved.

import numpy as np
from joblib import Parallel, delayed
#import numba
#from numba import jit
#from numba.experimental import jitclass
import os

#For transcription, transcribing region can be specified in the function PolycombPropensity
#For enzyme limitation, specify that in the funtion PolycombPropensity (delta = ps). 
#Variable-outconfF (function Stocha) saves Epigenetic configurations as a function of time and nucleosome position. 
#Variable transconfig saves mRNA, GFP,GFP* levels as a fucntion of time       

class Painter:
    def __init__(self, outputDir):
        self.N          = int(100)               #Total Number of Nucleosomes.
        self.delT       = 1                 #Time resolution
        self.T          = 300               #End Time
        self.snaps      = np.int32(self.T / self.delT)  #Number of timepoints
        self.Tcyc       = 50                #Replication/Generation time
        self.ts          = 20000            # Specifies the time at which painter unbinds
        
        self.i_start    = 48               # specifying the painter region (start)
        self.i_end      = 52               # specifying the painter region (end)
        self.gamma      = 1                # Contact scaling exponent
        
        self.K0         = 0.1              # Nucleosome turnover
        self.Kdme       = 0                # Demethylation rate
        self.Edme       = 0                # Demethylation efficiency
        self.Kme        = 2*self.K0             # Acetylation, Methylation rate
        self.Eme        = 0.6              # Spreading efficiency
        self.r          = 0                # Allosteric boost
        self.delta      = 0.0              # Reader Writer Mechanism

        self.Ntot        = 10               # Total number of HMEs
        self.rconst      = 1                # Ratio of binding-unbinding rates of HMEs
        self.limt        = 1                # Number of Trajectories or cells 

        #TRANSCRIPTION
        self.a0          =  1               # Maximal transcription rate. 
        self.d           =  0.05            # Steepness of transcriptional switch  
        self.mc          = -0.4             # Critical magenetization value at which the switch happens
        self.mRNA        = 0
        self.GFP         = 0
        self.GFP_m       = 0
        self.m_decay     = 1/10             # mRNA decay rate 
        self.g_decay     = 1/30             # protein decay rate
        self.m_gfp       = 1                # mRNA translation rate
        self.g_mat       = 1/10             # protein maturation rate

        self.outputDir   = outputDir
        self.outFile     = os.path.join(self.outputDir, "Cell_1.res")

    def paint(self):
            M      = np.zeros(self.N, dtype=np.float64)
            DM     = np.zeros(self.N, dtype=np.float64)
            #M     = M  + 0.01 
            M[self.i_start:self.i_end+1]  =  1    
            self.Mij  = np.zeros((self.N, self.N))
            self.Lij  = np.zeros((self.N, self.N))
            self.Piu  = np.zeros(self.N, dtype=np.float64)
            self.Piu1  = np.zeros(self.N, dtype=np.float64)
            self.Oij  = np.zeros((self.N, self.N))
            self.Oij1  = np.zeros((self.N, self.N))
            self.Rui  = np.zeros(self.N, dtype=np.float64)                              
            self.Rui1  = np.zeros(self.N, dtype=np.float64)
            for k in range(0,self.N):
                for l in range(0,self.N):
                    if l == k:
                        self.Mij[l][k] = 0
                        self.Lij[l][k] = 0
                        self.Oij[l][k] = 0
                        self.Oij1[l][k] = 0
                    else:
                        self.Mij[l][k]  = 1/(abs(k-l)**self.gamma)
                        self.Lij[l][k]  = DM[l]*self.Mij[l][k]
                        self.Oij[l][k]  = M[l]*self.Mij[l][k]
            #            Oij1[l][k] = M1[l]*Mij[l][k]
            mij_sum = self.Mij.sum(axis=1)
            lij_sum = self.Lij.sum(axis=0)
            oij_sum = self.Oij.sum(axis=0)
            oij_sum1 = self.Oij1.sum(axis=0)
            for i in range(0,self.N):
                            self.Piu[i] = np.float64(self.K0 + self.Kdme*(DM[i] + self.Edme*(lij_sum[i])))
                            self.Rui[i] = np.float64(self.Kme*M[i] + self.Kme*self.Eme*oij_sum[i]) 
                            self.Piu1[i] = np.float64(self.K0 + self.Kdme*(DM[i] + self.Edme*(lij_sum[i])))
                            self.Rui1[i] = np.float64(0)       
                                            
            self.fin_out = Parallel(n_jobs=-1)(delayed(self.Stocha)(iii)for iii in range(0,self.limt))                
            np.savetxt(self.outFile, self.fin_out[0]) #Saves the first trajectory in a text file   

    
    def updation(self,num,population):
                c0        = population[0].copy()
                ci        = population[1].copy()
        
                l   = num
                if l < self.N:               # 0 to  -1
                    c0[l] = 0
                    ci[l] = 1
        
                elif l < 2*self.N:  # 0 to  1
                    ci[l-self.N] = 0
                    c0[l-self.N] = 1
                    #w_check = 2
                elif l ==2*self.N:               # 0 to  -1
                    population[2,0] =  population[2,0] + 1#mRNA birth
                elif l ==2*self.N+1:               # 0 to  -1
                    population[2,0] =  population[2,0] - 1 #mRNA decay    
                elif l ==2*self.N+2:               # 0 to  -1
                    population[2,1] =  population[2,1] + 1 #Translation 
                    #population[0,0] = population[0,0] - 1 #mRNA translated
                elif l ==2*self.N+3:               # 0 to  -1
                    population[2,2] = population[2,2] + 1 #GFP Maturation
                    population[2,1] = population[2,1] - 1
                elif l ==2*self.N+4:               # 0 to  -1
                    population[2,1] = population[2,1] - 1 #GFP decay      
                elif l ==2*self.N+5:               # 0 to  -1
                    population[2,2] = population[2,2] - 1 #Mature GFP decay     
            
                population[0,:] = c0.copy()
                population[1,:] = ci.copy()  
            
                return population
              
            
    def sample_discrete(self,probs):              
             
            q = np.random.rand()
            i = 0
            p_sum = 0
            while p_sum < q:
                p_sum = p_sum + probs[i]
                i = i + 1
            return i - 1
        
           
    def Gillespie_Polycomb(self,population,time):
        
            props = self.PolycombPropensity(population,time)
            props_sum = props.sum()
            #Time
            time = np.random.exponential(1.0 / props_sum)
            rxn_probs = props / props_sum
            rxn_probs = rxn_probs.flatten()
            rxn      = self.sample_discrete(rxn_probs)
            return rxn, time
        

    def RandAssign1D(self,N):   
        
            permut = np.array([0,-1])
            nucconfig = np.zeros(N, dtype=np.float64)
            for i in range(0,N):
                nucconfig[i] = np.random.choice(permut)
            return nucconfig

  
    def PolycombPropensity(self,population,time): 
            
            c0        = population[0].copy()
            ci        = population[1].copy()
            propi0      = np.zeros(self.N, dtype=np.float64)
            prop0i      = np.zeros(self.N, dtype=np.float64)
            props      =  np.zeros((3,self.N), dtype=np.float64)
            propi0      = self.Piu.copy()                                   
            rs          = np.zeros(self.N, dtype=np.float64)          
            dels        = np.zeros(self.N, dtype=np.float64)
            Nsite       = np.sum(ci)
            alpha       = self.Ntot 
    
            if Nsite > 0:
                        ps          = ((1 + Nsite/alpha + self.rconst/alpha) -  np.sqrt((1 + Nsite/alpha + self.rconst/alpha)**2 - 4*Nsite/alpha))/(2*(Nsite/alpha))
            else:
                        ps = 0
            #delta      = ps        #For enzyme limitation 
            if time <= self.ts:
                propi0      = self.Piu.copy()                           #-1 to  0 # I to U
                if self.r !=0 or self.delta !=0:
                    for i in range(0,self.N):
                        for j in range(0,self.N):
                            if ci[j] == 1 and j != i:
                                #MMij[i][j] = 1/(abs(i-j)**gamma)
                                rs[i] = rs[i] + self.Oij[j][i]
                                dels[i] = dels[i] + self.Mij[j][i]
                        prop0i[i]    = np.float64(self.Rui[i] + self.Kme*self.Eme*self.r*(rs[i]) + self.Kme*self.Eme*self.delta*(dels[i]) + self.Kme*self.Eme*self.delta*self.r*(dels[i]))
                else:
                   prop0i    = self.Rui.copy()
            else:
                propi0      = self.Piu1.copy()                               #-1 to  0 # I to U
                if self.r !=0 or self.delta !=0:
                    for i in range(0,self.N):
                        for j in range(0,self.N):
                            if ci[j] == 1 and j != i:
                                #MMij[i][j] = 1/(abs(i-j)**gamma)
                                dels[i] = dels[i] + self.Mij[j][i]
                        prop0i[i]    = np.float64(self.Rui1[i] +  self.Kme*self.Eme*self.delta*(dels[i]) + self.Kme*self.Eme*self.delta*self.r*(dels[i]))
                else:
                    prop0i    = self.Rui1.copy()  
            mRNA  = population[2,0]
            GFP   = population[2,1]
            GFP_m = population[2,2]
            props[0,:] = prop0i*c0
            props[1,:] = propi0*ci

            #Transcription
            mT         = -1*(np.sum(ci[48:53])/5) #specifies the transcribing region
            aT         = self.a0*(np.tanh((mT-self.mc)/self.d)-np.tanh((-1-self.mc)/self.d))/(np.tanh((1-self.mc)/self.d)-np.tanh((-1-self.mc)/self.d)) + 0.2   
            props[2,0] = aT
            props[2,1] = mRNA*self.m_decay
            props[2,2] = mRNA*self.m_gfp
            props[2,3] = GFP*self.g_mat
            props[2,4] = GFP*self.g_decay
            props[2,5] = GFP_m*self.g_decay
            return props
        
    
    def Polycomb_Replication(self,population):     
    
            c0        = population[0].copy()
            ci        = population[1].copy()        
            for j in range(0,self.N):
                rc = np.random.uniform(0.0,1.0)
                if rc >= 0.5:
                        c0[j] = 1
                        ci[j] = 0
            population[0,:] = c0.copy()
            population[1,:] = ci.copy()
            return population

        
    def Stocha(self,traject):

                initialconf = self.RandAssign1D(self.N)
                c0          = np.zeros(self.N, dtype=np.float64)
                ci          = np.zeros(self.N, dtype=np.float64)
                for q in range(0,self.N):
                    if initialconf[q] == -1:
                        ci[q] = np.float64(1)
                    else:
                        c0[q] = np.float64(1)
                outconfF        = np.zeros((self.snaps+100,self.N), dtype = np.float64)
                population      = np.zeros((3,self.N), dtype=np.float64)
                population[0,:] = c0.copy()
                population[1,:] = ci.copy()
                transconfig      = np.zeros((self.N**3,3), dtype = np.int64)
                t           = 0                              # Start time
                ii           = 0                              #Counter
                repli        = 0                              #Counter Replication
                dTcyc        = self.Tcyc
                deltaT       = 0         #snaps
                while t < self.T:
        
        
                    num, dt = self.Gillespie_Polycomb(population,t)
                    prevpopula    = population.copy()
                    population    = self.updation(num,population)
                    t             = t + dt
                    config        = prevpopula[0]*0 + prevpopula[1]*-1 
                    if t >= dTcyc:                       # DNA Replication
                        t             = dTcyc 
                        repli         = repli + 1
                        dTcyc         = dTcyc + self.Tcyc
                        population    = self.Polycomb_Replication(prevpopula)
                        config        = population[0]*0 + population[1]*(-1)
                    while t > deltaT:
                        deltaT              = deltaT + self.delT
                        outconfF[ii]        = config.copy()
                        transconfig[ii]     = prevpopula[2,0:3]
                        ii                  = ii + 1
                return transconfig[0:self.snaps]                         #Variable-outconfF saves Epigenetic configurations as a function of time and nucleosome position. 
                                                                         #Variable transconfig saves mRNA, GFP,GFP* levels as a fucntion of time       


sim = Painter(outputDir='/Users/amith/Documents/')
sim.paint()
