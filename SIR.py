#!/usr/bin/env python
# coding: utf-8

# In[185]:


import matplotlib.pyplot as plt
import statistics as stats
import pandas as pd
import numpy as np
import math
import random

class BetaSimulation:
    
    def __init__(self, i, demographic_dist, interactions_range, nIterations):
        self.initial_infected_pop = i
        self.transmission_probability_nomask = 0.174
        self.transmission_probability_mask = 0.031
        
        # extracts from array containing the respective age group percentages (c = children, t = teen, a = adult, s = senior)
        self.c_percent = demographic_dist[0]
        self.t_percent = demographic_dist[1]
        self.a_percent = demographic_dist[2]
        self.s_percent = demographic_dist[3]
        
        # same as above but for the range of interactions per day
        self.c_range = interactions_range[0]
        self.t_range = interactions_range[1]
        self.a_range = interactions_range[2]
        self.s_range = interactions_range[3]
        
        # how many simulations that are ran
        self.iterations = nIterations
        self.results = []
    
    def run(self):
        for i in range(self.iterations):
            # sum of all beta values by demographic
            result = self.getExpectedBeta(self.c_percent, self.c_range[0], self.c_range[1]) 
            result += self.getExpectedBeta(self.t_percent, self.t_range[0], self.t_range[1]) 
            result += self.getExpectedBeta(self.a_percent, self.a_range[0], self.a_range[1]) 
            result += self.getExpectedBeta(self.s_percent, self.s_range[0], self.s_range[1])
            self.results.append(result)
        
        # total interactions multiplied by transmission probability with masking
        return (stats.mean(self.results) * self.transmission_probability_nomask)
    
    def getExpectedBeta(self, percent, min, max):
        return (self.initial_infected_pop*percent)*random.uniform(min,max)

# Denmark used as an example
population = 5813302
ndays = 135
dt = 1 #time step in days
gamma = (1.0/14.0) #recovery rate - assumed to be constant

S = np.zeros(ndays) #susceptible population
I = np.zeros(ndays) #infected population
R = np.zeros(ndays) #removed population
t = np.arange(ndays)*dt


I[0] = 1/population
S[0] = 1 - I[0]
R[0] = 0.00

# gets beta value by using BetaSimulation class
beta = BetaSimulation(1, [0.1642, 0.1233, 0.5134, 0.1991], [[0.0,2.0],[0.0,3.0],[0.0,5.0],[0.0,2.0]], 1000)
b = beta.run()

for i in range(ndays-1):
    S[i+1] = S[i] - b*S[i]*I[i]*dt
    I[i+1] = I[i] + (b*S[i]*I[i]-gamma*I[i])*dt
    R[i+1] = R[i] + gamma*I[i]*dt
    
    
fig2 = plt.figure(2)
fig2.clf()
                     
#plt.plot(t,S,'r',lw=3, label='Susceptible') - can uncomment to see complete SIR, due to large population of Denmark it has been excluded
plt.plot(t,I,'g',lw=3, label='Infected')
plt.plot(t,R,'b',lw=3, label='Removed')
fig2.legend(); plt.xlabel('Days'); plt.ylabel('Fraction of Population')


# In[123]:


arr = [0,0,0,0,0,0,0,0,0,0,0,1,2,10,13,13,18,30,60,85,144,177,223,324,422,618,768,919,990,1057,1138,1256,1381,1516,1687,1815,1922,2093,2302,2594,2849,3182,3447,3689,4060,4487,4833,5317,5737,6024,6279,6578,7233,7772,8331,8519,8813,9180,9532,9920,10409,10752,11357,11803,12104,12480,13125,13770,14442,15135,15080,15305,15767,16229,16612,16773,16819,16564,16223,16345,16662,17018,17064,17111,16975,16800,16860,17150,17411,17589,17673,17408,16979,16877,17155,17502,17405,17215,16860,16288,16305,16751,16987,17019,16993,16789,16522,16871,17510,18079,18460,18861,18856,18617,18785,19443,20425,20964,21595,21970,21700,22027,22977,24004,24832,25233,25321,25104,25502,26601,27809,28342,28746,28727,28369,28663,29206,29243,29030,28678,28003,27172,26643,26459,26315,25712,24644,23659,22436,21574,21468,21070,20129,18956,17653,16554,15987,15898,15394,14305,12869,11728,10566,9882,9750,9325,8819,8392,8008,7352,7153,7171,7345,7445,7292,7218,6922,6810,7121,7394,7445,7502,7460,7239,7222, 7426,7646,7753,7754,7694,7489,7525,7705,7878,7797,7675,7504,7294,7153,7286,7334,7287,7124,6917,6604,6529,6692,6810,6647,6495,6338,6100,6094,6323]
length = len(arr)
print(length)
print(arr)


# In[ ]:




