#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import statistics as stats
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
    
    
fig = plt.figure(0)
fig.clf()
                     
#plt.plot(t,S,'r',lw=3, label='Susceptible') - can uncomment to see complete SIR, due to large population of Denmark it has been excluded
plt.plot(t,I,'g',lw=3, label='Infected')
plt.plot(t,R,'b',lw=3, label='Removed')
fig.legend(); plt.xlabel('Days'); plt.ylabel('Fraction of Population')
