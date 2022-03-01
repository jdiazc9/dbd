# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:28:37 2020
@author: Juan
"""



# %%
### INITIALIZATION

# reset variables etc.
from IPython import get_ipython
get_ipython().magic('reset -sf')

# libraries
import community_simulator
import community_simulator.usertools
import community_simulator.visualization
import numpy as np
import pandas as pd
# import matplotlib
#import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette()
# %matplotlib inline

#import random
#import math
#import copy
import statsmodels.api as sm

# check if running on Windows (if so, set parallel=False in Community object)
import os
if os.name == 'nt':
    par = False
else:
    par = True
    


# %%
### USER-DEFINED FUNCTIONS

# function to sample species from each family's pool
def sampleSpecies(how_many_families,how_many_species,assumptions,mode='random'):
    
    # 'mode' controls whether families are picked randomly or in order,
    # if mode = 'sequential' the first n families will be picked,
    # if mode = 'random' n families will be picked at random
    # (n meaning how_many_families)
    
    # total number of families and species, etc
    n_families = len(assumptions['SA'])
    n_species = sum(assumptions['SA'])
    sp_ranges = np.concatenate(([0],np.cumsum(assumptions['SA'])))
    
    # randomly choose families
    if mode == 'random':
        f = np.random.choice(range(n_families),size=how_many_families,replace=False)
    elif mode =='sequential':
        f = np.array([i for i in range(how_many_families)])
    species_to_sample = [sp_i for f_i in f for sp_i in list(range(sp_ranges[f_i],sp_ranges[f_i+1]))]
    sampled_species = np.random.choice(species_to_sample,size=how_many_species*how_many_families,replace=False)
    
    # set abundances of sampled species to 1, all others to 0
    N = np.zeros((n_species,1))
    N[sampled_species,:] = 1
    
    return N

# stabilize communities through serial passaging
def stabilizeCommunities(plate):
    
    # parameters
    dilution_factor = 1/100
    n_transfers = 20
    time_per_cycle = 1
    # extinction_threshold = 1e-4 # we don't use this because we don't use StadyState()
    
    # run stabilization protocol
    Ntraj, Rtraj = plate.RunExperiment(f = np.diag(np.ones(plate.n_wells))*dilution_factor,
                                       T = time_per_cycle,
                                       npass = n_transfers,
                                       compress_resources=False,
                                       compress_species=True)
    
    return([plate.N,plate.R])



# %%

# number of simulations to run
nSimul = 100

# initialize parameter arrays
nFamilies = np.random.randint(3,8 + 1,size=nSimul)
nSpeciesPerFamily = np.random.randint(150,300 + 1,size=nSimul)
nSampledSpeciesPerFamily = np.random.randint(20,100 + 1,size=nSimul)
nResourcesPerClass = np.random.randint(3,15 + 1,size=nSimul)
leakageFraction = np.random.uniform(0.5,0.8,size=nSimul)
preferenceStrength = np.random.uniform(0.0,1.0,size=nSimul)
byproductSparsity = 10**np.random.uniform(-1,1,size=nSimul)
metabolism = np.random.choice(['common','specific'],size=nSimul)

# run simulations
r_fit = [np.nan for i in range(nSimul)]
p_fit = [np.nan for i in range(nSimul)]
for nsim in range(nSimul):
    
    # check progress
    print(nsim + 1)
    
    ### RUN MODEL
    
    ### MODEL DEFINITION
    
    # general assumptions
    assumptions = community_simulator.usertools.a_default.copy()
    assumptions['n_wells'] = 100
    
    # specify model assumptions based on function input
    assumptions['S'] = nSampledSpeciesPerFamily[nsim] # number of species per family sampled at initialization
    
    assumptions['SA'] = [nSpeciesPerFamily[nsim]] * nFamilies[nsim] # [100, 100, 100] # number of species per specialist family
    assumptions['Sgen'] = 0 # 30 # number of generalists
    assumptions['MA'] = [nResourcesPerClass[nsim]] * nFamilies[nsim]
    assumptions['l'] = leakageFraction[nsim] # leakage fraction
    
    assumptions['response'] = 'type I'
    assumptions['regulation'] = 'energy' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
    assumptions['sampling'] = 'Gamma' # 'Binary', 'Gaussian', 'Uniform' or 'Gamma' (sampling of the matrix c)
    assumptions['supply'] = 'off'
    assumptions['R0_food'] = 1000
    assumptions['m'] = 0 # turn off mortality
    
    assumptions['q'] = preferenceStrength[nsim] #0.9 # preference strength (0 for generalist and 1 for specialist)
    assumptions['c0'] = 0.0 # background consumption rate in binary model
    assumptions['c1'] = 1.0 # specific consumption rate in binary model
    assumptions['sigc'] = 1 #3 # standard deviation of sum of consumption rates for Gaussian and Gamma models
    
    assumptions['sparsity'] = byproductSparsity[nsim] #0.05 # variability in secretion fluxes among resources (must be less than 1)  
    assumptions['fs'] = 0.45 #0.45 # fraction of secretion flux to resources of the same type as the consumed one
    assumptions['fw'] = 0.45 #0.45 # fraction of secretion flux to waste resources
    assumptions['metabolism'] = metabolism[nsim] # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
    
    # parameters
    params = community_simulator.usertools.MakeParams(assumptions)
    
    
    
    ### MAKE FAMILY-STRUCTURED MATRICES
    
    # consumer preference: unstructured case (all species are generalists and share the same metabolic matrix)
    a = assumptions.copy()
    a['q'] = 0
    a['metabolism'] = 'common'
    c = community_simulator.usertools.MakeMatrices(a)[0]
    
    # add family structure: scale uptake rates of non-preferred resources by 1-q
    cf = c.copy()
    for i in range(len(assumptions['SA'])):
        family_i = 'F' + str(i)
        class_i = 'T' + str(i)
        
        cf.loc[family_i] = 1 - assumptions['q']
        cf.loc[family_i,class_i] = 1
    
    c = c*cf
    
    params['c'] = c
    
    # if metabolism is 'specific', make family-specific (rather than species-specific) D matrices
    if assumptions['metabolism'] == 'specific':
        D = []
        for i in range(len(assumptions['MA'])):
            D = D + [community_simulator.usertools.MakeMatrices(a)[1]] * assumptions['SA'][i]
        
        params['D'] = D
    
    
    
    ### MAKE DYNAMICS
    # dynamics
    def dNdt(N,R,params):
        return community_simulator.usertools.MakeConsumerDynamics(assumptions)(N,R,params)
    def dRdt(N,R,params):
        return community_simulator.usertools.MakeResourceDynamics(assumptions)(N,R,params)
    
    dynamics = [dNdt,dRdt]
    
    
    
    ### MAKE INITIAL STATE
    
    # default
    N0, R0 = community_simulator.usertools.MakeInitialState(assumptions)
    
    # use sampleSpecies function (see above)
    for i in range(assumptions['n_wells']):
        N0.iloc[:,i] = sampleSpecies(np.random.randint(1,len(assumptions['SA'])+1),
                                     assumptions['S'],
                                     assumptions,
                                     mode='sequential')
    
    
    
    ### SEED AND STABILIZE COMMUNITIES
    
    # initialize plate
    plate = community_simulator.Community((N0,R0),
                                          dynamics,
                                          params,
                                          parallel=par,
                                          scale=1e6)
    
    # stabilize plates
    N, R = stabilizeCommunities(plate)
    
    
    
    ### GET CORRELATION
    
    # remove extinct wells, proceed only if there is more than one non-extinct well (return nan otherwise)
    not_extinct_wells = np.where(N.sum() > 0)[0]
    
    if not_extinct_wells.size == 0:
        
        r = [np.nan, np.nan]
        
    else:
        
        # keep non-extinct wells only
        N = N.iloc[:,not_extinct_wells]
        
        # get relative abundances
        F = N/N.sum()
    
        # total number of families
        n_families_total = len(assumptions['SA'])
        
        # get number of families and number of species per family for each community
        rel_abundance_threshold = 0
        
        n_families = [np.nan] * assumptions['n_wells']
        n_species_focal = [np.nan] * assumptions['n_wells']
        
        family_indexes = [0] + np.cumsum(assumptions['SA']).tolist()
        
        # run through communities
        for i in range(assumptions['n_wells']):
        
            # run through families
            species_per_family = np.array([0] * n_families_total)
            for j in range(n_families_total):
                F_ij = np.asarray(F.iloc[family_indexes[j]:family_indexes[j+1],i])
                species_per_family[j] = sum(F_ij > rel_abundance_threshold)
                
            # how many families persist in the stabilized community?
            n_families[i] = sum(species_per_family[1:]>0)
            
            # on average, how many species are in each of the persistent families?
            n_species_focal[i] = species_per_family[0]
            
        # get correlation and p-value between number of families and average number of species per family
        r = [np.corrcoef(n_families,n_species_focal)[0,1],
             sm.OLS(n_families,n_species_focal).fit().pvalues[0]]

    # save simulation results (fit to n_species vs. n_species_per_family)
    r_fit[nsim] = r[0]
    p_fit[nsim] = r[1]
    
    # if correlation was above/below a threshold, save data of n_species_focal vs n_families
    if r[0] > 0.3:
        data_out = pd.DataFrame({'n_families':n_families,
                                 'n_species_focal':n_species_focal})
        data_out.to_csv('dbd_rpos.txt',sep='\t',index=False)
    
    if r[0] < -0.6:
        data_out = pd.DataFrame({'n_families':n_families,
                                 'n_species_focal':n_species_focal})
        data_out.to_csv('dbd_rneg.txt',sep='\t',index=False)                                                    
                                                                
# save summary data from all simulations
data_out = pd.DataFrame({'nFamilies':nFamilies,
                         'nSpeciesPerFamily':nSpeciesPerFamily,
                         'nSampledSpeciesPerFamily':nSampledSpeciesPerFamily,
                         'nResourcesPerClass':nResourcesPerClass,
                         'leakageFraction':leakageFraction,
                         'preferenceStrength':preferenceStrength,
                         'byproductSparsity':byproductSparsity,
                         'metabolism':metabolism,
                         'r':r_fit,
                         'p':p_fit})
data_out.to_csv('dbd_corr_vs_params.txt',sep='\t',index=False)








