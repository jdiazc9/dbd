# -*- coding: utf-8 -*-
"""
Created on Mon Nov 7 18:28:37 2022
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
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette()
# %matplotlib inline

import random
#import math
#import copy
#import statsmodels.api as sm
import time
t0 = time.time()

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

# function to sample species based on experimental distribution
def sampleSpecies_exp(sample_size,assumptions,p):
    
    n_species = sum(assumptions['SA'])
    sp_index = np.arange(0,n_species)
    chosen_sp = np.random.choice(sp_index,size=sample_size,p=p)
    chosen_sp = np.unique(chosen_sp)
    
    N = np.zeros((n_species,1))
    N[chosen_sp,:] = 1
    
    return N

# make consumer preference matrix (slightly different from the built-in MakeMatrices function from the Community Simulator but essentially the same)
def makeCmatrix(assumptions):

    # consumer preference: unstructured case (all species are generalists and share the same metabolic matrix, this will be tuned later)
    M = np.sum(assumptions['MA'])
    T = len(assumptions['MA'])
    S = np.sum(assumptions['SA'])+assumptions['Sgen']
    F = len(assumptions['SA'])
    resource_names = ['R'+str(k) for k in range(M)]
    type_names = ['T'+str(k) for k in range(T)]
    family_names = ['F'+str(k) for k in range(F)]
    consumer_names = ['S'+str(k) for k in range(S)]
    resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])],
                      resource_names]
    consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                      +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
    
    c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)
    
    c_mean = (assumptions['muc']/M)
    c_var = (assumptions['sigc']**2/M)
    thetac = c_var/c_mean
    kc = c_mean**2/c_var
    c.iloc[:][:] = np.random.gamma(kc,scale=thetac,size=(S,M))

    # add family structure: scale uptake rates of non-preferred resources by 1-q
    cf = pd.DataFrame(np.ones((S,M)),columns=resource_index,index=consumer_index) - assumptions['q']
    sp_ranges = np.concatenate(([0],np.cumsum(assumptions['SA'])))
    res_ranges = np.concatenate(([0],np.cumsum(assumptions['MA'])))
    for i in range(len(assumptions['SA'])):   
        cf.iloc[sp_ranges[i]:sp_ranges[i+1],res_ranges[i]:res_ranges[i+1]] = 1
    
    c = c*cf
    
    return c

# stabilize communities through serial passaging
def stabilizeCommunities(plate):
    
    # parameters
    dilution_factor = 1/100
    n_transfers = 10
    time_per_cycle = 0.1
    # extinction_threshold = 1e-4 # we don't use this because we don't use StadyState()
    
    # run stabilization protocol
    for i in range(n_transfers-1):
        print('Transfer '+str(i))
        plate.Propagate(time_per_cycle,compress_species=True,compress_resources=False)
        plate.N = plate.N * dilution_factor
        plate.R = plate.R * dilution_factor
        plate.R.iloc[0,:] = plate.R0.iloc[0,:]
    plate.Propagate(time_per_cycle,compress_species=True,compress_resources=False)
    
    return([plate.N,plate.R])

# stabilize communities through serial passaging (explicit, doesn't parallelize but shows steps, good if we need to track progress)
def stabilizeCommunities_explicit(plate):
    
    # parameters
    dilution_factor = 1/100
    n_transfers = 10
    time_per_cycle = 0.1
    # extinction_threshold = 1e-4 # we don't use this because we don't use StadyState()
    
    wells = list(plate.N.columns)
    
    # run stabilization protocol
    for well in wells:
        print('Well '+well)
        for i in range(n_transfers):
            print('   Transfer '+str(i))
            ttraj,Ntraj,Rtraj = plate.TestWell(T=time_per_cycle,compress_resources=False,compress_species=True,well_name=well,ns=100,show_plots=False)
            plate.N[well] = Ntraj[99,:].T * dilution_factor # dilute consumers
            plate.R[well] = Rtraj[99,:].T * dilution_factor # dilute resources
            plate.R[well][0] = plate.R0[well][0] # replenish external resource
        # final transfer
        ttraj,Ntraj,Rtraj = plate.TestWell(T=time_per_cycle,compress_resources=False,compress_species=True,well_name=well,ns=100,show_plots=False)
        plate.N[well] = Ntraj[99,:].T
        plate.R[well] = Rtraj[99,:].T
            
    return([plate.N,plate.R])

# heatmaps of consumer preference and metabolic matrices
def plotMatrices(c,D):
    
    fig,ax=plt.subplots(1,2)
    
    sns.heatmap(c,cmap='Greys',vmin=0,cbar=False,ax=ax[0])
    ax[0].set_title('Consumer preference\nmatrix')
    ax[0].set_xlabel('Resources')
    ax[0].set_ylabel('Consumers')
    
    sns.heatmap(D,cmap='Greys',vmin=0,cbar=False,ax=ax[1])
    ax[1].set_title('Metabolic\nmatrix')
    ax[1].set_xlabel('Consumed resources')
    ax[1].set_ylabel('Secreted resources')
    
    fig.set_size_inches(15,7)
    fig.tight_layout(pad=3.0)
    
    return fig, ax

# stack plot (for community compositions)
def myStackPlot(N,max_wells=10):

    # check number of wells and species
    n_species = N.shape[0]
    
    # if there are more wells than the maximum to plot, trim N
    if N.shape[1] > max_wells:
        N = N.iloc[:,:max_wells]
    n_wells = N.shape[1]
    
    # normalize each column to get fractions
    N = N/N.sum()
    
    # numpy data frame to list of lists (list of as many lists as species in total)
    N = N.values.tolist()
    
    # colormap (adjusted to number of species)
    random.seed(0)
    cmap = plt.get_cmap("gist_stern")
    cmap = cmap(np.linspace(0,1,num=n_species)).tolist()
    random.shuffle(cmap)
    
    # position in the x axis
    xpos = np.arange(n_wells).tolist()
    xticklabs = np.add(xpos,[1]).tolist()
    
    # define cumulative sum of abundances
    cumsum = np.zeros(n_wells).tolist()
    
    # initialize figure
    fig, ax = plt.subplots(1)
    
    # plot bars - iterate across species and stack
    for i in range(n_species):
        ax.bar(xpos,N[i],color=cmap[i],bottom=cumsum)
        cumsum = np.add(N[i],cumsum).tolist()
        
    ax.set_xticks(xpos)
    ax.set_xticklabels(xticklabs)
    ax.set_ylabel('Fraction')
    ax.set_ylim(0,1)
        
    return fig, ax

# rank plot
def myRankPlot(N):
    
    # abundance threshold
    u = 1/10000
    
    # check number of wells
    n_wells = N.shape[1]
    
    # normalize to get fractions
    N = N/N.sum()
        
    # initialize figure
    fig, ax = plt.subplots(1)
    
    # plot
    for i in range(n_wells):
        
        # get ith column of N, remove elements below threshold and sort it
        Ni = N.iloc[:,i].to_numpy()
        Ni = Ni[Ni>u]
        Ni.sort()
        Ni = Ni[::-1]
        
        # plot
        ax.plot(Ni)
    
    ax.set_yscale("log")
    ax.set_ylabel('Relative abundance')
    ax.set_xlabel("Rank")
    ax.set_ylim(u,1e0)
    ax.set_aspect(5)
        
    return fig, ax

# get correlation between number of focal species and number of non-focal families
def getDBDcor(N,focal_families):
    
    abundance_threshold = 1/10000
    
    # get number of species in focal family
    Nr = N/N.sum()
    n_fSpecies = [np.sum(Nr[Nr.columns[i]][focal_families[i]] > abundance_threshold) for i in range(len(Nr.columns))]
    
    # get number of non-focal families
    fam_abundances = N.groupby(level=0).sum()
    fam_rabundances = fam_abundances/fam_abundances.sum()
    n_nfFamilies = [np.sum(fam_rabundances[fam_rabundances.columns[i]] > abundance_threshold) - 1 for i in range(len(fam_rabundances.columns))]
    
    # get correlation
    corr = np.corrcoef(n_nfFamilies,n_fSpecies)[0,1]
    
    return corr



# %%

# number of simulations to run
nSimul = 100

# initialize parameter arrays
preferenceStrength = np.random.uniform(0.0,1.0,size=nSimul)
byproductSparsity = 10**np.random.uniform(-1,1,size=nSimul)
metabolism = np.random.choice(['common','specific'],size=nSimul)

# distribution of species abundances in empirical data
distrib_exp = pd.read_csv('../data/inoc.csv')
sp_per_fam_exp = pd.read_csv('../data/sp_per_fam.csv')

# run simulations
r_inoc = [np.nan for i in range(nSimul)]
r_stable = [np.nan for i in range(nSimul)]
for nsim in range(nSimul):
    
    # check progress
    print('--- RUN '+str(nsim + 1)+' ---')
    
    ### RUN MODEL
    
    ### MODEL DEFINITION
    
    # general assumptions
    assumptions = community_simulator.usertools.a_default.copy()
    assumptions['n_wells'] = 100
    
    # specify model assumptions based on function input
    assumptions['S'] = 100 # number of species per family sampled at initialization
    
    assumptions['SA'] = list(sp_per_fam_exp['nSpecies']) # [100, 100, 100] # number of species per specialist family
    assumptions['Sgen'] = 0 # 30 # number of generalists
    assumptions['MA'] = [5] * len(assumptions['SA'])
    assumptions['l'] = 0.8 # leakage fraction
    
    assumptions['response'] = 'type I'
    assumptions['regulation'] = 'energy' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
    assumptions['sampling'] = 'Gamma' # 'Binary', 'Gaussian', 'Uniform' or 'Gamma' (sampling of the matrix c)
    assumptions['supply'] = 'off'
    assumptions['R0_food'] = 10000
    assumptions['m'] = 0 # turn off mortality
    
    assumptions['q'] = preferenceStrength[nsim] #0.9 # preference strength (0 for generalist and 1 for specialist)
    assumptions['c0'] = 0.0 # background consumption rate in binary model
    assumptions['c1'] = 1.0 # specific consumption rate in binary model
    assumptions['muc'] = 10 # mean of sum of consumption rates for Gaussian and Gamma models
    assumptions['sigc'] = 1 #3 # standard deviation of sum of consumption rates for Gaussian and Gamma models
    
    assumptions['sparsity'] = byproductSparsity[nsim] #0.05 # variability in secretion fluxes among resources (must be less than 1)  
    assumptions['fs'] = 0.05 #0.45 # fraction of secretion flux to resources of the same type as the consumed one
    assumptions['fw'] = 0.05 #0.45 # fraction of secretion flux to waste resources
    assumptions['metabolism'] = metabolism[nsim] # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
    
    # parameters
    params = community_simulator.usertools.MakeParams(assumptions)
    
    
    
    ### MAKE FAMILY-STRUCTURED MATRICES
    
    # family-structured consumer matrix
    params['c'] = makeCmatrix(assumptions)
    
    # if metabolism is 'specific', make family-specific D matrices
    if assumptions['metabolism'] == 'specific':
        a = assumptions.copy()
        a['metabolism'] = 'common'
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
    
    # use sampleSpecies function for initial species abundances (see above)
    for i in range(assumptions['n_wells']):
        N0.iloc[:,i] = sampleSpecies_exp(np.random.randint(20,200 + 1,size=1), # sample 10 to 1000 individuals
                                         assumptions, p=distrib_exp['Relative_Abundance.avg'])
    
    
    ### SEED AND STABILIZE COMMUNITIES
    
    # initialize plate
    plate = community_simulator.Community((N0,R0),
                                          dynamics,
                                          params,
                                          parallel=par,
                                          scale=1e6)
    
    # save initial state
    plate.N.to_csv('../simul_runs/Run'+str(nsim)+'_T0.txt',sep='\t',index=True)
    
    # stabilize plates
    N, R = stabilizeCommunities_explicit(plate)
    
    # save final communities
    plate.N.to_csv('../simul_runs/Run'+str(nsim)+'_T20.txt',sep='\t',index=True)
    
    
    
    ### PRINT ELAPSED TIME
    print(time.time() - t0)



### SAVE DATA
                                                                  
# save summary data from all simulations
data_out = pd.DataFrame({'preferenceStrength':preferenceStrength,
                         'byproductSparsity':byproductSparsity,
                         'metabolism':metabolism})
data_out.to_csv('../simul_params/dbd_corr_vs_params.txt',sep='\t',index=False)
    
    
    
    
    
    
    
    
