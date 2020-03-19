#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 18:22:32 2020

@author: youssouf
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm

##################################################################################################
## Keys
##################################################################################################
state = 'state'
previous_state = 'previous_state'
resume_time = 'resume_time'
exposition_time = 'exposition_time'
exposed_at = 'exposed_at'
infected_at = 'infected_at'
##################################################################################################
## States :
##################################################################################################
susceptible = 'S'
infected = 'I'
recovered = 'R'
dead = 'D'
exposed = 'E'
quaranted = 'Q'
###################################################################################################
## Population parameters :
###################################################################################################
population_size = 4000000
exposition_size = 100
###################################################################################################
## Growth factor :
###################################################################################################
growth_factor = 2.5
###################################################################################################
## Graph of population :
###################################################################################################
graph = nx.fast_gnp_random_graph(population_size, exposition_size/population_size)
###################################################################################################
## Time parameters :
###################################################################################################
periods = 210
mean_recovery_time = 15
mean_incubation_time = 10
mean_exposition_time = 12
gsc = 1.5
gsh = int(mean_recovery_time/gsc)
##################################################################################################
## Rates :
##################################################################################################
exposition_rate = growth_factor / exposition_size
death_rate = 0.01
##################################################################################################
## Start time of the infection:
##################################################################################################
start_time = 0
##################################################################################################
## Number of simulations :
##################################################################################################
number_of_simulations = 1000
##################################################################################################
## Compute the evolution of cases
##################################################################################################
def evolution(infection_rate):
    total_dcases = np.zeros(periods)

    for _ in tqdm(range(number_of_simulations)):

        for i in list(graph.nodes):
            graph.nodes[i][state] = susceptible

        # number of infected daily cases
        dcases = np.zeros(periods)

        for t in tqdm(range(periods)):
            if t >= start_time and t < start_time + 1:
                graph.nodes[0][state] = infected
                graph.nodes[0][infected_at] = t
                graph.nodes[0][resume_time] = np.random.gamma(shape=gsh, scale=gsc)

            ### get the previous state
            for i in list(graph.nodes):
                graph.nodes[i][previous_state] = graph.nodes[i][state]

            ### Change the state
            for i in list(graph.nodes):
                # if the node was exposed
                if graph.nodes[i][previous_state] == exposed:
                    # if the exposition time has been passed
                    if t - graph.nodes[i][exposed_at] >= graph.nodes[i][exposition_time]:
                        u = np.random.rand()
                        # the node will be infected
                        if u <= infection_rate:
                            graph.nodes[i][state] = infected
                            dcases[t] += 1
                            graph.nodes[i][infected_at] = t
                            graph.nodes[i][resume_time] = np.random.gamma(shape=gsh, scale=gsc)
                        # the node will be susceptible again
                        else:
                            graph.nodes[i][state] = susceptible
                # if the node is infected
                if graph.nodes[i][previous_state] == infected:
                    # the node has been resumed
                    if t - graph.nodes[i][infected_at] >= graph.nodes[i][resume_time]:
                        u = np.random.rand()
                        # the node is dead
                        if u <= death_rate:
                            graph.nodes[i][state] = dead
                        # the node has been recovered
                        else:
                            graph.nodes[i][state] = recovered
                if graph.nodes[i][previous_state] == susceptible:
                    sp = 1.
                    for j in list(graph.adj[i]):
                        if graph.nodes[j][previous_state] == infected:
                            sp *= exposition_rate
                    u = np.random.rand()
                    if u <= sp:
                        graph.nodes[i][state] = susceptible
                    else:
                        graph.nodes[i][state] = exposed
                        graph.nodes[i][exposed_at] = t
                        graph.nodes[i][exposition_time] = np.random.exponential(scale=mean_exposition_time)
                        dcases[t] += 1
        total_dcases += dcases

    mean_dcases = total_dcases * 1./number_of_simulations
    return mean_dcases
##################################################################################################

##################################################################################################
## Start with low infection rate due to the extreme quarantine measures
##################################################################################################
infection_rate = 0.2
##################################################################################################
##################################################################################################
##################################################################################################
extreme_quarantine_dcases = evolution(infection_rate)
##################################################################################################

##################################################################################################
## Second with medium infection rate due to the medium quarantine measures
##################################################################################################
infection_rate = 0.5
##################################################################################################
##################################################################################################
##################################################################################################
medium_quarantine_dcases = evolution(infection_rate)
##################################################################################################

##################################################################################################
## Finally with high infection rate due to the non quarantine measures
##################################################################################################
infection_rate = 0.95
##################################################################################################
##################################################################################################
##################################################################################################
low_quarantine_dcases = evolution(infection_rate)
##################################################################################################

##################################################################################################
## Plot the results
##################################################################################################
plt.figure(figsize=(16, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)


times = range(0, periods)

lot_size = 40000
# yticks = 0.01 * np.arange(0, lot_size+step, step)
# number_of_lots = int(population_size/lot_size)
# plt.ylim(0, 0.01*lot_size)
plt.xlim(0, periods)
# plt.yticks(yticks, [str(y) for y in yticks], fontsize=14)
plt.xticks(fontsize=14)
#for y in yticks:
#    plt.plot(times, [0.01*y] * len(times), "--", lw=0.5, color="black", alpha=0.3)
for dcases in [low_quarantine_dcases, medium_quarantine_dcases, extreme_quarantine_dcases]:
    plt.plot(times,
            dcases*lot_size*1./population_size,
            lw=2.5)

# plt.text(periods/2, lot_size + step, "Number of daily cases per " + str(lot_size) , fontsize=17, ha="center")

plt.savefig("simulation_dcases.png", bbox_inches="tight")
