# Authors: Gabriela Miyazato Szini and Andrea Titton - March 2021

import time
import numpy as np
import multiprocessing as mp

from numba import jit
from itertools import permutations, combinations
from utils.stats import perm_np, s, diff_tetrad, double_diff, s_ij_comb, gen_dgp, combination_fix_ij, s_ij_comb

np.random.seed(110792)
verbose = True

#------------------------- Initialization ------------------------------------#

N = 40 # number of nodes in the network
sim = 10 # number of simulations

#---------------------- Loop for simulations ---------------------------------#

def simulation():

    #-------------------- Data Generating Process-----------------------------#
    n_dyads, n_combinations, x_ij, u_ij = gen_dgp(N)
        
    #-------------------- Obtaining U-statistic ------------------------------#
    tetrads =  perm_np(N, 4)  

    verbose and print("Setting up the U-statistic...")
    
    X = diff_tetrad(x_ij, tetrads)
    U = diff_tetrad(u_ij, tetrads)
    u_stat = np.mean(X * U)
    
    #--------- Calculating the variance of the U-statistic -------------------#
    
    # To calculate the variance we need to obtain \Delta_2.
    # Here we follow a consistent estimator proposed in Graham (2018).
    # To get the \Delta_2 value we need the average of \bar{s}_ij \bar{s}_ij'
    
    verbose and print("Computing Graham estimator...")
    
    s_ij_bar = np.zeros(n_dyads)

    
    for l, dyad in enumerate(permutations(range(N), 2)):
    
        verbose and print(f"Dyad {l+1} / {int(n_dyads)}", end="\r")
    
        i, j = dyad

        # for a given dyad, we need to average the scores of the combinations that contain i and j, so here we take the combinations
        comb = combination_fix_ij(i,j,N)

        # for a given combination, the score consists of an average of its permutations
        s_ij = s_ij_comb(comb, i, j, n_combinations, u_ij, x_ij)

        bar = np.mean(s_ij)
        s_ij_bar[l] = bar
    
    # Delta_2 is simply the average of \bar{s}_ij \bar{s}_ij' over all dyads
    
    delta_2 = np.mean(np.square(s_ij_bar))
    
    verbose and print("\n...done!")
    
    return u_stat, delta_2

def run(_): return simulation()


if __name__ == "__main__":
    
    start = time.time()
    
    with mp.Pool(8) as p:
    
        result = p.map(run, range(sim))

    end = time.time()
    
    print(result)

    with open('results_N40_secondsim_T10.txt', 'w') as fp:
        fp.write('\n'.join('%s %s' % x for x in result))
    
    print(f"Simulation with N={N} and T={sim} took {end-start:.2f} seconds")
