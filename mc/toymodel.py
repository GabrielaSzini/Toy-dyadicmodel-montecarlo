# Authors: Gabriela Miyazato Szini and Andrea Titton - March 2021

import time
import numpy as np
import multiprocessing as mp

from itertools import permutations, combinations

np.random.seed(110792)
verbose = False

#------------------------- Initialization ------------------------------------#

N = 30 # number of nodes in the network
n_dyads = N*(N-1)  # number of directed dyads
sim = 10 # number of simulations


#---------------------------- Utils-------------------------------------------#


def perm_np(to: int, size: int) -> np.array:
    """
    Generates a matrix of "to x size" of permutations up "to"
    """

    it = permutations(range(N), size)
    return np.array(list(it))


def diff_tetrad(A: np.array, t: np.array) -> np.array:
    """
    Differencing out FE of tetrad of A, using tetrads t
    """

    first = A[t[:, 0], t[:, 1]] - A[t[:, 0], t[:, 2]]
    second = A[t[:, 3], t[:, 1]] - A[t[:, 3], t[:, 2]]

    return first - second


def double_diff(mat, i, j, k, l):
    return (mat[i, j] - mat[i, k]) - (mat[l, j] - mat[l, k])


def s(U, X, tetrad):
    """
    Computes s given a tetrad
    """

    return double_diff(X, *tetrad) * double_diff(U, *tetrad)

#---------------------- Loop for simulations ---------------------------------#

def simulation():
    
    
    n_dyads = N * (N - 1)
    n_pairs = n_dyads // 2
    n_combinations = (N-2)*(N-3) // 2
    
    #------------------- Data Generating Process -----------------------------#
    verbose and print("Setting up DGP...")
    
    u_ij = np.random.normal(0, 1, (N, N))  # drawing the error terms
    np.fill_diagonal(u_ij, 0)
    
    # drawing the FE likewise Jochmans (2016)
    a_i = np.random.beta(2, 2, size=(N, 1)) - 1/2
    x_ij = - np.abs(a_i - a_i.T)
    
    
    #------------------- Obtaining U-statistic -------------------------------#
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
    
    s_ij_bar = np.zeros((n_pairs, 2))

    
    for l, dyad in enumerate(combinations(range(N), 2)):
    
        verbose and print(f"Dyad {l+1} / {int(n_dyads/2)}", end="\r")
    
        i, j = dyad
    
        # for a given dyad, we need to average the scores of the combinations that contain i and j, so here we take the combinations
        comb = filter(
            lambda c: (i in c) and (j in c),
            combinations(range(N), 4)
        )
    
        s_ij = np.zeros(n_combinations)  # vector to store those scores
    
        # for a given combination, the score consists of an average of its permutations
        for m, c in enumerate(comb):
    
            s_ij[m] = np.mean([
                s(x_ij, u_ij, perm) for perm in permutations(c)
            ])
        
        bar = np.mean(s_ij)
        s_ij_bar[l, :] = [bar, bar]
    
    # Delta_2 is simply the average of \bar{s}_ij \bar{s}_ij' over all dyads
    s_ij_bar = s_ij_bar.reshape(-1)
    delta_2 = np.mean(np.square(s_ij_bar))
    
    verbose and print("\n...done!")
    
    return u_stat, delta_2

def run(_): return simulation()


if __name__ == "__main__":
    
    start = time.time()
    
    with mp.Pool(4) as p:
    
        result = p.map(run, range(sim))

    end = time.time()
    
    print(result)
    
    print(f"Simulation with N={N} and T={sim} took {end-start:.2f} seconds")
