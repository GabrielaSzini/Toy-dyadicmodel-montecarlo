import numpy as np
import numba as numba
from numba import jit
from itertools import permutations, combinations


def perm_np(N: int, size: int) -> np.array:
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

@jit(nopython=True)
def s(U, X, tetrad):
    """
    Computes s given a tetrad
    """
    return double_diff(X, *tetrad) * double_diff(U, *tetrad)

def s_ij_comb(i, j, k, l, U, X):
    """
    Computes the summand in s_ij_bar
    """
    sij_comb = 1/24 * (
        (double_diff(X, i, j, k, l)) * U[i,j] +
        (double_diff(X, i, k, j, l)) * (- U[i,j]) +
        (double_diff(X, k, j, l, i)) * (- U[i,j]) +
        (double_diff(X, l, k, j, i)) * U[i,j] +
        (double_diff(X, k, l, j, i)) * U[i,j] +
        (double_diff(X, l, j, k, i)) * (- U[i,j]) +
        (double_diff(X, i, j, l, k)) * U[i,j] +
        (double_diff(X, i, l, j, k)) * (- U[i,j]) 
    )
    return sij_comb

@jit(nopython=True)
def gen_dgp(N):
    """
    Generates the data generating process, with number of dyads, number of 
    combinations, the u-statistic, the transformed X and U
    """
    n_dyads = N * (N - 1)
    n_combinations = (N-2)*(N-3) // 2
    # Data Generating Process 
    # verbose and print("Setting up DGP...")
    u_ij = np.random.normal(0, 1, (N, N))  # drawing the error terms
    np.fill_diagonal(u_ij, 0) 
    # drawing the FE likewise Jochmans (2016)
    a_i = np.random.beta(2, 2, size=(N, 1)) - 1/2
    x_ij = - np.abs(a_i - a_i.T)
    return n_dyads, n_combinations, x_ij, u_ij

def combination_fix_ij(i,j,N):
    """
    Generates the combinations of 2 indices given i, and j fixed
    """
    comb = combinations((n for n in range(N) if (n != i) and (n != j)),
        2
    )
    return np.array(list(comb))

def s_ij(comb, i, j, n_combinations):      #can use @jit if I manage to use it for double_diff and s_ij_comb
    """
    Generates the scores for a fixed i, j, and all possible combinations of k, l 
    """
    s_ij = np.zeros(n_combinations) 
    for c in np.arange(0,n_combinations):
        s_ij[c] = s_ij_comb(i, j, comb[c][0], comb[c][1], u_ij, x_ij)
    return s_ij
        
    
