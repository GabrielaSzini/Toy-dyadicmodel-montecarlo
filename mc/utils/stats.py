import numpy as np

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