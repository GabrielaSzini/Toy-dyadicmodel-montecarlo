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