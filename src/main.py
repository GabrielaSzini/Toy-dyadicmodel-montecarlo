from itertools import permutations, combinations, product
from collections import Counter

# checking number of combinations/permutations

N = [1, 2, 3, 4, 5]
positive = set([0, 3])
negative = set([1, 2])

def variance_order(i, j, k, l):
    return [(i, j), (i, k), (l, j), (l, k)]

def compare_orders(o1, o2):
    same = []
    for i, tup1 in enumerate(o1):
        for j, tup2 in enumerate(o2):
            if tup1 == tup2:
                f, s = tup1
                sign = (i in positive and j in positive) or (i in negative and j in negative)
                same.append((f"{f}{s}", sign))
    return same

def N_var_order(N, m):
    orders = list(range(1, 5 + 1))
    for comb in combinations(orders, m):
        for i, j, k, l in permutations(comb):
            yield variance_order(i, j, k, l)


x1u1 = variance_order(1,2,3,4)
lhs = [f"{st[0]}{st[1]}" for st in x1u1]
common = [compare_orders(x1u1, other) for other in N_var_order(5, 4)]

    
with open("compinations_results.txt", "w+") as file:
    for i, row1 in enumerate(N_var_order(5,4)):
        row2 = common[i]

        this_row = str(row1) + str(row2) + "\n"
        file.write(this_row)


x2u2 = variance_order(5,2,3,4)
lhs2 = [f"{st[0]}{st[1]}" for st in x2u2]
common2 = [compare_orders(x2u2, other) for other in N_var_order(5, 4)]

    
with open("compinations_results_2.txt", "w+") as file:
    for i, row1 in enumerate(N_var_order(5,4)):
        row2 = common2[i]

        this_row = str(row1) + str(row2) + "\n"
        file.write(this_row)