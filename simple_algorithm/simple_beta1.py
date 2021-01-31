import numpy as np

class Simplex(object):
    def __init__(self, obj, max_mode=False):
        self.max_mode = max_mode  # default is solve min LP, if want to solve max lp,should * -1
        self.mat = np.array([[0] + obj]) * (-1 if max_mode else 1)

    def add_constraint(self, a, b):
        self.mat = np.vstack([self.mat, [b] + a])

    def solve(self):
        m, n = self.mat.shape  # m - 1 is the number slack variables we should add
        temp, B = np.vstack([np.zeros((1, m - 1)), np.eye(m - 1)]), list(range(n - 1, n + m - 1))  # add diagonal array
        mat = self.mat = np.hstack([self.mat, temp])  # combine them!
        while mat[0, 1:].min() < 0:
            col = np.where(mat[0, 1:] < 0)[0][0] + 1  # use Bland's method to avoid degeneracy. use mat[0].argmin() ok?
            row = np.array([mat[i][0] / mat[i][col] if mat[i][col] > 0 else 0x7fffffff for i in
                            range(1, mat.shape[0])]).argmin() + 1  # find the theta index
            if mat[row][col] <= 0: return None  # the theta is ∞, the problem is unbounded
            mat[row] /= mat[row][col]
            ids = np.arange(mat.shape[0]) != row
            mat[ids] -= mat[row] * mat[ids, col:col + 1]  # for each i!= row do: mat[i]= mat[i] - mat[row] * mat[i][col]
            B[row] = col
        return mat[0][0] * (1 if self.max_mode else -1), {B[i]: mat[i, 0] for i in range(1, m) if B[i] < n}

"""
       minimize -x1 - 14x2 - 6x3
       st
        x1 + x2 + x3 <=4
        x1 <= 2aaaa
        x3 <= 3
        3x2 + x3 <= 6
        x1 ,x2 ,x3 >= 0
       answer :-32
    """
t = Simplex([-1, -14, -6])
t.add_constraint([1, 1, 1], 4)
t.add_constraint([1, 0, 0], 2)
t.add_constraint([0, 0, 1], 3)
t.add_constraint([0, 3, 1], 6)
print(t.solve())
print(t.mat)
