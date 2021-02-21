from scipy.optimize import linprog
import numpy as np
from math import floor, ceil
import copy


class Node(object):
    def __init__(self, x_bounds=[], freeze_var_list=[], index=0, upper_or_lower=0):
        self._x_bounds = x_bounds
        self._freeze_var_list = freeze_var_list
        self._index = index
        self._upper_or_lower = upper_or_lower

        print("create Node: {}".format(index))
        print('')

    def freeze_var(self, index, val):
        self._x_bounds[index][0] = val
        self._x_bounds[index][1] = val
        self._freeze_var_list.append(index)

    def set_lp_res(self, res):
        self._res = res
        s = ""
        for l in range(len(self._res['x'])):
            if l in self._freeze_var_list:
                s += "[" + str(self._res['x'][l]) + "]"
            else:
                s += " " + str(self._res['x'][l])
        print("x: ", s)

    def check_integer_val_solved(self, m):
        return True if m == len(self._freeze_var_list) else False


class BbAlgorithm(object):
    def __init__(self, c, a_ub, b_ub, x_b, integer_val):
        self.c = c
        self.a_ub = a_ub
        self.b_ub = b_ub
        self.x_b = x_b
        self._integer_val = integer_val
        self.best_solution = float('inf')
        self.best_node = None
        self.nodes = []
        self.nodes_solution = []

    def solve_lp(self, cur_x_b):
        return linprog(self.c, A_ub=self.a_ub, b_ub=self.b_ub, bounds=cur_x_b)

    # def set_lp_res(self, res):
    #     self._res = res
    #     s = ""
    #     for l in range(len(self._res['x'])):
    #         if l in self._freeze_var_list:
    #             s += "[" + str(self._res['x'][l]) + "]"
    #         elif l in self._integer_var:
    #             s += "\'" + str(self._res['x'][l]) + "\' "
    #         else:
    #             s += " " + str(self._res['x'][l])
    #     print("x: ", s)

    def check_fessible(self, res):
        if res['status'] == 0:
            return True
        elif res['status'] == 2:
            return False
        else:
            raise ("Problem Unbounded")

    def add_node(self, node):
        res = self.solve_lp(node._x_bounds)
        if self.check_fessible(res) and res['fun'] < self.best_solution:
            node.set_lp_res(res)
            self.nodes_solution.append(res['fun'])
            self.nodes.append(node)
            if node.check_integer_val_solved(len(self._integer_val)):
                self.best_solution = res['fun']
                self.best_node = node
                print("----------------current solution-------------------")
                print("x: ", node._res['x'])
                print("z: ", node._res['fun'])
                print("---------------------------------------------------\n")
            print("==> Add node to tree list: ", node._index)
            print("==> current nodes: ", self.nodes_solution)
            print("")
            return True
        else:
            print("==> Node infeasible: ", node._index)
            print("==> current nodes: ", self.nodes_solution)
            print("")
            return False

    def del_higher_val_node(self, z_s):
        del_list = []
        for i in range(len(self.nodes_solution)):
            if self.nodes_solution[i] >= z_s:
                del_list.append(i)
        s = ""
        for i in del_list:
            s += " " + str(self.nodes[i]._index)
        print("Remove nodes: ", s)
        self.nodes = list(np.delete(self.nodes, del_list))
        self.nodes_solution = list(np.delete(self.nodes_solution, del_list))
        print("Current nodes: ", self.nodes_solution)
        print("")

    def del_item(self, index):
        print("Remove node: ", self.nodes[index]._index)
        self.nodes = list(np.delete(self.nodes, index))
        self.nodes_solution = list(np.delete(self.nodes_solution, index))
        print("Current nodes: ", self.nodes_solution)
        print("")

    def check_bounds(self, temp_x_b, index, u_or_l):
        if u_or_l == 1:
            if self.x_b[index][0] is not None and temp_x_b[index][0] is None:
                return False
            elif self.x_b[index][0] is None and temp_x_b[index][0] is not None:
                return True
            elif self.x_b[index][0] is not None and temp_x_b[index][0] is not None:
                return False if(self.x_b[index][0] > temp_x_b[index][0]) else True
        elif u_or_l == 2:
            if self.x_b[index][1] is not None and temp_x_b[index][1] is None:
                return False
            elif self.x_b[index][1] is None and temp_x_b[index][1] is not None:
                return True
            elif self.x_b[index][1] is not None and temp_x_b[index][1] is not None:
                return False if(self.x_b[index][1] < temp_x_b[index][1]) else True
        else:
            print("Error of bounds")
            exit()

    def run(self):
        print("####################### start B & B #####################\n")
        node_count = 0
        node = Node(copy.deepcopy(self.x_b), [], node_count)
        node_count += 1
        res = self.solve_lp(self.x_b)

        lower = floor(res['x'][self._integer_val[0]])
        upper = lower + 1

        lower_node = Node(copy.deepcopy(self.x_b), [], node_count, 1)
        lower_node.freeze_var(self._integer_val[0], lower)
        self.add_node(lower_node)
        node_count += 1

        upper_node = Node(copy.deepcopy(self.x_b), [], node_count, 2)
        upper_node.freeze_var(self._integer_val[0], upper)
        self.add_node(upper_node)
        node_count += 1

        while len(self.nodes) > 0:
            index = np.argmin(self.nodes_solution)
            x_b = self.nodes[index]._x_bounds
            freeze_list = self.nodes[index]._freeze_var_list
            res = self.nodes[index]._res
            freeze_var_index = len(freeze_list)

            lower = floor(res['x'][self._integer_val[freeze_var_index]])
            upper = lower + 1
            lower_node = Node(copy.deepcopy(x_b), copy.deepcopy(freeze_list), node_count, 1)
            lower_node.freeze_var(self._integer_val[freeze_var_index], lower)
            self.add_node(lower_node)
            node_count += 1

            upper_node = Node(copy.deepcopy(x_b), copy.deepcopy(freeze_list), node_count, 2)
            upper_node.freeze_var(self._integer_val[freeze_var_index], upper)
            self.add_node(upper_node)
            node_count += 1

            self.del_item(index)
            self.del_higher_val_node(self.best_solution)
            print("########################################")
        print("")
        print("######################### Best Solution #######################")
        print(self.best_node._res)


if __name__ == "__main__":
    integer_val = [0, 1, 3]
    c = [-78, -77, -90, -97, -31]
    a_ub = [
        [11, 4, -41, 44, 7],
        [-87, 33, 24, 14, -13],
        [61, 69, 69, -57, 23]
    ]
    b_ub = [82, 77, 87]
    x_bounds = [[0, None] for _ in range(len(c))]
    bb_algorithm = BbAlgorithm(c, a_ub, b_ub, x_bounds, integer_val)
    bb_algorithm.run()