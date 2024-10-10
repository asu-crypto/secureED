#!/usr/bin/env python3

import sys
sys.path.insert(0, '../secure-edit-distance')

from preprocess import load_D_box
from ways import readWays, readWaysMultiProcMin

from Compiler.types import sint, sintbit, Array, Matrix
from Compiler.library import for_range_opt, print_ln, for_range_parallel
from Compiler.util import min as secmin

'''
This is all evaluation code for filling all boxes without Ukkonen. 
Alternatively, Vanegas et al.'s code can be used.
'''

def all_equality(a, b):
    E = Matrix(len(a), len(b), sintbit)
    for i in range(0, len(a)):
    # @for_range_opt(0, len(a))
    # def _(i):
        for j in range(0, len(b)):
        # @for_range_opt(0, len(b))
        # def _(j):
            E[i][j] = a[i].not_equal(b[j])
    # print_ln("%s", E.reveal())
    return E

def boxFind(D, E, ways, i, j):
    '''Find possible costs to reach a cell'''
    minval = sint(len(E))
    for pi in range(1, len(ways)):
        w = ways[pi]
        if w[0] == 4090:
            break
        value = w[0]
        value += D[i+w[1]+1][j+w[2]+1]
        for pj in range(3, len(w), 2):
            if w[pj] < 0:
                break
            vi = i+w[pj]
            vj = j+w[pj+1]
            value += E[vi][vj]
        minval = minval.min(value)
    return minval

def boxFindRetAll(D, E, ways, i, j):
    '''Find possible costs to reach a cell'''
    values = Array(len(ways)-1, sint)
    values.assign_all(len(E))
    for pi in range(1, len(ways)):
        w = ways[pi]
        if w[0] == 4090:
            break
        value = w[0]
        value += D[i+w[1]+1][j+w[2]+1]
        # print_ln("%s", D[i+w[1]+1][j+w[2]+1].reveal())
        # print("Adding", i+w[1]+1, j+w[2]+1)
        for pj in range(3, len(w), 2):
            if w[pj] < 0:
                break
            vi = i+w[pj]
            vj = j+w[pj+1]
            value += E[vi][vj]
        values[pi-1] = value
    # print("gg")
    return values

def fillBox(D, E, i, j, ww):
    for x in range(len(ww)):
        ix = ww[x][0][0]
        iy = ww[x][0][1]
        D[i+ix+1][j+iy+1] = boxFind(D, E, ww[x], i, j)

def boxFill(a, b, tau, m, n):

    # global g1, g2, g3
    E = all_equality(a, b)
    D = load_D_box(E, m, n, n+m-1, m-1)
    waysInBox = readWays(tau)
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            fillBox(D, E, i, j, waysInBox)
    # print("Tau, GG", tau, g1, g2, g3)
    return D[m][n]

def boxFillMultiProcMin(a, b, tau, m, n):

    E = all_equality(a, b)
    D = load_D_box(E, m, n, n+m-1, m-1)
    _, waysInBox = readWaysMultiProcMin(tau)
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            mvees = Matrix(len(waysInBox), len(waysInBox[0])-1, sint)
            for x in range(len(waysInBox)):
                mvees[x] = boxFindRetAll(D, E, waysInBox[x], i, j)
            @for_range_parallel(len(waysInBox), len(waysInBox))
            def _(x):
                ix = waysInBox[x][0][0]
                iy = waysInBox[x][0][1]
                D[i+ix+1][j+iy+1] = secmin(mvees[x])
    # print_ln("D %s", D.reveal())
    return D[m][n]