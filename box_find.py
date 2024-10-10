#!/usr/bin/env python3

from Compiler.types import sint, Array, Matrix
from Compiler.library import for_range_opt
from Compiler.util import min as secmin

def boxFindRetAll(D, E, ways, i, j, p, t):
    '''Find possible costs to reach a cell and return all'''
    values = Array(len(ways)-1, sint)
    values.assign_all(len(E))
    ct = 0
    for pi in range(1, len(ways)):
        ct+=1
        w = ways[pi]
        # Implementation note: Many 4090 values added for uniform array length (for multiproc)
        # Can skip all those.
        if w[0] == 4090:
            break
        # If destination exceeds diagonal boundary, skip.
        if j+w[2] < -p+i+w[1] or j+w[2] >= -p+i+t+w[1]:
            continue
        value = w[0]
        value += D[i+w[1]+1][j+w[2]+1]
        for pj in range(3, len(w), 2):
            # Implementation note: Many -1 values added for uniform array length (for multiproc)
            # Can skip all those.
            if w[pj] < 0:
                break
            vi = i+w[pj]
            vj = j+w[pj+1]
            # If path exceeds diagonal boundary, skip.
            if vj < -p+vi or vj >= -p+t+vi:
                value = len(E)
                break
            vj = p+(vj-vi)
            value += E[vi][vj]
        values[pi-1] = value
    return values[:ct]

def boxFindRetMin(D, E, ways, i, j, p):
    '''
    Find possible costs to reach a cell and return min.
    @see boxFindRetAll for notes
    '''
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
            vj = p+(vj-vi)
            value += E[vi][vj]
        minval = minval.min(value)
    return minval

def boxFindRetAll_with_ab(D, a, b, ways, i, j):
    '''
    Find possible costs to reach a cell and return all.
    @see boxFindRetAll for notes
    '''
    values = Array(len(ways)-1, sint)
    values.assign_all(len(a))
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
            value += a[vi].not_equal(b[vj])
        values[pi-1] = value
    return values

def boxFindMultiProc(D, E, ways, i, j, p):
    '''Find possible costs to reach a cell and return min using multiprocessing'''
    values = Array(len(ways)-1, sint)
    values.assign_all(len(E))
    @for_range_opt(1, len(ways), budget=10)
    def _(pi):
        w = ways[pi]
        value = w[0]
        value += D[i+w[1]+1][j+w[2]+1]
        for pj in range(3, len(w), 2):
            vi = (i+w[pj]).max(0)
            vj = (j+w[pj+1]).max(0)
            vj = p+(vj-vi)
            value += E[vi][vj]
        values[pi-1] = value
    return secmin(values)