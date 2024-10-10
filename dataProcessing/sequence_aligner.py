#!/usr/bin/env python3

import sys
sys.setrecursionlimit(5000)

def sequence_aligner(a, b, D):

    a_aligned = ''
    b_aligned = ''

    def backtrack(i, j, ct):
        nonlocal a_aligned, b_aligned
        if i == 0 and b == 0:
            return
        if i == 0:
            a_aligned += '_'*j
            b_aligned += b[:j][::-1]
            return
        if j == 0:
            a_aligned += a[:i][::-1]
            b_aligned += '_'*i
            return
        if D[i-1][j-1] <= D[i-1][j] and D[i-1][j-1] <= D[i][j-1]:
            a_aligned += a[i-1]
            b_aligned += b[j-1]
            if D[i-1][j-1] != D[i][j]:
                ct+=1
            backtrack(i-1, j-1, ct)
        elif D[i-1][j] < D[i][j-1]:
            a_aligned += a[i-1]
            b_aligned += '_'
            backtrack(i-1, j, ct+1)
        else:
            a_aligned += '_'
            b_aligned += b[j-1]
            backtrack(i, j-1, ct+1)

    backtrack(len(a), len(b), 0)
    return a_aligned[::-1], b_aligned[::-1]