#!/usr/bin/env python3

import sys
sys.path.insert(0, '../secure-edit-distance')

from Compiler.util import min as secmin
from Compiler.library import for_range

from preprocess import load_D_box, iterateTauedDiagonal

'''
Easier alternative to avoid any mistakes/inefficiencies related to filling boxes.
Paths for tau = 3, 4 are hardcoded here. Not scalable to different taus.
To scale for different taus, use ukkonen_fill
'''

def ukkonenFill4Tau(E, t, m, n, tt=0):

    tau = 4
    p = (t-(n-m)) // 2 # Position of leading diagonal
    kval = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)

    howmany = (p//tau)*tau
    p = pp
    thiss = n-tau

    @for_range(0, m-tau, tau)
    def _(i):
        start = (kval+i > 0).if_else(i-howmany, 0)
        end = (kval+i+t < thiss).if_else(kval+i+t, thiss)
        @for_range(start, end, tau)
        def _(j):

            terms_i1_j4 = [
                4 + D[i+2][j+1],
                1 + D[i+1][j+5],
                3 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)],
                0 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i2_j4 = [
                4 + D[i+3][j+1],
                2 + D[i+1][j+5],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                0 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i3_j4 = [
                4 + D[i+4][j+1],
                3 + D[i+1][j+5],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+4)-(i+3)]
            ]
            terms_i3_j4_part2 = [
                3 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                0 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+4] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i4_j1 = [
                1 + D[i+5][j+1],
                4 + D[i+1][j+2],
                3 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                0 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)]
            ]

            terms_i4_j2 = [
                2 + D[i+5][j+1],
                4 + D[i+1][j+3],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                0 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                3 + D[i+1][j+2] + E[i+4][p+(j+2)-(i+4)],
                3 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)]
            ]

            terms_i4_j3 = [
                3 + D[i+5][j+1],
                4 + D[i+1][j+4],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                0 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+4][j+1] + E[i+4][p+(j+3)-(i+4)]
            ]
            terms_i4_j3_part2 = [
                2 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                2 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+1][j+3] + E[i+4][p+(j+3)-(i+4)],
                3 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)],
                3 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            terms_i4_j4 = [
                4 + D[i+5][j+1],
                4 + D[i+1][j+5],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                0 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)]
            ]
            terms_i4_j4_part2 = [
                3 + D[i+4][j+1] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+3)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                3 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+1][j+4] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+1][j+4] + E[i+3][p+(j+4)-(i+3)],
                3 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]
            
            D[i+1+1][j+4+1] = secmin(terms_i1_j4)
            D[i+2+1][j+4+1] = secmin(terms_i2_j4)
            # Spliting large arrays into two makes secmin efficient (rounds reduction)
            D[i+3+1][j+4+1] = secmin(terms_i3_j4).min(secmin(terms_i3_j4_part2))
            D[i+4+1][j+1+1] = secmin(terms_i4_j1)
            D[i+4+1][j+2+1] = secmin(terms_i4_j2)
            D[i+4+1][j+3+1] = secmin(terms_i4_j3).min(secmin(terms_i4_j3_part2))
            D[i+4+1][j+4+1] = secmin(terms_i4_j4).min(secmin(terms_i4_j4_part2))

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())

    return D[m][n]

def ukkonenFill4Tau_opt_version(E, t, m, n, tt=0):

    tau = 4
    p = (t-(n-m)) // 2 # Position of leading diagonal
    kval = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    p = pp

    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            if j < kval or j >= kval+t:
                continue

            terms_i1_j4 = [
                4 + D[i+2][j+1],
                1 + D[i+1][j+5],
                3 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)],
                0 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i2_j4 = [
                4 + D[i+3][j+1],
                2 + D[i+1][j+5],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                0 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i3_j4 = [
                4 + D[i+4][j+1],
                3 + D[i+1][j+5],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+4)-(i+3)]
            ]
            terms_i3_j4_part2 = [
                3 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                0 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+4] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i4_j1 = [
                1 + D[i+5][j+1],
                4 + D[i+1][j+2],
                3 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                0 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)]
            ]

            terms_i4_j2 = [
                2 + D[i+5][j+1],
                4 + D[i+1][j+3],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                0 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                3 + D[i+1][j+2] + E[i+4][p+(j+2)-(i+4)],
                3 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)]
            ]

            terms_i4_j3 = [
                3 + D[i+5][j+1],
                4 + D[i+1][j+4],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                0 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+4][j+1] + E[i+4][p+(j+3)-(i+4)]
            ]
            terms_i4_j3_part2 = [
                2 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                2 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+1][j+3] + E[i+4][p+(j+3)-(i+4)],
                3 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)],
                3 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            terms_i4_j4 = [
                4 + D[i+5][j+1],
                4 + D[i+1][j+5],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                0 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)]
            ]
            terms_i4_j4_part2 = [
                3 + D[i+4][j+1] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+3)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                3 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+1][j+4] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+1][j+4] + E[i+3][p+(j+4)-(i+3)],
                3 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]
            
            D[i+1+1][j+4+1] = secmin(terms_i1_j4)
            D[i+2+1][j+4+1] = secmin(terms_i2_j4)
            # Spliting large arrays into two makes secmin efficient (rounds reduction)
            D[i+3+1][j+4+1] = secmin(terms_i3_j4).min(secmin(terms_i3_j4_part2))
            D[i+4+1][j+1+1] = secmin(terms_i4_j1)
            D[i+4+1][j+2+1] = secmin(terms_i4_j2)
            D[i+4+1][j+3+1] = secmin(terms_i4_j3).min(secmin(terms_i4_j3_part2))
            D[i+4+1][j+4+1] = secmin(terms_i4_j4).min(secmin(terms_i4_j4_part2))
        kval+=tau

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())

    return D[m][n]

def ukkonenFill3Tau(E, t, m, n, tt=0):

    tau = 3
    p = (t-(n-m)) // 2 # Position of leading diagonal
    kval = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    howmany = (p//tau)*tau
    p = pp
    thiss = n-tau

    @for_range(0, m-tau, tau)
    def _(i):
        start = (kval+i > 0).if_else(i-howmany, 0)
        end = (kval+i+t < thiss).if_else(kval+i+t, thiss)
        @for_range(start, end, tau)
        def _(j):

            term_i1_j3 = [
                3 + D[i+2][j+1],
                1 + D[i+1][j+4],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)],
                0 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            term_i2_j3 = [
                3 + D[i+3][j+1],
                2 + D[i+1][j+4],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                0 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            term_i3_j1 = [
                1 + D[i+4][j+1],
                3 + D[i+1][j+2],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                0 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)]
            ]

            term_i3_j2 = [
                2 + D[i+4][j+1],
                3 + D[i+1][j+3],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                0 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                2 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)]
            ]

            term_i3_j3 = [
                3 + D[i+4][j+1],
                3 + D[i+1][j+4],
                0 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                1 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]
            
            D[i+1+1][j+3+1] = secmin(term_i1_j3)
            D[i+2+1][j+3+1] = secmin(term_i2_j3)
            D[i+3+1][j+1+1] = secmin(term_i3_j1)
            D[i+3+1][j+2+1] = secmin(term_i3_j2)
            D[i+3+1][j+3+1] = secmin(term_i3_j3)

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())

    return D[m][n]

def ukkonenFill3Tau_opt_version(E, t, m, n, tt=0):

    tau = 3
    p = (t-(n-m)) // 2 # Position of leading diagonal
    kval = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    p = pp

    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            if j < kval or j >= kval+t:
                continue

            term_i1_j3 = [
                3 + D[i+2][j+1],
                1 + D[i+1][j+4],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)],
                0 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            term_i2_j3 = [
                3 + D[i+3][j+1],
                2 + D[i+1][j+4],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                0 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            term_i3_j1 = [
                1 + D[i+4][j+1],
                3 + D[i+1][j+2],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                0 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)]
            ]

            term_i3_j2 = [
                2 + D[i+4][j+1],
                3 + D[i+1][j+3],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                0 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                2 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)]
            ]

            term_i3_j3 = [
                3 + D[i+4][j+1],
                3 + D[i+1][j+4],
                0 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                1 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]
            
            D[i+1+1][j+3+1] = secmin(term_i1_j3)
            D[i+2+1][j+3+1] = secmin(term_i2_j3)
            D[i+3+1][j+1+1] = secmin(term_i3_j1)
            D[i+3+1][j+2+1] = secmin(term_i3_j2)
            D[i+3+1][j+3+1] = secmin(term_i3_j3)
        kval+=tau

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())

    return D[m][n]

def ukkonenFill1Tau(E, t, m, n, tt=0):

    p = (t-(n-m)) // 2 # Position of leading diagonal
    kval = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    howmany = p
    p = pp

    @for_range(0, m)
    def _(i):
        start = (kval+i > 0).if_else(i-howmany, 0)
        end = (kval+i+t < n).if_else(kval+i+t, n)
        @for_range(start, end)
        def _(j):
            vals = [
                D[i][j] + E[i][p+(j-i)],
                D[i][j+1] + 1,
                D[i+1][j] + 1
            ]
            D[i+1][j+1] = secmin(vals)
    return D[m][n]

# opt version for 1 tau is same as ukkonenRegular in ukkonen.py

# The Diagonal versions are also opt versions
def ukkonenFill4TauDiagonal(E, t, m, n, tt=0):
    num_boxes_filled = 0
    tau = 4
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    p = pp
    si, sj = 0, 0
    ei, ej = 0, 0
    for diag in range(0, ((((m-1)//tau)-1)*2)+1):
        i, j = si, sj
        l = (si-ei+1)
        fills = []
        for x in range(0, l, tau):

            terms_i1_j4 = [
                4 + D[i+2][j+1],
                1 + D[i+1][j+5],
                3 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)],
                0 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i2_j4 = [
                4 + D[i+3][j+1],
                2 + D[i+1][j+5],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                0 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                1 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i3_j4 = [
                4 + D[i+4][j+1],
                3 + D[i+1][j+5],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+4)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                0 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+4] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                2 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]

            terms_i4_j1 = [
                1 + D[i+5][j+1],
                4 + D[i+1][j+2],
                3 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                0 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)]
            ]

            terms_i4_j2 = [
                2 + D[i+5][j+1],
                4 + D[i+1][j+3],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                0 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                1 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                3 + D[i+1][j+2] + E[i+4][p+(j+2)-(i+4)],
                3 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)],
                3 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)],
                3 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)]
            ]

            terms_i4_j3 = [
                3 + D[i+5][j+1],
                4 + D[i+1][j+4],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                0 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+4][j+1] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                2 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                2 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+1][j+3] + E[i+4][p+(j+3)-(i+4)],
                3 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)],
                3 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                3 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)],
            ]

            terms_i4_j4 = [
                4 + D[i+5][j+1],
                4 + D[i+1][j+5],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                0 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)] + E[i+4][p+(j+2)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+3)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+2)-(i+4)],
                3 + D[i+4][j+1] + E[i+4][p+(j+1)-(i+4)],
                3 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)] + E[i+4][p+(j+3)-(i+4)],
                1 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+4][p+(j+4)-(i+4)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+4][p+(j+4)-(i+4)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+3][p+(j+4)-(i+3)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+1][j+4] + E[i+4][p+(j+4)-(i+4)],
                3 + D[i+1][j+4] + E[i+3][p+(j+4)-(i+3)],
                3 + D[i+1][j+4] + E[i+2][p+(j+4)-(i+2)],
                3 + D[i+1][j+4] + E[i+1][p+(j+4)-(i+1)]
            ]
            num_boxes_filled+=1
            fills.append(((i+1+1, j+4+1), terms_i1_j4))
            fills.append(((i+2+1, j+4+1), terms_i2_j4))
            fills.append(((i+3+1, j+4+1), terms_i3_j4))
            fills.append(((i+4+1, j+1+1), terms_i4_j1))
            fills.append(((i+4+1, j+2+1), terms_i4_j2))
            fills.append(((i+4+1, j+3+1), terms_i4_j3))
            fills.append(((i+4+1, j+4+1), terms_i4_j4))
            # D[i+1+1][j+4+1] = secmin(terms_i1_j4)
            # D[i+2+1][j+4+1] = secmin(terms_i2_j4)
            # D[i+3+1][j+4+1] = secmin(terms_i3_j4)
            # D[i+4+1][j+1+1] = secmin(terms_i4_j1)
            # D[i+4+1][j+2+1] = secmin(terms_i4_j2)
            # D[i+4+1][j+3+1] = secmin(terms_i4_j3)
            # D[i+4+1][j+4+1] = secmin(terms_i4_j4)

            i-=tau
            j+=tau

        # @for_range(0, len(fills))
        # def _(i):
        for i in range(len(fills)):
            D[fills[i][0][0]][fills[i][0][1]] = secmin(fills[i][1])

        si, sj, ei, ej = iterateTauedDiagonal(si, sj, ei, ej, m, n, k, t, tau)

    print("Number of boxes", num_boxes_filled)
    return D[m][n]

def ukkonenFill3TauDiagonal(E, t, m, n, tt=0):
    num_boxes_filled = 0
    tau = 3
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    p = pp
    si, sj = 0, 0
    ei, ej = 0, 0
    for diag in range(0, ((((m-1)//tau)-1)*2)+1):
        i, j = si, sj
        l = (si-ei+1)
        fills = []
        for x in range(0, l, tau):

            terms_i1_j3 = [
                3 + D[i+2][j+1],
                1 + D[i+1][j+4],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)],
                0 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            terms_i2_j3 = [
                3 + D[i+3][j+1],
                2 + D[i+1][j+4],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                0 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                1 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]

            terms_i3_j1 = [
                1 + D[i+4][j+1],
                3 + D[i+1][j+2],
                2 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)],
                0 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)]
            ]

            terms_i3_j2 = [
                2 + D[i+4][j+1],
                3 + D[i+1][j+3],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)],
                0 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                1 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                2 + D[i+1][j+2] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)],
                2 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)]
            ]

            terms_i3_j3 = [
                3 + D[i+4][j+1],
                3 + D[i+1][j+4],
                0 + D[i+1][j+1] + E[i+1][p+(j+1)-(i+1)] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+2][j+1] + E[i+2][p+(j+1)-(i+2)] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+2)-(i+3)],
                2 + D[i+3][j+1] + E[i+3][p+(j+1)-(i+3)],
                1 + D[i+1][j+2] + E[i+2][p+(j+2)-(i+2)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+3][p+(j+3)-(i+3)],
                1 + D[i+1][j+2] + E[i+1][p+(j+2)-(i+1)] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+3] + E[i+3][p+(j+3)-(i+3)],
                2 + D[i+1][j+3] + E[i+2][p+(j+3)-(i+2)],
                2 + D[i+1][j+3] + E[i+1][p+(j+3)-(i+1)]
            ]
            num_boxes_filled+=1
            fills.append(((i+1+1, j+3+1), terms_i1_j3))
            fills.append(((i+2+1, j+3+1), terms_i2_j3))
            fills.append(((i+3+1, j+1+1), terms_i3_j1))
            fills.append(((i+3+1, j+2+1), terms_i3_j2))
            fills.append(((i+3+1, j+3+1), terms_i3_j3))

            i-=tau
            j+=tau
        
        for i in range(len(fills)):
            D[fills[i][0][0]][fills[i][0][1]] = secmin(fills[i][1])

        si, sj, ei, ej = iterateTauedDiagonal(si, sj, ei, ej, m, n, k, t, tau)

    print("Number of boxes", num_boxes_filled)
    return D[m][n]