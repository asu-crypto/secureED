#!/usr/bin/env python3

import sys
sys.path.insert(0, '../secure-edit-distance')

from Compiler.types import cint, sint, Array
from Compiler.util import min as secmin
from Compiler.library import print_ln

from preprocess import load_D_ukk, iterateDiagonal, iterateDiagonal_with_p

def ukkonenRegular(E, t, m, n, tt=0):
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_ukk(m, n, t)
    compareCount, minCount = 0, 0
    for i in range(m):
        for j in range(n):
            if j < k or j >= k+t:
                continue
            minCount+=1
            D[i+1][j+1] = secmin([D[i][j]+E[i][pp+(j-i)], D[i][j+1]+1, D[i+1][j]+1])
        k+=1
    return D[m][n], compareCount, minCount

def ukkonenDiagonal(E, t: int, m: int, n: int, tt=0):
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_ukk(m, n, t)
    compareCount, minCount = 0, 0
    si, sj = 0, 0
    ei, ej = 0, 0
    for _ in range((m+n)-1):
        vi, vj = si, sj
        l = si-ei+1
        for x in range(l):
            minCount+=1
            D[vi+1][vj+1] = secmin([D[vi][vj]+E[vi][pp+(vj-vi)], D[vi][vj+1]+1, D[vi+1][vj]+1])
            vi-=1
            vj+=1
        si, sj, ei, ej = iterateDiagonal(si, sj, ei, ej, m, n, k, t)
    return D[m][n], compareCount, minCount

def ukkonenDiagonal_with_ab(a, b, t: int, m: int, n: int):
    '''Uses a, b instead of equality matrix E'''
    p = (t-(n-m)) // 2
    k = -p
    D = load_D_ukk(m, n, t)
    compareCount, minCount = 0, 0
    si, sj = 0, 0
    ei, ej = 0, 0
    for _ in range((m+n)-1):
        vi, vj = si, sj
        l = si-ei+1
        for x in range(l):
            compareCount+=1
            minCount+=1
            D[vi+1][vj+1] = secmin([D[vi][vj]+a[vi].not_equal(b[vj]), D[vi][vj+1]+1, D[vi+1][vj]+1])
            vi-=1
            vj+=1
        si, sj, ei, ej = iterateDiagonal(si, sj, ei, ej, m, n, k, t)
    return D[m][n], compareCount, minCount

# Shape optimization to reduce t as we progress. Works well but leaks some info
def ukkonenDiagonalShape(E, t: int, m: int, n: int, tt=0):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_ukk(m, n, t)
    compareCount, minCount = 0, 0
    expected_leading = [] # Need to fill based on mvees (MP-SPDZ limitations)
    actual_leading = [] # Need to fill based on run with no shaping (MP-SPDZ limitations)
    si, sj = 0, 0
    ei, ej = 0, 0
    flagct = 0
    for tester in range((m+n)-1):
        vi, vj = si, sj
        l = si-ei+1
        flag = False
        for x in range(l):
            minCount+=1
            D[vi+1][vj+1] = secmin([D[vi][vj]+E[vi][pp+(vj-vi)], D[vi][vj+1]+1, D[vi+1][vj]+1])
            if vi == vj:
                compareCount+=1
                if actual_leading[vi]+1 < expected_leading[vi]:
                    flag = True
            vi-=1
            vj+=1
        if flag:
            flagct+=1
            expected_leading[(tester//2)+1:] = [val-2 for val in expected_leading[(tester//2)+1:]]
            t -= 2
            k += 1
        si, sj, ei, ej = iterateDiagonal(si, sj, ei, ej, m, n, k, t)

    vv = Array(m, sint)
    for i in range(1, len(D)):
        vv[i-1] = D[i][i]
    print_ln('\n%s', vv.reveal())
    print("Final T", t)
    return D[m][n], compareCount, minCount

def ukkonenDiagonalMultiProc(a, b, t: int, m: int, n: int):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    D = load_D_ukk(m, n, t)

    si, sj = 0, 0
    ei, ej = 0, 0
    iVals = Array(t, cint)
    jVals = Array(t, cint)
    for _ in range((m+n)-1):
        
        vi, vj = si, sj
        l = si-ei+1
        # iVals = [vi-p for p in range(l)]
        # jVals = [vj+p for p in range(l)]
        for x in range(l):
            iVals[x], jVals[x] = vi, vj
            vi-=1
            vj+=1
        
        # @for_range_multithread(2, l, l)
        # @for_range_parallel(l, l)
        # @for_range(l)
        def _(x):
            vi, vj = iVals[x], jVals[x]
            comp = a[vi].not_equal(b[vj])
            D[vi+1][vj+1] = secmin([D[vi][vj]+comp, D[vi][vj+1]+1, D[vi+1][vj]+1])
            # D[vi+1][vj+1] = (D[vi][vj]+comp).min(D[vi][vj+1]+1).min(D[vi+1][vj]+1)

        si, sj, ei, ej = iterateDiagonal(si, sj, ei, ej, m, n, k, t)

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]

def ukkonenDiagonalShape_opt(E, t: int, m: int, n: int, tt, modval):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_ukk(m, n, t)
    compCount, minCount = 0, 0
    si, sj = 0, 0
    ei, ej = 0, 0
    cumulative = Array(t, sint)
    actual_mvees = Array(((m+n)-1) // modval, sint)
    vv, sub = 0, 0
    expected = [] # Need to fill based on mvees (MP-SPDZ limitations)
    actual = [] # Need to fill based on run with no shaping (MP-SPDZ limitations)
    for iter in range((m+n)-1):
        vi, vj = si, sj
        l = si-ei+1
        for x in range(l):
            minCount+=1
            D[vi+1][vj+1] = secmin([D[vi][vj]+E[vi][pp+(vj-vi)], D[vi][vj+1]+1, D[vi+1][vj]+1])
            if iter != 0 and (iter % modval == 0 or iter % modval == modval-1):
                cumulative[p+(x*2)] = D[vi+1][vj+1]
            vi-=1
            vj+=1
        si, sj, ei, ej, p = iterateDiagonal_with_p(si, sj, ei, ej, p, m, n, k, t)
        
        if iter!=0 and iter % modval == 0:
            minCount+=1
            actual_mvees[vv] = secmin(cumulative)
            if actual[vv]+1 < expected[vv]-sub:
                print("hit", iter)
                t -= 2
                k += 1
                sub += 2
            vv+=1

    # print_ln('\nACVEES%s', actual_mvees.reveal())
    print("Final T", t)
    return D[m][n], compCount, minCount