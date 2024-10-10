#!/usr/bin/env python3

import sys
sys.path.insert(0, '../secure-edit-distance')

from Compiler.util import min as secmin
from Compiler.types import sint, Array
from Compiler.library import print_ln

from ways import readWays
from preprocess import load_D_box, iterateTauedDiagonal, iterateTauedDiagonal_with_p
from box_find import boxFindRetMin, boxFindRetAll_with_ab, boxFindRetAll

def fillBox_with_ab(D, a, b, i, j, p, ww, k, t):
    for x in range(len(ww)):
        ix = ww[x][0][0]
        iy = ww[x][0][1]
        # If start cell beyond diagonal boundary, skip.
        if j+iy < k+ix or j+iy >= k+t+ix:
            continue
        D[i+ix+1][j+iy+1] = secmin(boxFindRetAll_with_ab(D, a, b, ww[x], i, j))

def fillBox(D, E, i, j, p, ww, k, t):
    cc = 0 # For logging number of minimums
    for x in range(len(ww)):
        ix = ww[x][0][0]
        iy = ww[x][0][1]
        # If start cell beyond diagonal boundary, skip.
        if j+iy < k+ix or j+iy >= k+t+ix:
            continue
        cc+=1
        D[i+ix+1][j+iy+1] = secmin(boxFindRetAll(D, E, ww[x], i, j, p, t))
    return cc

def ukkonenFill(E, t: int, tau: int, m: int, n: int, tt=0):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    minCount = 0 # for logging
    waysInBox = readWays(tau)
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            if j < k or j >= k+t:
                continue
            minCount+= fillBox(D, E, i, j, pp, waysInBox, k, t)
        k+=tau

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n], minCount

def ukkonenFillShape(E, t: int, tau: int, m: int, n: int, tt=0):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    minCount, compCount = 0, 0
    waysInBox = readWays(tau)
    expected_leading = [] # Need to fill based on mvees (MP-SPDZ limitations)
    actual_leading = [] # Need to fill based on run with no shaping (MP-SPDZ limitations)
    for i in range(0, m-tau, tau):
        flag = False
        for j in range(0, n-tau, tau):
            if j < k or j >= k+t:
                continue
            minCount += fillBox(D, E, i, j, pp, waysInBox, k, t)
            if i == j:
                compCount+=1
                if actual_leading[i]+1 < expected_leading[i]:
                    flag = True
        if flag:
            expected_leading[i+1:] = [v-2 for v in expected_leading[i+1:]]
            t -= 2
            k += 1
        k+=tau

    vv = Array(m, sint)
    for i in range(1, len(D)):
        vv[i-1] = D[i][i]
    print_ln('%s', vv.reveal())
    print("Final T", t)
    return D[m][n], compCount, minCount

def ukkonenFillDiagonal(E, t, tau, m, n, tt=0):
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    waysInBox = readWays(tau)
    si, sj = 0, 0
    ei, ej = 0, 0
    for diag in range(0, ((((m-1)//tau)-1)*2)+1):
        vi, vj = si, sj
        l = (si-ei+1)
        for x in range(0, l, tau):
            fillBox(D, E, vi, vj, pp, waysInBox, vi+k, t)
            vi-=tau
            vj+=tau
        si, sj, ei, ej = iterateTauedDiagonal(si, sj, ei, ej, m, n, k, t, tau)
    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]

def reverseIterDiagonal_with_p(si, sj, ei, ej, p, k, t):
    if si <= -k or (k+si+1) > sj:
        si-=1
        p+=1
    else:
        sj-=1
        p-=1
    if ej <= -k or (k+ei+t) <= (ej+1):
        ej-=1
    else:
        ei-=1
    return si, sj, ei, ej, p

def loadCumulative(D, cumulative, p, si, sj, ei, ej, k, t):
    vi, vj = si, sj
    l = (si-ei+1)
    for x in range(l):
        cumulative[p+(x*2)] = D[vi+1][vj+1]
        vi-=1
        vj+=1
    si, sj, ei, ej, p = reverseIterDiagonal_with_p(si, sj, ei, ej, p, k, t)
    vi, vj = si, sj
    l = (si-ei+1)
    for x in range(l):
        cumulative[p+(x*2)] = D[vi+1][vj+1]
        vi-=1
        vj+=1
    # print_ln('THISS %s', cumulative.reveal())

def ukkonenFillDiagonalShape_opt(E, t: int, tau: int, m: int, n: int, tt, modval):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_box(E, m, n, t, p, pp)
    cumulative = Array(t, sint)
    cumulative.assign_all(max(m, n))
    minCount, compCount = 0, 0
    waysInBox = readWays(tau)
    si, sj = 0, 0
    ei, ej = 0, 0
    diag_ct = (((m-1)//tau)*2)-1
    actual_mvees = Array(diag_ct*tau // modval, sint)
    vv, sub = 0, 0
    expected = [] # Need to fill based on mvees (MP-SPDZ limitations)
    actual = [] # Need to fill based on run with no shaping (MP-SPDZ limitations)
    for diag in range(0, diag_ct):
        vi, vj = si, sj
        l = (si-ei+1)
        for x in range(0, l, tau):
            minCount += fillBox(D, E, vi, vj, pp, waysInBox, vi+k, t)
            vi-=tau
            vj+=tau
        si, sj, ei, ej, p = iterateTauedDiagonal_with_p(si, sj, ei, ej, m, n, k, t, tau, p)

        if ((diag*tau)+tau) % (modval) == 0:
            loadCumulative(D, cumulative, p, si, sj, ei, ej, k, t)
            minCount+=1
            actual_mvees[vv] = secmin(cumulative)
            if actual[vv]+1 < expected[vv]-sub:
                print("hit", diag)
                t -= 2
                k += 1
                sub += 2
            vv+=1
    
    # print_ln('\nACVEES%s', actual_mvees.reveal())
    print("Final T", t)
    return D[m][n], compCount, minCount

# Match optimization: Did not perform better than C compiler's compiletime optimization

def ukkonenFillMatchAll(E, leaks: list, t: int, tau: int, m: int, n: int, tt=0):

    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)

    waysInBox = readWays(tau)
    waysIfMatch = readWays(tau, True)
    match_ct, box_ct = 0, 0
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):

            if j < k or j >= k+t:
                continue

            if not leaks[box_ct]:
                match_ct+=1
                fillBox(D, E, i, j, pp, waysIfMatch)
            else:
                fillBox(D, E, i, j, pp, waysInBox)

            box_ct+=1
        k+=tau

    print("Matched Boxes: %d / %d" % (match_ct, box_ct))
    print("Percent: {:.3f} %".format(((match_ct/box_ct)*100)))

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]

def ukkonenFillMatchLeading(E, leaks, t, tau, m, n, tt=0):

    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)

    waysInBox = readWays(tau)
    waysIfMatch = readWays(tau, True)
    v, match_ct, box_ct = 0, 0, 0
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):

            if j < k or j >= k+t:
                continue

            match_flag = False
            if i == j:
                if not leaks[v]:
                    match_ct+=1
                    match_flag = True
                v+=1
            
            if match_flag:
                fillBox(D, E, i, j, pp, waysIfMatch)
            else:
                fillBox(D, E, i, j, pp, waysInBox)

            box_ct+=1
        k+=tau

    print("Matched Boxes: %d / %d" % (match_ct, box_ct))
    print("Percent: {:.3f} %".format(((match_ct/box_ct)*100)))

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]