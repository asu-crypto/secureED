#!/usr/bin/env python3

import sys
sys.path.insert(0, '../secure-edit-distance')

from Compiler.GC.types import sbitint, sbitintvec
from Compiler.types import sint, Array, Matrix

from utils import min_multidim, min_pair
from preprocess import iterateDiagonal_with_p

def find_t_opt(a, b, t: int, modval: int):
    '''
    Demo function for optimization for t.
    Does not have storage/logging.
    '''

    m, n = len(a), len(b)
    lead_diag = (t-(n-m)) // 2
    k = -lead_diag
    p = lead_diag

    cumulative = Matrix(t, 2, sint)
    for i in range(t):
        cumulative[i][0] = abs(p-i) # Add cost to reach that diagonal from current
        cumulative[i][1] = i # Indices secret shared
    
    min_val_cumulated = sint(0)
    si, sj = 0, 0 # Start i, Start j
    ei, ej = 0, 0 # End i, End j

    for iter in range((m+n)-1):

        vi, vj = si, sj
        l = si-ei+1
        for x in range(l):
            cumulative[p+(x*2)][0] = cumulative[p+(x*2)][0] + a[vi].not_equal(b[vj])
            vi-=1
            vj+=1
        
        si, sj, ei, ej, p = iterateDiagonal_with_p(si, sj, ei, ej, p, m, n, k, t)

        if iter!=0 and iter % modval == 0:
            mv = min_multidim(cumulative) # Minimum based on 0th value
            min_val_cumulated = min_val_cumulated + mv[0]
            for i in range(t): 
                cumulative[i][0] = abs(mv[1]-i) # Add cost to reach that diagonal from current
            
    mv = min_multidim(cumulative)
    # Add cost to get back to leading diagonal
    min_val_cumulated = min_val_cumulated + mv[0] + abs(mv[1]-lead_diag)

    return min_val_cumulated

def find_opt_diag(E, t: int, modval: int):

    m, n = len(E), len(E)
    lead_diag = (t-(n-m)) // 2
    k = -lead_diag
    p = lead_diag

    cumulative = Matrix(t, 2, sint)
    for i in range(t):
        cumulative[i][0] = abs(p-i)
        cumulative[i][1] = i

    # Could store cumulated minval if shaping is going to be used
    mvees = Array(((m+n)-1) // modval, sint)
    vv = 0

    min_val_cumulated = sint(0)
    compCount, minCount = 0, 0 # For logging
    si, sj = 0, 0
    ei, ej = 0, 0
    for diagIter in range((m+n)-1):
        vi, vj = si, sj
        l = si-ei+1
        for x in range(l):
            compCount+=1
            cumulative[p+(x*2)][0] = cumulative[p+(x*2)][0] + E[vi][lead_diag+(vj-vi)]
            vi-=1
            vj+=1
        si, sj, ei, ej, p = iterateDiagonal_with_p(si, sj, ei, ej, p, m, n, k, t)

        if diagIter!=0 and diagIter % modval == 0:
            minCount+=1
            fa = min_multidim(cumulative)
            min_val_cumulated = fa[0] + min_val_cumulated
            mvees[vv] = min_val_cumulated
            vv+=1
            for i in range(t):
                cumulative[i][0] = abs(fa[1]-i)
                # cumulative[i][1] = i
            
    # print_ln("C %s", cumulative.reveal())
    minCount+=1
    fa = min_multidim(cumulative)
    min_val_cumulated = fa[0] + min_val_cumulated + abs(fa[1]-lead_diag)
    # print_ln('\nMVEES %s', mvees.reveal())

    return min_val_cumulated, compCount, minCount, mvees

def find_opt_diag_gc(E, t: int, modval: int, bitlengths):
    
    m, n = len(E), len(E)
    lead_diag = (t-(n-m)) // 2
    k = -lead_diag
    p = lead_diag

    gc_i = sbitint.get_type(bitlengths[0]) # Bit length for under max(m,n) values
    gc_t = sbitint.get_type(bitlengths[1]) # Bit length for under t values
    gc_vec = sbitintvec.get_type(bitlengths[0]) # Bit length for under max(m,n) values
    
    cumulative = [0 for _ in range(t)]
    ctarray = [0 for _ in range(t)]
    for i in range(t):
        cumulative[i] = gc_i(abs(p-i))
        ctarray[i] = gc_t(i)
    
    # Could store cumulated minval if shaping is going to be used
    mvees = [0 for _ in range(((m+n)-1) // modval)]
    vv = 0

    min_val_cumulated = gc_i(0)
    compCount, minCount = 0, 0 # For logging
    si, sj = 0, 0
    ei, ej = 0, 0

    for diagIter in range((m+n)-1):
        
        vi, vj = si, sj
        l = si-ei+1
        
        # Parallelized addition
        c1, c2 = [], []
        for x in range(l):
            c1.append(cumulative[p+(x*2)])
            c2.append(E[vi][lead_diag+(vj-vi)])
            vi-=1
            vj+=1
        c3 = (gc_vec(c1) + gc_vec(c2)).elements()
        vi, vj = si, sj
        for x in range(l):
            compCount+=1
            cumulative[p+(x*2)] = c3[x]
            vi-=1
            vj+=1

        # Regular addition
        # for x in range(l):
        #     compCount+=1
        #     cumulative[p+(x*2)] = cumulative[p+(x*2)] + E[vi][diag_ct+(vj-vi)]
        #     vi-=1
        #     vj+=1

        si, sj, ei, ej, p = iterateDiagonal_with_p(si, sj, ei, ej, p, m, n, k, t)

        if diagIter!=0 and diagIter % modval == 0:
            minCount+=1
            fa = min_pair(cumulative, ctarray)
            min_val_cumulated += fa[0]
            mvees[vv] = min_val_cumulated
            vv+=1
            for i in range(t):
                cumulative[i] = gc_i(abs(fa[1]-i))
            
    # print_ln("C %s", cumulative.reveal())
    minCount+=1
    fa = min_pair(cumulative, ctarray)
    min_val_cumulated += fa[0] + abs(fa[1]-lead_diag)
    # print_ln('\nMVEES %s', mvees.reveal())

    return min_val_cumulated, compCount, minCount, mvees