#!/usr/bin/env python3

import sys
sys.path.insert(0, '../secure-edit-distance')

from Compiler.util import min as secmin
from Compiler.GC.types import sbitint, sbitintvec
from Compiler.library import for_range

from preprocess import load_D_ukk_gc, iterateDiagonal

def ukkonenRegular_gc(E, t, m, n, tt, bitlengths):
    gc_vec = sbitintvec.get_type(bitlengths[0])
    gc_vec_bool = sbitintvec.get_type(2)
    gc_bool = sbitint.get_type(2)
    p = (t-(n-m)) // 2 # Position of leading diagonal
    kval = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_ukk_gc(m, n, t, bitlengths[0])
    compareCount, minCount = 0, 0 # For logging
    @for_range(0, m)
    def _(i):
        start = (kval+i > 0).if_else(i-p, 0)
        end = (kval+i+t < n).if_else(kval+i+t, n)
        @for_range(start, end)
        def _(j):
            pd = gc_vec([D[i][j], D[i][j+1], D[i+1][j]])
            qd = gc_vec_bool([E[i][pp+(j-i)], gc_bool(1), gc_bool(1)])
            D[i+1][j+1] = secmin((pd + qd).elements())
    return D[m][n], compareCount, minCount

def ukkonenRegular_gc_compiled(E, t, m, n, tt, bitlengths):
    gc_vec = sbitintvec.get_type(bitlengths[0])
    gc_vec_bool = sbitintvec.get_type(2)
    gc_bool = sbitint.get_type(2)
    p = (t-(n-m)) // 2 # Position of leading diagonal
    kval = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_ukk_gc(m, n, t, bitlengths[0])
    compareCount, minCount = 0, 0 # For logging
    for i in range(m):
        for j in range(n):
            if j < kval or j >= kval+t:
                continue
            minCount+=1
            pd = gc_vec([D[i][j], D[i][j+1], D[i+1][j]])
            qd = gc_vec_bool([E[i][pp+(j-i)], gc_bool(1), gc_bool(1)])
            D[i+1][j+1] = secmin((pd + qd).elements())
        kval+=1
    return D[m][n], compareCount, minCount

def ukkonenDiagonal_gc(E, t: int, m: int, n: int, tt, bitlengths):
    
    gc_vec = sbitintvec.get_type(bitlengths[0])
    gc_vec_bool = sbitintvec.get_type(2)
    gc_bool = sbitint.get_type(2)
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2
    D = load_D_ukk_gc(m, n, t, bitlengths[0])

    compareCount, minCount = 0, 0
    si, sj = 0, 0
    ei, ej = 0, 0
    for _ in range((m+n)-1):

        vi, vj = si, sj
        l = si-ei+1

        # Parallelized addition
        c1, c2 = [], []
        for x in range(l):
            c1.append(D[vi][vj])
            c1.append(D[vi][vj+1])
            c1.append(D[vi+1][vj])
            c2.append(E[vi][pp+(vj-vi)])
            c2.append(gc_bool(1))
            c2.append(gc_bool(1))
            vi-=1
            vj+=1
        c3 = (gc_vec(c1) + gc_vec_bool(c2)).elements()

        # Parallelized minimum
        # g1, g2, g3 = [], [], []
        # for x in range(l):
        #     g1.append(c3[x*3])
        #     g2.append(c3[(x*3)+1])
        #     g3.append(c3[(x*3)+2])
        # g1, g2, g3 = gc_vec(g1), gc_vec(g2), gc_vec(g3)
        # compare = gc_vec(g1.min(g2)).elements()
        # inv_compare = gc_vec(g2 < g1)
        # intermediate = (g1*compare + g2*inv_compare)
        # compare = gc_vec(g3 <= intermediate)
        # inv_compare = gc_vec(intermediate < g3)
        # minvals = (g3*compare + intermediate*inv_compare).elements()

        vi, vj = si, sj
        for x in range(l):
            minCount+=1
            D[vi+1][vj+1] = secmin([c3[x*3], c3[(x*3)+1], c3[(x*3)+2]])
            # D[vi+1][vj+1] = compare[x]
            vi-=1
            vj+=1

        # Regular addition and regular minimum
        # for x in range(l):
        #     D[vi+1][vj+1] = secmin([D[vi][vj]+E[vi][pp+(vj-vi)], D[vi][vj+1]+1, D[vi+1][vj]+1])
        #     vi-=1
        #     vj+=1
            
        si, sj, ei, ej = iterateDiagonal(si, sj, ei, ej, m, n, k, t)

    return D[m][n], compareCount, minCount

def ukkonenDiagonal_gc_with_ab(a, b, t: int, m: int, n: int, bitlengths):
    
    gc_vec = sbitintvec.get_type(bitlengths[0])
    gc_vec_bool = sbitintvec.get_type(2)
    gc_bool = sbitint.get_type(2)
    
    p = (t-(n-m)) // 2
    k = -p
    D = load_D_ukk_gc(m, n, t, bitlengths[0])

    compareCount, minCount = 0, 0 # For logging
    si, sj = 0, 0
    ei, ej = 0, 0
    for _ in range((m+n)-1):

        vi, vj = si, sj
        l = si-ei+1

        # Parallelized addition
        c1, c2 = [], []
        for x in range(l):
            c1.append(D[vi][vj])
            c1.append(D[vi][vj+1])
            c1.append(D[vi+1][vj])
            compareCount+=1
            c2.append(gc_bool(a[vi].not_equal(b[vj])))
            c2.append(gc_bool(1))
            c2.append(gc_bool(1))
            vi-=1
            vj+=1
        c3 = (gc_vec(c1) + gc_vec_bool(c2)).elements()
        vi, vj = si, sj
        for x in range(l):
            minCount+=1
            D[vi+1][vj+1] = secmin([c3[x*3], c3[(x*3)+1], c3[(x*3)+2]])
            vi-=1
            vj+=1

        si, sj, ei, ej = iterateDiagonal(si, sj, ei, ej, m, n, k, t)

    return D[m][n], compareCount, minCount
