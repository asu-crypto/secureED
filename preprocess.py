#!/usr/bin/env python3

from Compiler.GC.types import sbitint, sbitintvec, sbits
from Compiler.types import cint, sint, sintbit, Array, Matrix
from Compiler.library import for_range_opt, for_range_parallel, print_ln
from Compiler.util import min as secmin

sb = sbits.get_type(1)

def xor_nucleotids(N1, N2):
    xor = Array(2, sb)
    xor[0] = N1[0] ^ N2[0]
    xor[1] = N1[1] ^ N2[1]
    return xor

def or_nucleotid(N):
    return N[0].bit_or(N[1])

def equal_nucleotids(N1, N2):
    return or_nucleotid(xor_nucleotids(N1, N2)).bit_not()

def bitEquality(a, b, t):
    p = ((t-(len(b)-len(a))) // 2)
    k = -p
    E = Matrix(len(a), t+4, sintbit)
    E.assign_all(1)
    for i in range (0, len(a)):
        for j in range(max(k+i, 0), min(k+t+i, len(b))):
            E[i][p+(j-i)] = sintbit(equal_nucleotids(a[i], b[j]).bit_not())
    # print_ln('%s', E.reveal())
    return E, 0

def equalityParallel(a, b, t: int):
    '''Compute the equality matrix parallely'''
    p = ((t-(len(a)-len(b))) // 2)
    k = -p
    E = Matrix(len(a), t, sintbit)
    @for_range_opt(0, len(a))
    def _(i):
        @for_range_opt((k+i).max(0), (k+t+i).min(len(b)))
        def _(j):
            E[i][p+(j-i)] = a[i].not_equal(b[j])
    return E

def equality(a, b, t: int):
    '''Compute the equality matrix'''
    p = ((t-(len(b)-len(a))) // 2)
    k = -p
    E = Matrix(len(a), t, sintbit)
    E.assign_all(1)
    compareCount = 0
    for i in range (0, len(a)):
        for j in range(max(k+i, 0), min(k+t+i, len(b))):
            compareCount+=1
            E[i][p+(j-i)] = a[i].not_equal(b[j])
    return E, compareCount

def equality_gc(a, b, t: int):
    '''Compute the equality matrix parallelized'''
    p = ((t-(len(b)-len(a))) // 2)
    k = -p
    gc_bool = sbitint.get_type(2) # Bit length for 0/1
    gc_vec_bool = sbitintvec.get_type(2)
    E = [[gc_bool(0) for _ in range(t)] for _ in range(len(a))]
    compareCount = 0
    # Parallelized comparison
    seq1 = [a[i] for i in range(len(a)) for _ in range(max(k+i, 0), min(k+t+i, len(b)))]
    seq2 = [b[j] for i in range(len(a)) for j in range(max(k+i, 0), min(k+t+i, len(b)))]
    comp_array = gc_vec_bool(seq1).not_equal(gc_vec_bool(seq2)).elements()
    vv = 0
    for i in range (0, len(a)):
        for j in range(max(k+i, 0), min(k+t+i, len(b))):
            compareCount+=1
            # E[i][p+(j-i)] = gc_bool(a[i].not_equal(b[j])) # Regular comparison
            E[i][p+(j-i)] = gc_bool(comp_array[vv])
            vv+=1
    return E, compareCount

def load_D_box(E, m: int, n: int, t: int, p: int, pp):
    mv = max(m, n)
    D = Matrix(m+1, n+1, sint)
    D.assign_all(sint(mv))
    for i in range(0, m+1):
        D[i][0] = i
    for j in range(1, n+1):
        D[0][j] = j
    for i in range(0, p+1):
        D[i+1][1] = secmin([D[i][0]+E[i][pp-i], D[i][1]+1, D[i+1][0]+1])
    for j in range(1, t-p):
        D[1][j+1] = secmin([D[0][j]+E[0][pp+j], D[1][j]+1, D[0][j+1]+1])
    #print_ln('%s', D.reveal())
    return D

def load_D_box_with_ab(a, b, m: int, n: int, t: int, p: int):
    mv = max(m, n)
    D = Matrix(m+1, n+1, sint)
    D.assign_all(sint(mv))
    for i in range(0, m+1):
        D[i][0] = i
    for j in range(1, n+1):
        D[0][j] = j
    for i in range(0, p+1):
        D[i+1][1] = (D[i][0]+a[i].not_equal(b[0])).min(D[i][1]+1).min(D[i+1][0]+1)
    for j in range(1, t-p):
        D[1][j+1] = (D[0][j]+a[0].not_equal(b[j])).min(D[1][j]+1).min(D[0][j+1]+1)
    #print_ln('%s', D.reveal())
    return D

def load_D_ukk(m: int, n: int, t):
    '''Loads the initial D matrix'''
    mv = max(m, n)
    D = Matrix(m+1, n+1, sint)
    D.assign_all(mv)
    for i in range(0, t):
        D[i][0] = sint(i)
    for j in range(1, t):
        D[0][j] = sint(j)
    return D

def load_D_ukk_gc(m: int, n: int, t, bitlen):
    '''Loads the initial D matrix'''
    mv = max(m, n)
    si32 = sbitint.get_type(bitlen)
    # D = Matrix(m+1, n+1, si32)
    # for i in range(m+1):
    #     for j in range(n+1):
    #         D[i][j] = si32(mv)
    D = [[si32(mv) for _ in range(n+1)] for _ in range(m+1)]
    for i in range(0, t):
        D[i][0] = si32(i)
    for j in range(1, t-1):
        D[0][j] = si32(j)
    return D

def iterateDiagonal(si, sj, ei, ej, m, n, k, t):
    if si == m-1 or (k+si+1) > sj:
        sj+=1
    else:
        si+=1
    if ej == n-1 or (k+ei+t) <= (ej+1):
        ei+=1
    else:
        ej+=1
    return si, sj, ei, ej

def iterateDiagonal_with_p(si, sj, ei, ej, p, m, n, k, t):
    if si == m-1 or (k+si+1) > sj:
        sj+=1
        p+=1
    else:
        si+=1
        p-=1
    if ej == n-1 or (k+ei+t) <= (ej+1):
        ei+=1
    else:
        ej+=1
    return si, sj, ei, ej, p

def iterateTauedDiagonal(si, sj, ei, ej, m, n, k, t, tau):
    if si == m-1-tau or (k+si+tau) > sj:
        sj+=tau
    else:
        si+=tau
    if ej == n-1-tau or (k+ei+t) <= (ej+tau):
        ei+=tau
    else:
        ej+=tau
    return si, sj, ei, ej

def iterateTauedDiagonal_with_p(si, sj, ei, ej, m, n, k, t, tau, p):
    if si == m-1-tau or (k+si+tau) > sj:
        sj+=tau
        p+=tau
    else:
        si+=tau
        p-=tau
    if ej == n-1-tau or (k+ei+t) <= (ej+tau):
        ei+=tau
    else:
        ej+=tau
    return si, sj, ei, ej, p

# Ideas from the past: (Not used anymore)

def load_leaks(E, t: int, tau: int, m: int, n: int):
    p = ((t-(n-m)) // 2)
    k = -p
    leakedComparisons = Array((t//tau)*(n//tau) + m, sintbit)
    v = 0
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            if j < k or j >= k+t:
                continue
            test = sintbit(0)
            for tx in range(1, tau+1):
                test = test | E[i+tx][j+tx-max(0, (-p+i+tx))]
            leakedComparisons[v] = test
            v+=1
        k+=tau
    return leakedComparisons

def load_leaks_multi_proc(E, t: int, tau: int, m: int, n: int):
    p = ((t-(n-m)) // 2)
    k = -p
    kd = p % tau
    leakedComparisons = Matrix(m, t+tau, sintbit)
    kstart = Array(m, cint)
    kend = Array(m, cint)
    for i in range(0, m, tau):
        kstart[i] = 0 if k+i < 0 else k+i+kd
        kend[i] = min(n-tau, k+i+t)
    @for_range_opt(0, m-tau, tau)
    def _(i):
        @for_range_opt(kstart[i], kend[i], tau)
        def __(j):
            test = sintbit(0)
            for tx in range(1, tau+1):
                test = test | E[i+tx][j+tx-((-p+i+tx).max(0))]
            leakedComparisons[i][j-((-p+i).max(0))] = test
    return leakedComparisons

def load_leaks_leading(a, b, tau: int, m: int, n: int):
    leakedComparisons = Array(max(m, n)//tau, sintbit)
    v = 0
    for i in range(0, m-tau, tau):
        test = sintbit(0)
        for tx in range(1, tau+1):
            test = test | a[i+tx].not_equal(b[i+tx])
        leakedComparisons[v] = test
        v+=1
    return leakedComparisons

def cumulated_leaks_leading(a, b, tau: int, m: int, n: int):
    leakedComparisons = Array(max(m, n)//tau, sint)
    v = 0
    for i in range(-1, m-tau, tau):
        test = sintbit(0)
        for tx in range(1, tau+1):
            test = test | a[i+tx].not_equal(b[i+tx])
        # Compares [-1] for first. Beware when parallelising.
        leakedComparisons[v] = leakedComparisons[v-1] + test
        v+=1
    return leakedComparisons

def load_mx_leading(a, b, tau: int, m: int, n: int):
    mx = sint(0)
    for i in range(0, m-tau, tau):
        test = sintbit(0)
        for tx in range(1, tau+1):
            test = test | a[i+tx].not_equal(b[i+tx])
        mx += test
    return mx

def find_opt_endonly(a, b, t: int, m: int, n: int):
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    cumulative = Array(t, sint)
    for i in range(m):
        start = max(k+i, 0)
        idx = max(p,0)
        for j in range(start, min(k+t+i, len(b))):
            cumulative[j+idx-start] = cumulative[j+idx-start] + a[i].not_equal(b[j])
        p-=1
    p = (t-(n-m)) // 2
    for i in range(t):
        cumulative[i] = cumulative[i]+(2*abs(p-i))
    # print_ln("%s", cumulative.reveal())
    return secmin(cumulative)
