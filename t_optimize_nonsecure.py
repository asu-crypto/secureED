#!/usr/bin/env python3

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

def find_opt_nonsecure(a, b, t: int, modval: int):

    m, n = len(a), len(b)
    lead_diag = (t-(n-m)) // 2
    end_diag = lead_diag+(n-m)
    k = -lead_diag
    p = lead_diag

    cumulative = [0 for _ in range(t)]
    for i in range(t):
        cumulative[i] = abs(lead_diag-i)

    min_val_cumulated = 0
    compCount, minCount = 0, 0 # For logging
    si, sj = 0, 0
    ei, ej = 0, 0
    for diagIter in range((m+n)-1):
        vi, vj = si, sj
        l = si-ei+1
        for x in range(l):
            compCount+=1
            comp = 0 if a[vi] == b[vj] else 1
            cumulative[p+(x*2)] = cumulative[p+(x*2)] + comp
            vi-=1
            vj+=1
        si, sj, ei, ej, p = iterateDiagonal_with_p(si, sj, ei, ej, p, m, n, k, t)

        if diagIter != 0 and diagIter % modval == 0 and diagIter != (m+n)-1:
            minCount+=1
            mv = min(cumulative)
            midx = cumulative.index(mv)
            min_val_cumulated = mv + min_val_cumulated
            for i in range(t):
                cumulative[i] = abs(midx-i)
            
    minCount+=1
    for i in range(t):
        cumulative[i] = abs(end_diag-i)
    mv = min(cumulative)
    min_val_cumulated = min_val_cumulated + mv
    return min_val_cumulated, compCount, minCount