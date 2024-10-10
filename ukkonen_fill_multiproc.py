
from preprocess import load_D_box
from ways import readWays, readWaysMultiProcMin, readWaysMultiProcAll
from box_find import boxFindRetMin, boxFindRetAll_with_ab, boxFindMultiProc, boxFindRetAll

from Compiler.util import min as secmin
from Compiler.types import cint, sint, Array, Matrix
from Compiler.library import for_range, for_range_multithread, for_range_opt, for_range_parallel

'''
Uses MP-SPDZ's runtime multiprocessing. Was found to be inefficient for our use case compared to 
compiletime multiprocessing done by the C compiler. Retaining for any future updates.
'''

def fillBoxMultiProcMin(D, i, j, mvees, idxes, ww_size):
    # @for_range(ww_size)
    @for_range_parallel(ww_size, cint(ww_size))
    # @for_range_multithread(1, ww_size, ww_size)
    # @for_range_opt(cint(0), cint(ww_size))
    # for x in range(ww_size):
    def f(x):
        # print("ONE", gg, idxes[x][0])
        ix = idxes[x][0]
        iy = idxes[x][1]
        D[i+ix+1][j+iy+1] = secmin(mvees[x])
        # print_ln("%s", x)

def ukkonenFillMultiProcMin(E, t, tau, m, n, tt=0):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)

    waysInBox = readWays(tau)
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            
            if j < k or j >= k+t:
                continue
            
            mvees = Matrix(len(waysInBox), len(waysInBox[0])-1, sint)
            for x in range(len(waysInBox)):
                mvees[x] = boxFindRetAll_with_ab(D, E, waysInBox[x], i, j)

            # print_ln("%s %s", mvees, mvees1)

            # minvals = list(pool.map(secmin, mvees))
            a = Array(len(waysInBox), sint)
            # for x in range(len(waysInBox)):
            @for_range_parallel(len(waysInBox), len(waysInBox))
            def _(x):
                # a[x] = secmin(mvees[x])
            # for x in range(len(waysInBox)):
            # @for_range(len(waysInBox))
            # @for_range_parallel(len(waysInBox), len(waysInBox))
            # def _(x):
                # ix = idxes[x][0]
                # iy = idxes[x][1]
                ix = waysInBox[x][0][0]
                iy = waysInBox[x][0][1]
                D[i+ix+1][j+iy+1] = secmin(mvees[x])
            # fillBoxMultiProcMin(D, i, j, mvees, idxes, len(waysInBox))
            
        k+=tau

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]

def ukkonenFillMultiProcMinMatchAll(E, leaks: list, t, tau, m, n, tt):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)

    idxes, waysInBox = readWaysMultiProcMin(tau)
    _, waysIfMatch = readWaysMultiProcMin(tau, True)
    match_ct, box_ct = 0, 0
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            
            if j < k or j >= k+t:
                continue

            if not leaks[box_ct]:
                match_ct+=1
                mvees = Matrix(len(waysInBox), len(waysIfMatch[0])-1, sint)
                for x in range(0, len(waysInBox)):
                    mvees[x] = boxFindRetAll(D, E, waysIfMatch[x], i, j, pp)
            else:
                mvees = Matrix(len(waysInBox), len(waysInBox[0])-1, sint)
                for x in range(0, len(waysInBox)):
                    mvees[x] = boxFindRetAll(D, E, waysInBox[x], i, j, pp)
            
            fillBoxMultiProcMin(D, E, i, j, mvees, idxes, len(waysInBox))
            
            box_ct+=1
        k+=tau

    print("Matched Boxes: %d / %d" % (match_ct, box_ct))
    print("Percent: {:.3f} %".format(((match_ct/box_ct)*100)))

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]

def ukkonenFillMultiProcMinMatchLeading(E, leaks: list, t, tau, m, n, tt=0):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)

    idxes, waysInBox = readWaysMultiProcMin(tau)
    _, waysIfMatch = readWaysMultiProcMin(tau, True)
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
                mvees = Matrix(len(waysInBox), len(waysIfMatch[0])-1, sint)
                for x in range(0, len(waysInBox)):
                    mvees[x] = boxFindRetAll(D, E, waysIfMatch[x], i, j, pp)
            else:
                mvees = Matrix(len(waysInBox), len(waysInBox[0])-1, sint)
                for x in range(0, len(waysInBox)):
                    mvees[x] = boxFindRetAll(D, E, waysInBox[x], i, j, pp)
            
            fillBoxMultiProcMin(D, E, i, j, mvees, idxes, len(waysInBox))
            
            box_ct+=1
        k+=tau

    print("Matched Boxes: %d / %d" % (match_ct, box_ct))
    print("Percent: {:.3f} %".format(((match_ct/box_ct)*100)))

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]

def fillBoxMultiProcAll(D, E, i, j, p, ww):
    # @for_range_opt(0, len(ww))
    @for_range_parallel(len(ww), len(ww))
    def f(x):
        ix = ww[x][0][0]
        iy = ww[x][0][1]
        D[i+ix+1][j+iy+1] = boxFindMultiProc(D, E, ww[x], i, j, p)

def ukkonenFillMultiProcAll(E, t, tau, m, n, tt=0):
    
    p = (t-(n-m)) // 2 # Position of leading diagonal
    k = -p # Mapping of first element in array
    pp = p
    if tt!=0:
        pp = (tt-(n-m)) // 2 # E might be created with different t
    D = load_D_box(E, m, n, t, p, pp)
    
    waysInBox = readWaysMultiProcAll(tau)
    for i in range(0, m-tau, tau):
        for j in range(0, n-tau, tau):
            if j < k or j >= k+t:
                continue
            fillBoxMultiProcAll(D, E, i, j, pp, waysInBox)
        k+=tau

    # for i in range(len(D)):
    #     print_ln('%s', D[i].reveal())
    return D[m][n]