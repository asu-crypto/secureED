#!/usr/bin/env python3

def fetch_paths(x, y):

    ways = []
    ways.append([y, x, 0])
    ways.append([x, 0, y])

    def possible_minimal_paths(w, shift, dtrav, rtrav, i, j):
        if i == x and j == y:
            ways.append(w)
            return
        dist = min(x-i, y-j)
        if j != 0 and dtrav == 0:
            if shift > 0:
                tw = w.copy()
                tw[0]+=1
                possible_minimal_paths(tw, shift-1, 0, 2, i+1, j)
            elif dist >= 3:
                tw = w.copy()
                tw[0]+=1
                possible_minimal_paths(tw, shift-1, 0, 2, i+1, j)
        if i != 0 and rtrav == 0:
            if shift < 0:
                tw = w.copy() 
                tw[0]+=1
                possible_minimal_paths(tw, shift+1, 2, 0, i, j+1)
            elif dist >= 3:
                tw = w.copy()
                tw[0]+=1
                possible_minimal_paths(tw, shift+1, 2, 0, i, j+1)
        if i < x and j < y:
            tw = w.copy() 
            tw.append(i+1)
            tw.append(j+1)
            possible_minimal_paths(tw, shift, max(0, dtrav-1), max(rtrav-1, 0), i+1, j+1)
    for i in range(0, x):
        shift = (x-i)-y
        possible_minimal_paths([0, i, 0], shift, 0, 0, i, 0)
    for j in range(1, y):
        shift = x-(y-j)
        possible_minimal_paths([0, 0, j], shift, 0, 0, 0, j)

    return ways

def minimal_ways_in_box(tau):
    
    path_ct = {2: 8, 3: 16, 4: 36, 5:85, 6:213, 7:549, 8:1447, 9:3881, 10:10549}
    waysInBox = [[[-4000 for _ in range((tau*2)+3)] for _ in range(path_ct[tau])] for _ in range((tau*2)-1)]
    v = 0
    for x in range(1, tau):
        ww = fetch_paths(x, tau)
        waysInBox[v][0][:2] = [x, tau]
        for s in range(path_ct[tau]-1):
            if s >= len(ww):
                waysInBox[v][s+1][0] = 4090
            else:
                waysInBox[v][s+1][:len(ww[s])] = ww[s]
        v+=1
    for y in range(1, tau+1):
        ww = fetch_paths(tau, y)
        waysInBox[v][0][:2] = [tau, y]
        for s in range(path_ct[tau]-1):
            if s >= len(ww):
                waysInBox[v][s+1][0] = 4090
            else:
                waysInBox[v][s+1][:len(ww[s])] = ww[s]
        v+=1
    return waysInBox

def fetch_paths_if_match(x, y):

    ways = []

    def possible_minimal_paths(w, shift, dtrav, rtrav, i, j):
        if i == x and j == y:
            ways.append(w)
            return
        dist = min(x-i, y-j)
        if j != 0 and dtrav == 0:
            if shift > 0:
                tw = w.copy()
                tw[0]+=1
                possible_minimal_paths(tw, shift-1, 0, 2, i+1, j)
            elif dist >= 3 and i!=j:
                tw = w.copy()
                tw[0]+=1
                possible_minimal_paths(tw, shift-1, 0, 2, i+1, j)
        if i != 0 and rtrav == 0:
            if shift < 0:
                tw = w.copy() 
                tw[0]+=1
                possible_minimal_paths(tw, shift+1, 2, 0, i, j+1)
            elif dist >= 3 and i!=j:
                tw = w.copy()
                tw[0]+=1
                possible_minimal_paths(tw, shift+1, 2, 0, i, j+1)
        if i < x and j < y:
            tw = w.copy()
            tw.append(i+1)
            tw.append(j+1)
            possible_minimal_paths(tw, shift, max(0, dtrav-1), max(rtrav-1, 0), i+1, j+1)
    
    if x == y:
        ways.append([0, 0, 0])
        return ways
    ways.append([abs(x-y), 0, 0])
    if y < x:
        for i in range(1, x):
            shift = (x-i)-y
            possible_minimal_paths([0, i, 0], shift, 0, 0, i, 0)
        ways.append([y, x, 0])
    if x < y:
        for j in range(1, y):
            shift = x-(y-j)
            possible_minimal_paths([0, 0, j], shift, 0, 0, 0, j)
        ways.append([x, 0, y])

    return ways

def minimal_ways_in_match_box(tau):

    path_ct = {2: 3, 3: 5, 4: 10, 5:22, 6:49, 7:111, 8:258, 9:617, 10:1515}
    waysInBox = [[[-4000 for _ in range((tau*2)+3)] for _ in range(path_ct[tau])] for _ in range((tau*2)-1)]
    v = 0
    for x in range(1, tau):
        ww = fetch_paths_if_match(x, tau)
        waysInBox[v][0][:2] = [x, tau]
        for s in range(path_ct[tau]-1):
            if s >= len(ww):
                waysInBox[v][s+1][0] = 4090
            else:
                waysInBox[v][s+1][:len(ww[s])] = ww[s]
        v+=1
    for y in range(1, tau+1):
        ww = fetch_paths_if_match(tau, y)
        waysInBox[v][0][:2] = [tau, y]
        for s in range(path_ct[tau]-1):
            if s >= len(ww):
                waysInBox[v][s+1][0] = 4090
            else:
                waysInBox[v][s+1][:len(ww[s])] = ww[s]
        v+=1
    return waysInBox