#!/usr/bin/env python3

def fetch_paths(x, y):

    ways = []
    ways.append([[x, 0], [], y])
    ways.append([[0, y], [], x])

    def possible_minimal_paths(initIdx, compareIndices, addValue, shift, dtrav, rtrav, i, j):
        if i == x and j == y:
            ways.append([initIdx, compareIndices, addValue])
            return
        dist = min(x-i, y-j)
        if j != 0 and dtrav == 0:
            if shift > 0:
                possible_minimal_paths(initIdx, compareIndices, addValue+1, shift-1, 0, 2, i+1, j)
            elif dist >= 3:
                possible_minimal_paths(initIdx, compareIndices, addValue+1, shift-1, 0, 2, i+1, j)
        if i != 0 and rtrav == 0:
            if shift < 0:
                possible_minimal_paths(initIdx, compareIndices, addValue+1, shift+1, 2, 0, i, j+1)
            elif dist >= 3:
                possible_minimal_paths(initIdx, compareIndices, addValue+1, shift+1, 2, 0, i, j+1)
        if i < x and j < y:
            possible_minimal_paths(initIdx, compareIndices+[[i+1, j+1]], addValue, shift, max(0, dtrav-1), max(rtrav-1, 0), i+1, j+1)

    for i in range(0, x):
        shift = (x-i)-y
        possible_minimal_paths([i, 0], [], 0, shift, 0, 0, i, 0)
    for j in range(1, y):
        shift = x-(y-j)
        possible_minimal_paths([0, j], [], 0, shift, 0, 0, 0, j)

    return ways

def fetch_paths_if_match(x, y):

    ways = []

    def possible_minimal_paths(initIdx, compareIndices, addValue, shift, i, j):
        if i == x and j == y:
            ways.append([initIdx, compareIndices, addValue])
            return
        if shift > 0 and j != 0:
            possible_minimal_paths(initIdx, compareIndices, addValue+1, shift-1, i+1, j)
        if shift < 0 and i != 0:
            possible_minimal_paths(initIdx, compareIndices, addValue+1, shift+1, i, j+1)
        if i < x and j < y:
            possible_minimal_paths(initIdx, compareIndices+[[i+1, j+1]], addValue, shift, i+1, j+1)

    if x == y:
        ways.append([[0, 0], [], 0])
        return ways
    ways.append([[0,0], [], abs(x-y)])
    if y < x:
        for i in range(1, x):
            shift = (x-i)-y
            possible_minimal_paths([i, 0], [], 0, shift, i, 0)
        ways.append([[x, 0], [], y])
    if x < y:
        for j in range(1, y):
            shift = x-(y-j)
            possible_minimal_paths([0, j], [], 0, shift, 0, j)
        ways.append([[0, y], [], x])

    return ways

def minimal_ways_in_box(tau, match=False):

    waysInBox = dict()
    for x in range(1, tau):
        if match:
            waysInBox[(x, tau)] = fetch_paths_if_match(x, tau)
        else:
            waysInBox[(x, tau)] = fetch_paths(x, tau)
    for y in range(1, tau+1):
        if match:
            waysInBox[(tau, y)] = fetch_paths_if_match(tau, y)
        else:
            waysInBox[(tau, y)] = fetch_paths(tau, y)

    return waysInBox