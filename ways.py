#!/usr/bin/env python3

import json

from Compiler.types import cint, Matrix, MultiArray

def readWays(tau, match=False):
    '''Read precomputed ways in dependency graph of box'''
    jsonString = ''
    filename = '../secure-edit-distance/Common-Data/waysMatch.json' if match else '../secure-edit-distance/Common-Data/ways.json'
    with open(filename) as f:
        jsonString = f.readline()
    ways = json.loads(jsonString)
    if str(tau) in ways:
        return ways[str(tau)]
    return []

def readWaysMultiProcMin(tau, match=False):
    tt = str(tau)
    jsonString = ''
    filename = '../secure-edit-distance/Common-Data/waysMatch.json' if match else '../secure-edit-distance/Common-Data/ways.json'
    with open(filename) as f:
        jsonString = f.readline()
    ways = json.loads(jsonString)
    if tt in ways:
        idxes = Matrix(len(ways[tt]), 2, cint)
        for x in range(len(ways[tt])):
            idxes[x][0] = ways[tt][x][0][0]
            idxes[x][1] = ways[tt][x][0][1]
        return idxes, ways[tt]
    return [], []

def readWaysMultiProcAll(tau, match=False):
    filename = '../secure-edit-distance/Common-Data/waysMatch.json' if match else '../secure-edit-distance/Common-Data/ways.json'
    tt = str(tau)
    jsonString = ''
    with open(filename) as f:
        jsonString = f.readline()
    ways = json.loads(jsonString)
    if tt in ways:
        ww = MultiArray([len(ways[tt]), len(ways[tt][0]), len(ways[tt][0][0])], cint)
        ww[:] = ways[tt]
        return ww
    return []