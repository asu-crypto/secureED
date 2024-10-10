#!/usr/bin/env python3

from Compiler.GC.types import result_conv
from Compiler.util import tree_reduce, if_else

def if_else_pair(condition, x, y, x2, y2):
        diff = x ^ y
        diff2 = x2 ^ y2
        this = result_conv(x, y)(condition * (diff) ^ y)
        that = result_conv(x2, y2)(condition * (diff2) ^ y2)
        return [this, that]

def tree_reduce_pair(function, seq, seq2):
    seq, seq2 = list(seq), list(seq2)
    assert len(seq) > 0
    n = len(seq)
    if n == 1:
        return [seq[0], seq2[0]]
    else:
        reduced = [[], []]
        for i in range(n//2):
            fa = function(seq[2*i], seq2[2*i], seq[2*i+1], seq2[2*i+1])
            reduced[0].append(fa[0])
            reduced[1].append(fa[1])
        return tree_reduce_pair(function, reduced[:][0] + seq[n//2*2:], reduced[:][1] + seq2[n//2*2:])

def min_pair(x, xp, y=None, yp=None):
    if y is None:
        return tree_reduce_pair(min_pair, x, xp)
    else:
        return if_else_pair((x < y), x, y, xp, yp)
    
def min_multidim(x, y=None):
    if y is None:
        return tree_reduce(min_multidim, x)
    else:
        return if_else(x[0] < y[0], x, y)