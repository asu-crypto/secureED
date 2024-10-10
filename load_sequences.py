#!/usr/bin/env python3

from Compiler.types import sint, Array, Matrix
from Compiler.GC.types import sbitint, sbits

def tau_adjusted(v: int, tau: int):
    if v % tau == 0:
        v += 1
    elif v % tau > 1:
        v += (tau - v%tau) + 1
    return v

def load_sequences(m: int, n: int, tau: int):

    m_adj, n_adj = m, n
    if tau > 1:
        m_adj, n_adj = tau_adjusted(m, tau), tau_adjusted(n, tau)
    
    a = Array(m_adj, sint)
    b = Array(n_adj, sint)

    for i in range(m):
        a[i] = sint(sint.get_input_from(0))
    for i in range(n):
        b[i] = sint(sint.get_input_from(1))

    return a, b

def load_sequences_gc(m: int, n: int, tau: int):

    m_adj, n_adj = m, n
    if tau > 1:
        m_adj, n_adj = tau_adjusted(m, tau), tau_adjusted(n, tau)

    si4 = sbitint.get_type(4) # For garbled circuits, use sbitint of length 4
    a = [si4(0) for _ in range(m_adj)] # No advantage of using runtime `Array` type
    b = [si4(0) for _ in range(n_adj)]

    for i in range(m):
        a[i] = si4(sbitint.get_input_from(0))
    for i in range(n):
        b[i] = si4(sbitint.get_input_from(1))

    return a, b

def load_sequences_shuffled(m: int, n: int, tau: int):

    m_adj = max(m, n)
    if tau > 1:
        m_adj = tau_adjusted(max(m, n), tau)
    
    temp = Matrix(m_adj, 2, sint)
    a = Array(m_adj, sint)
    b = Array(m_adj, sint)

    for i in range(m):
        temp[i][0] = sint.get_input_from(0)
    for i in range(m):
        temp[i][1] = sint.get_input_from(1)

    temp.secure_shuffle()

    for i in range(m):
        a[i] = temp[i][0]
        b[i] = temp[i][1]

    return a, b

def load_sequences_bits(m: int, n: int, tau: int):

    sb = sbits.get_type(1)

    m_adj, n_adj = m, n
    if tau > 1:
        m_adj, n_adj = tau_adjusted(m, tau), tau_adjusted(n, tau)
    
    a = Matrix(m_adj, 2, sb)
    b = Matrix(n_adj, 2, sb)

    for i in range(m):
        a[i][0] = sb.get_input_from(0)
        a[i][1] = sb.get_input_from(0)
    for i in range(n):
        b[i][0] = sb.get_input_from(1)
        b[i][1] = sb.get_input_from(1)

    return a, b