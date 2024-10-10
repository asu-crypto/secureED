#!/usr/bin/env python3

import sys
from pathlib import Path

from data_loader import fetch_data_idash, fetch_aligned_idash_only_sequences

def generate_player_data(p1, p2, double=False):

    data = fetch_data_idash()

    if p1 > len(data) or p2 > len(data):
        print("Not found! The dataset's range is 0 to", len(data)-1)
        return

    if double:
        for i in range(len(data)):
            data[i] = data[i] + data[i]

    mp = {'A': '0', 'C': '1', 'G': '2', 'T': '3', '_': '4'}
    s1 = [mp[d] for d in data[p1]]
    s2 = [mp[d] for d in data[p2]]

    s1 = ' '.join(s1)
    with open('Player-Data/Input-P0-0', 'w') as f:
        f.write(s1)

    s2 = ' '.join(s2)
    with open('Player-Data/Input-P1-0', 'w') as f:
        f.write(s2)

def generate_player_data_bitops(p1, p2, double=False):

    data = fetch_data_idash()

    if p1 > len(data) or p2 > len(data):
        print("Not found! The dataset's range is 0 to", len(data)-1)
        return

    if double:
        for i in range(len(data)):
            data[i] = data[i] + data[i]

    mp = {'A': '0 0', 'C': '0 1', 'G': '1 0', 'T': '1 1'}
    s1 = [mp[d] for d in data[p1]]
    s2 = [mp[d] for d in data[p2]]

    s1 = ' '.join(s1)
    with open('Player-Data/Input-P0-0', 'w') as f:
        f.write(s1)

    s2 = ' '.join(s2)
    with open('Player-Data/Input-P1-0', 'w') as f:
        f.write(s2)

if __name__ == "__main__":
    p1, p2 = int(sys.argv[-2]), int(sys.argv[-1])
    do_bitops = int(sys.argv[-3])
    Path("Player-Data/").mkdir(parents=True, exist_ok=True)
    if do_bitops:
        # Pass True if you want to double the size (Eg: If sequence length 4k is needed on iDash data)
        generate_player_data_bitops(p1, p2, double=False) 
    else:
        # Pass True if you want to double the size (Eg: If sequence length 4k is needed on iDash data)
        generate_player_data(p1, p2, double=False)
