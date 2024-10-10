#!/usr/bin/env python3

import sys
import json
from pathlib import Path

from old_path_finding import minimal_ways_in_box as ways_in_box_old
from new_path_finding import minimal_ways_in_box as ways_in_box_new
from new_path_finding import minimal_ways_in_match_box as ways_in_box_match

def generate_ways_old(lim=10):

    ways = {}
    for i in range(2, lim+1):
        ways[i] = ways_in_box_old(i)

    with open('Common-Data/ways-old.json', 'w') as f:
        f.write(json.dumps(ways))

def generate_ways_new(lim=10):

    ways = {}
    for i in range(2, lim+1):
        ways[i] = ways_in_box_new(i)

    with open('Common-Data/ways.json', 'w') as f:
        f.write(json.dumps(ways))

def generate_ways_match(lim=10):

    ways = {}
    for i in range(2, lim+1):
        ways[i] = ways_in_box_match(i)

    with open('Common-Data/waysMatch.json', 'w') as f:
        f.write(json.dumps(ways))

if __name__ == "__main__":
    Path("Common-Data/").mkdir(parents=True, exist_ok=True)
    lim = int(sys.argv[-1])
    generate_ways_new(lim)
    generate_ways_match(lim)