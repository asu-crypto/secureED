#!/usr/bin/env python3

import sys
from dataProcessing.data_loader import fetch_data_idash
from dataProcessing.edit_dist import edit_distance_D

# Running for two genomes and fixed length:
if __name__ == "__main__":
    p1, p2 = int(sys.argv[-3]), int(sys.argv[-2])
    ln = int(sys.argv[-1])
    data = fetch_data_idash()
    if len(data) < ln:
        for i in range(len(data)):
            data[i] = data[i] + data[i]
    a, b = data[p1][:ln], data[p2][:ln]
    v, _ = edit_distance_D(a, b)
    print("RES:", v)

# Running for all genome pairs:
# if __name__ == "__main__":
#     data = fetch_data_idash()
#     for i in range(0, len(data)):
#         for j in range(i+1, len(data)):
#             if i == j:
#                 continue
#             v, _ = edit_distance_D(data[i], data[j])
#             print("ED", i, j, v)
