#!/usr/bin/env python3

import math
import random
import sys
sys.path.insert(0, '../secure-edit-distance')

from dataProcessing.data_loader import fetch_data_mydash
from dataProcessing.edit_dist import edit_distance_D
from t_optimize_nonsecure import find_opt_nonsecure

def root_mean_square_error(actual, estimate):
    rmse = 0
    match_ct = 0
    for i in range(len(actual)):
        if actual[i] != 0: # Ignore if two genomes are same
            rmse += (((estimate[i]-actual[i]) / actual[i]) ** 2)
        else:
            match_ct+=1
    print("Matches:", match_ct)
    rmse = rmse / (len(actual)-match_ct)
    rmse = math.sqrt(rmse)
    return rmse

def avg_error(actual, estimate):
    # err = 0
    # match_ct = 0
    # for i in range(len(a)):
    #     if actual[i] != 0: # Ignore if two genomes are same
    #         err += ((actual[i]-estimate[i]) / actual[i])
    #     else:
    #         match_ct+=1
    # print("Matches:", match_ct)
    # return err / (len(actual)-match_ct)
    s1 = sum(actual)
    s2 = sum(estimate)
    return (s2-s1)/s1

def replicate_if_needed(a_seq, b_seq, ln):
    genome_lengths = [len(d) for d in b_seq]
    shortest_genome_len = min(genome_lengths)
    shortest_genome_idx = genome_lengths.index(shortest_genome_len)
    print("Shortest sequence length:", shortest_genome_len)
    while ln > len(b_seq[shortest_genome_idx]): # If length exceeds iDASH data length, replicate
       print("Replicating..")
       for i in range(len(a_seq)):
            a_seq[i] = a_seq[i] + a_seq[i]
            b_seq[i] = b_seq[i] + b_seq[i]
    return a_seq, b_seq


if __name__ == "__main__":

    ln = int(sys.argv[-3]) # Sequence length
    t, modval = int(sys.argv[-2]), int(sys.argv[-1]) # Modval is segment length
    actuals, estimates = [], []

    # edit_distance_D takes a lot of time to run and is independent of segment length.
    # If you've already found the actual avg. ED for a sequence length, fill it here to
    # avoid unnecessary reruns when you exxpiremeng with different segment lengths - 
    actuals_known = False # Make True 
    actuals = [] # Replace with known truth values

    a_seq, b_seq = fetch_data_mydash()
    a_seq, b_seq = replicate_if_needed(a_seq, b_seq, ln)
    genome_ct = len(a_seq) # Right now, 60 for PGP data; Can lower if you don't want to eval all possible pairs.
    
    for i in range(genome_ct):
        a, b = a_seq[i], b_seq[i]
        if not actuals_known:
            actual, _ = edit_distance_D(a, b)
            actuals.append(actual)
        estimate, _, _ = find_opt_nonsecure(a, b, t, modval)
        estimates.append(estimate)
        if (i+1) % 10 == 0:
            print("Completed:", i+1, "/", genome_ct)

    if not actuals_known:
        print(actuals)
    print()
    print(estimates)
    print("\nNo. of pairs:", len(actuals))
    print("Actual avg. {:.1f}".format(sum(actuals)/len(actuals)))
    print("Approx avg. {:.1f}".format(sum(estimates)/len(estimates)))
    print("RMSE: {:.4f}".format(root_mean_square_error(actuals, estimates)))
    print("ERR: {:.4f}".format(avg_error(actuals, estimates)))