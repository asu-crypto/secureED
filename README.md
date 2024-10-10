# Secure Edit Distance

We use the [MP-SPDZ library](https://github.com/data61/MP-SPDZ) to perform secure edit distance and propose numerous optimizations as described in our paper for efficient performance in a multiparty computation setting using both garbled circuits (GC) as well as secret sharing (SS).

## Setup

Keep the MP-SPDZ library and this folder in the same parent folder. 

```
├── ParentFolder/
│   ├── MP-SPDZ/
│   ├── secure-edit-distance/ 
```

For the initial setup, navigate to MP-SPDZ and run the following -

```
Scripts/tldr.sh
make -j 8 yao
make -j 8 real-bmr-party.x
```

## Usage

Our implementation has been tested/demonstrated with the iDash Workshop's genome dataset. The first step is to convert the genomes into a format appropriate for MP-SPDZ.

After placing `iDash.txt` in iDASH-Data/, the player data (genomes) can be created in the appropriate format by running either of the following commands in the secure-edit-distance folder (first for GC; second for SS)-

```
./data_generator.sh gc 0 1
./data_generator.sh ss 0 1
```

Here, 0 and 1 are the genome IDs. The IDs can be 0 to 50. The files `Player-Data/Input-P0-0` and `Player-Data/Input-P1-0` would be generated and copied to `MP-SPDZ/Player-Data/Input-P0-0` and `MP-SPDZ/Player-Data/Input-P1-0`. Following this, we can compile and run the program. Each nucleotide in the player data is two 1-bit values for SS and one 2-bit integer for GC.

For the garbled circuits approach, set the parameters `m`, `n`, and `tt` and then compile the code with - 

```
./compile.py -B 16 -G -l ../secure-edit-distance/runner_gc
```

Then in two different terminal windows, run the following -
```
./yao-party.x -p 0 runner_gc
./yao-party.x -p 1 runner_gc
```

For the GC approach in a malicious setting, use real-bmr-party.x instead of yao-party.x

For the secret shared approach, set the parameters `m`, `n`, `tau` and `tt` and then compile the code with -

```
./compile.py -R 64 ../secure-edit-distance/runner_ss
```

Then in two different terminal windows, run the following -
```
./semi2k-party.x -N 2 -p 0 runner_ss
./semi2k-party.x -N 2 -p 1 runner_ss
```

For the SS approach in a malicious setting, disable `program.use_split` in `runner_ss`, compile the program and then use SPDz2k with the following commands on different terminals-

```
./spdz2k-party.x -p 0 runner_ss
```

## Implementation Notes

For any error/queries related to MP-SPDZ, do take a look at their detailed [documentation](https://mp-spdz.readthedocs.io/en/latest/index.html).

### Evaluation of Thresholding

To evaluate the accuracy of the threshold determination phase, for the iDash dataset, the existing `iDash.txt` file would suffice. For the PGP dataset, section pairs need to be generated using the reference genome and the edits. Place `U654E5F.vcf` in PGP-Data/ and `hg19_chr1.fa` in PGP-Data/ref/. Then, run the following in the secure-edit-distance folder which would generate `mydash1.txt` and `mydash2.txt` in PGP-Data with 60 pairs of length 1040.

```
python3 dataProcessing/process_vcf.py 1040 60
```

Now, the error of the thresholding algorithm can be determined with the following commands for iDash and PGP-

```
python3 approximate_estimate_idash *seq_len* *loose_upper_bound* *segment_length*
python3 approximate_estimate_mydash *seq_len* *loose_upper_bound* *segment_length*
```

For testing different segment lengths, you could specify the actual values found with edit_distance_D and avoid recomputation. Refer approximate_estimate files for further information.

To evaluate in a secret-shared setting, you could also run `eval_t_opt_secure.sh` and `eval_t_opt_secure2.sh` on different terminals after compiling the appropriate runner file.

### Scaling SS for different tau values

If you decide to test out a different tau value than 1, 3, 4, do go ahead and run the following command in secure-edit-distance -

```
python3 dataProcessing\generate_ways.py 10
```

This generates all possible ways for tau = 2 to 10. Then use ukkonenFill rather than ukkonenfillNTau in runner_ss. These steps need to be done only once unless you need a higher tau than specified above.

## Authors

Biodesign Institute - Arizona State University