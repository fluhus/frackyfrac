# FrackyFrac: a standalone UniFrac calculator

Make UniFrac accessible to everyone.

## System Requirements

None. Just download the binaries and get going!

## Usage

### FrackyFrac (`frcfrc`)

Calculates pairwise UniFrac for an input abundance table.

**Input format**

A whitespace-separated table (spaces or tabs). Rows are samples, columns are
species. The first row should contain species names, which should match the
names in the phylogenetic tree. Currently species names with spaces are not
allowed. The next rows should contain floating-point values for each species.
The values don't need to sum up to 1.

Example with 3 species and 2 samples:

```
species_A species_B species_C
0.1       0         0.3
0         1.5       0.4
```

**Output format**

One distance per line. The order is
(1,2),(1,3)...(1,n),(2,3),(2,4)...(2,n)...(n-1,n).

Example with 4 samples, the parentheses are not included in the output:

```
0.2  (1,2)
0.1  (1,3)
0.9  (1,4)
1    (2,3)
0.3  (2,4)
0.6  (3,4)
```

### TrickyTree (`trtr`)

Creates a phylogenetic tree for use with any UniFrac.
Uses min-hashing on k-mers to calculate the distances between samples.

**Input format**

Input are fasta files, one per species, which may be gzip'ed. The given file
paths may include glob characters (`*`, `?`, `[]`) which will be expanded.

**Output format**

Output is in [Newick format][nwk]. Leaf names are the base names of the input
fasta files. For example, `/path/to/species_1.fa` will have a leaf named
`species_1.fa`.

[nwk]: https://en.wikipedia.org/wiki/Newick_format

**Picking *K*'s value**

*K* is the length of k-mers that are used for calculating the distances. A
low *K* brings the distances close to 0, and a high *K* brings the distances
close to 1. An optimal *K* creates a wide spectrum of distances, which makes the
tree informative.

TrickyTree prints out the entropy of distances at the end of its process.
It is recommended to use the *K* that maximizes that entropy.
