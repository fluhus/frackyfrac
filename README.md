# FrackyFrac: a standalone UniFrac calculator

Make UniFrac accessible to everyone.

## Download

Go to [releases][rls] and download the relevant zip file.

Unzip the file and run the relevant executable from the terminal
(see instructions below).

## System Requirements

None.

If your operating system and CPU architecture are not listeted in the
[releases][rls], please let me know by opening an issue and I will do my best
to add it.

[rls]: https://github.com/fluhus/frackyfrac/releases

## Usage

Calculating UniFrac:

```
frcfrc -t my_genomes.tree -i my_abundances.tsv -o distances.txt
```

Creating a tree (optional):

```
trtr -k 21 -o my_genomes.tree species_1.fa species_2.fa species_3.fa
-- or --
trtr -k 21 -o my_genomes.tree "species_*.fa"
```

Consult the [wiki][wiki] for more details.

## Output

The output is the lower triangle of the distance matrix.

Example:

```python
import numpy as np

def load_frcfrc_distances(file:str, num_samples:int):
	distances = [float(x) for x in open(file)]
	mat = np.zeros([num_samples,num_samples])
	pos = np.tril_indices(num_samples,-1)
	mat[pos] = distances  # Lower triangle.
	mat[pos[::-1]] = distances  # Upper triangle.
	return mat
```

[wiki]: https://github.com/fluhus/frackyfrac/wiki

## Testing (for developers & reviewers)

See the `testdata` directory for tests of this implementation.
