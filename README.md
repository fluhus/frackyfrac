# FrackyFrac: a standalone UniFrac calculator

Make UniFrac accessible to everyone.

## Download

See [releases][rls].

## System Requirements

None.

If your operating system and CPU architecture are not listeted in the
[releases][rls], please let me know by opening an issue and I will do my best
to add it.

[rls]: https://github.com/fluhus/frackyfrac/releases

## Usage

Creating a tree (optional):

```
trtr -k 10 -o my_genomes.tree species_1.fa species_2.fa species_3.fa
-- or --
trtr -k 10 -o my_genomes.tree "species_*.fa"
```

Calculating UniFrac:

```
frcfrc -t my_genomes.tree -i my_abundances.txt -o distances.txt
```

Consult the [wiki][wiki] for more details.

[wiki]: https://github.com/fluhus/frackyfrac/wiki
