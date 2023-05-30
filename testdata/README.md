# Test Data

Simple data to test this implementation and compare it with QIIME's.

## Running a self-test

Compares FrackyFrac's results against manually calculated results.
Make sure the binaries are in your PATH.

Linux:

```sh
sh run.sh  # no output means OK
```

Windows:

```cmd
run.bat
```

## Running a QIIME comparison

Compares QIIME's implementation (SSU) against the manually calculated results,
to verify that the two implementations agree.
Make sure QIIME is installed and is in your PATH.

Linux only:

```
sh run_ssu.sh
```

## Files in this directory

| Name | Description |
| - | - |
| *.dense | Dense abundance table for each test |
| *.sparse | Sparse abundance table for each test |
| *.tree | Species tree for each test |
| *.want | Expected result for each test |
| *.biom.tsv | BIOM table for each test |

The tests are `wtd` for weighted UniFrac and `uwtd` for unweighted.
