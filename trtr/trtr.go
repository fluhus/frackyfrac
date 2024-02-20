// Command trtr creates phylogenetic trees from fasta files.
package main

import (
	"crypto/md5"
	"encoding/base32"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/sequtil"
	"github.com/fluhus/frackyfrac/common"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/ppln"
	"github.com/fluhus/gostuff/ptimer"
	"golang.org/x/exp/maps"
)

var (
	k    = flag.Uint("k", 0, "K-mer length, required")
	n    = flag.Uint("n", 10000, "Sketch length")
	fout = flag.String("o", "", "Path to output file (default stdout)")
	keep = flag.Bool("keep-temp", false, "Do not remove temporary files")
	nt   = flag.Int("t", 1, "Number of threads")
)

func main() {
	common.ExitIfError(parseArgs())

	files := expandFiles()
	if len(files) == 0 {
		common.ExitIfError(fmt.Errorf("no input files"))
	}
	fmt.Fprintln(os.Stderr, "Sketching", len(files), "files")

	tmp, err := os.MkdirTemp("", "trtr-")
	common.ExitIfError(err)
	if !*keep {
		defer os.RemoveAll(tmp)
	}
	fmt.Fprintln(os.Stderr, "Temp dir:", tmp)

	pt := ptimer.NewMessage("{} files sketched")
	var sketchFiles []string
	ppln.Serial(
		*nt,
		func(push func(string), stop func() bool) error {
			for _, file := range files {
				push(file)
			}
			return nil
		},
		func(file string, i, g int) (string, error) {
			fout := filepath.Join(tmp, "sketch_"+strhash(file)+".json.gz")
			return fout, sketchFile(file, fout)
		},
		func(f string) error {
			sketchFiles = append(sketchFiles, f)
			pt.Inc()
			return nil
		})
	pt.Done()

	fmt.Fprintln(os.Stderr, "Loading sketches")
	pt = ptimer.New()
	sketches, err := loadSketches(sketchFiles)
	common.ExitIfError(err)
	pt.Done()

	fmt.Fprintln(os.Stderr, "Building tree")
	pt = ptimer.New()
	tree := makeTree(sketches, baseNames(files))
	pt.Done()

	treeText, _ := tree.MarshalText()
	if *fout == "" {
		fmt.Printf("%s\n", treeText)
	} else {
		common.ExitIfError(os.WriteFile(*fout, treeText, 0o644))
	}
}

// Parses and checks arguments.
func parseArgs() error {
	if len(os.Args) == 1 {
		usage()
		os.Exit(1)
	}
	flag.Usage = usage
	flag.Parse()

	if *k == 0 {
		return fmt.Errorf("please provide a kmer length with -k")
	}

	return nil
}

// Expands the argument glob patterns to files, removing repetitions.
func expandFiles() []string {
	result := map[string]struct{}{}
	for _, pat := range flag.Args() {
		files, _ := filepath.Glob(pat)
		for _, file := range files {
			result[file] = struct{}{}
		}
	}
	keys := maps.Keys(result)
	sort.Strings(keys)
	return keys
}

// Calls f for each canonical kmer in the given reader.
func iterKmers(r *aio.Reader, k int, f func([]byte)) error {
	fqr := fasta.NewReader(r)
	var err error
	var fq *fasta.Fasta
	for fq, err = fqr.Read(); err == nil; fq, err = fqr.Read() {
		seq := fq.Sequence
		sequtil.CanonicalSubsequences(seq, k, func(kmer []byte) bool {
			f(kmer)
			return true
		})
	}
	if err != io.EOF {
		return err
	}
	return nil
}

// A convenience function for hashing strings.
func strhash(x string) string {
	h := md5.New()
	h.Write([]byte(x))
	return base32.StdEncoding.EncodeToString(h.Sum(nil))[:20]
}

// Returns the basenames of the given files.
func baseNames(files []string) []string {
	result := make([]string, len(files))
	for i, f := range files {
		result[i] = filepath.Base(f)
	}
	return result
}

// Prints usage help message.
func usage() {
	fmt.Fprintln(flag.CommandLine.Output(),
		`TrickyTree creates a phylogenetic tree for use with UniFrac.

Usage:
trtr [PARAMS] species1.fa species2.fa species3.fa ...

File names may be glob patterns with '*', '?', and '[abc123]'.

Params:`)
	flag.PrintDefaults()
}
