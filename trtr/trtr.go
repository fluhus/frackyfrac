// Command trtr creates phylogenetic trees from fasta files.
package main

import (
	"crypto/md5"
	"encoding/base32"
	"flag"
	"fmt"
	"iter"
	"maps"
	"os"
	"path/filepath"
	"regexp"
	"slices"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/sequtil"
	"github.com/fluhus/frackyfrac/common"
	"github.com/fluhus/gostuff/flagx"
	"github.com/fluhus/gostuff/iterx"
	"github.com/fluhus/gostuff/ppln"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
)

var (
	k    = flag.Uint("k", 0, "K-mer length, required")
	n    = flag.Uint("n", 10000, "Sketch length")
	fout = flag.String("o", "", "Path to output file (default stdout)")
	keep = flag.Bool("keep-temp", false, "Do not remove temporary files")
	nt   = flag.Int("t", 1, "Number of threads")
	gre  = flagx.Regexp("g", regexp.MustCompile(".*"), "Group genomes by this regex")
)

func main() {
	common.ExitIfError(parseArgs())

	files := expandFiles()
	if len(files) == 0 {
		common.ExitIfError(fmt.Errorf("no input files"))
	}
	fmt.Fprintln(os.Stderr, "Sketching", len(files), "files")

	tmp, err := os.MkdirTemp(".", "tmp.trtr.")
	common.ExitIfError(err)
	if !*keep {
		defer os.RemoveAll(tmp)
	}
	fmt.Fprintln(os.Stderr, "Temp dir:", tmp)

	pt := ptimer.NewMessage("{} files sketched")
	var sketchFiles []string
	var names []string
	ppln.Serial[[]*fasta.Fasta, [2]string](
		*nt,
		iterInputFiles(files),
		func(fas []*fasta.Fasta, i, g int) ([2]string, error) {
			name := string((*gre).Find(fas[0].Name))
			fout := filepath.Join(tmp, "sketch_"+strhash(name)+".json.gz")
			return [2]string{fout, name}, sketchFile(fas, fout)
		},
		func(a [2]string) error {
			sketchFiles = append(sketchFiles, a[0])
			names = append(names, a[1])
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
	tree := makeTree(sketches, names)
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
	result := sets.Set[string]{}
	for _, pat := range flag.Args() {
		files, _ := filepath.Glob(pat)
		result.Add(files...)
	}
	return slices.Sorted(maps.Keys(result))
}

// Iterates over grouped fasta entries, using the grouping regex.
func iterInputFiles(files []string) iter.Seq2[[]*fasta.Fasta, error] {
	return func(yield func([]*fasta.Fasta, error) bool) {
		for _, f := range files {
			for grp := range iterx.NewUnreader2(fasta.File(f)).GroupBy(func(oldT *fasta.Fasta, oldS error, newT *fasta.Fasta, newS error) bool {
				return slices.Equal((*gre).Find(oldT.Name), (*gre).Find(newT.Name))
			}) {
				if !yield(iterx.CollectErr(grp)) {
					return
				}
			}
		}
	}
}

func iterSeqs(fas []*fasta.Fasta) iter.Seq[[]byte] {
	return func(yield func([]byte) bool) {
		for _, fa := range fas {
			for seq := range sequtil.SubsequencesWith(fa.Sequence, "atcgATCG") {
				if !yield(seq) {
					return
				}
			}
		}
	}
}

// A convenience function for hashing strings.
func strhash(x string) string {
	h := md5.New()
	h.Write([]byte(x))
	return base32.StdEncoding.EncodeToString(h.Sum(nil))[:20]
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
