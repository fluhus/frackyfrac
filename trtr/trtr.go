// Command trtr creates phylogenetic trees from fasta files.
package main

import (
	"bytes"
	"crypto/md5"
	"encoding/base32"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/sequtil"
	"github.com/fluhus/frackyfrac/common"
	"github.com/fluhus/gostuff/aio"
	"golang.org/x/exp/maps"
)

var (
	k    = flag.Uint("k", 0, "K-mer length, required")
	n    = flag.Uint("n", 10000, "Sketch length")
	fout = flag.String("o", "", "Path to output file (default stdout)")
	keep = flag.Bool("keep-temp", false, "Do not remove temporary files")
)

func main() {
	// ezpprof.Start("/tmp/amitmit/profile")
	// defer ezpprof.Stop()
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

	tim := common.NewTimerMessasge("* files sketched")
	var sketchFiles []string
	for _, file := range files {
		f := filepath.Join(tmp, "sketch_"+strhash(file)+".json.gz")
		common.ExitIfError(sketchFile(file, f))
		sketchFiles = append(sketchFiles, f)
		tim.Inc()
	}
	tim.Done()

	fmt.Fprintln(os.Stderr, "Loading sketches")
	sketches, err := loadSketches(sketchFiles)
	common.ExitIfError(err)

	fmt.Fprintln(os.Stderr, "Building tree")
	tim = common.NewTimerMessasge("tree building")
	tree := makeTree(sketches, baseNames(files))
	tim.Done()

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
		os.Exit(0)
	}
	flag.Usage = usage
	flag.Parse()

	if *k == 0 {
		return fmt.Errorf("please provide a kmer length with -k")
	}

	return nil
}

// Expands the given glob patterns to files, removing repetitions.
func expandFiles() []string {
	result := map[string]struct{}{}
	for _, pat := range flag.Args() {
		files, _ := filepath.Glob(pat)
		for _, file := range files {
			result[file] = struct{}{}
		}
	}
	return maps.Keys(result)
}

// Calls f for each canonical kmer in the given reader.
func iterKmers(r *aio.Reader, k int, f func([]byte)) error {
	fqr := fasta.NewReader(r)
	var err error
	var fq *fasta.Fasta
	var rc []byte
	for fq, err = fqr.Read(); err == nil; fq, err = fqr.Read() {
		seq := fq.Sequence
		rc = sequtil.ReverseComplement(rc[:0], seq)
		nk := len(seq) - k + 1
		for i := 0; i < nk; i++ {
			kmer := seq[i : i+k]
			kmerRC := rc[len(rc)-i-k : len(rc)-i]
			if bytes.Compare(kmer, kmerRC) == 1 {
				kmer = kmerRC
			}
			f(kmer)
		}
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
