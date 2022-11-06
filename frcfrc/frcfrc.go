// Command frcfrc calculates UniFrac distances on input samples.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"time"

	"github.com/fluhus/biostuff/formats/newick"
	"github.com/fluhus/frackyfrac/common"
	"github.com/fluhus/frackyfrac/parser"
	"github.com/fluhus/gostuff/aio"
)

var (
	fin    = flag.String("i", "", "Path to input file (default stdin)")
	fout   = flag.String("o", "", "Path to output file (default stdout)")
	ftree  = flag.String("t", "", "Path to tree file, required")
	wgt    = flag.Bool("w", false, "Use weighted UniFrac (default unweighted)")
	sparse = flag.Bool("s", false, "Input is in sparse format")
	nt     = flag.Int("p", 1, "Number of threads")
)

func main() {
	common.ExitIfError(parseArgs())

	t := time.Now()
	fmt.Fprintln(os.Stderr, "Reading tree")
	tree, err := readTree()
	common.ExitIfError(err)

	fmt.Fprintln(os.Stderr, "Loading abundances")
	r, err := openInput()
	common.ExitIfError(err)
	var abnd []map[string]float64
	if *sparse {
		err = parser.ParseSparseAbundance(r, *nt, func(m map[string]float64) {
			abnd = append(abnd, m)
		})
	} else {
		err = parser.ParseAbundance(r, *nt, func(m map[string]float64) {
			abnd = append(abnd, m)
		})
	}
	common.ExitIfError(err)

	fmt.Fprintln(os.Stderr, "Validating")
	common.ExitIfError(validateSpecies(abnd, tree))

	fmt.Fprintln(os.Stderr, "Calculating distances")
	w, err := openOutput()
	common.ExitIfError(err)
	unifrac(abnd, tree, *wgt, func(f float64) error {
		_, err := fmt.Fprintln(w, f)
		return err
	})
	w.Close()
	fmt.Fprintln(os.Stderr, "Took", time.Since(t))
	fmt.Fprintln(os.Stderr, "Done")
}

// Parses and validates arguments.
func parseArgs() error {
	if len(os.Args) == 1 {
		usage()
		os.Exit(0)
	}

	flag.Usage = usage
	flag.Parse()
	if *ftree == "" {
		return fmt.Errorf("please provide a tree file with -t")
	}
	if *nt < 1 {
		return fmt.Errorf("bad number of threads: %d", *nt)
	}
	return nil
}

// Opens the output file, or stdout.
func openInput() (io.ReadCloser, error) {
	if *fin != "" {
		return aio.Open(*fin)
	} else {
		return os.Stdin, nil
	}
}

// Opens the input file, or stdin.
func openOutput() (io.WriteCloser, error) {
	if *fout != "" {
		return aio.Create(*fout)
	} else {
		return os.Stdout, nil
	}
}

// Reads the tree from the path in the argument.
func readTree() (*newick.Node, error) {
	f, err := aio.Open(*ftree)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return newick.NewReader(f).Read()
}

// Prints usage help message.
func usage() {
	fmt.Fprintln(flag.CommandLine.Output(),
		`FrackyFrac calculates UniFrac on the given abundance table.
Outputs one distance per line in the order (1,2),(1,3)...(1,n),(2,3)...(2,n).

Params:`)
	flag.PrintDefaults()
}
