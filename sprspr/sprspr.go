// Converts dense format tables to sparse format tables.
package main

import (
	"fmt"
	"io"
	"os"

	"github.com/fluhus/frackyfrac/common"
	"github.com/fluhus/frackyfrac/parser"
	"github.com/fluhus/gostuff/ptimer"
)

func main() {
	fmt.Fprintln(os.Stderr, usageMessage)
	common.ExitIfError(toSparse(os.Stdin, os.Stdout))
}

func toSparse(r io.Reader, w io.Writer) error {
	t := ptimer.New()
	err := parser.ParseAbundance(r, 2, func(m map[string]float64) {
		first := true
		for k, v := range m {
			if first {
				first = false
			} else {
				fmt.Fprint(w, "\t")
			}
			fmt.Fprintf(w, "%s:%g", k, v)
		}
		fmt.Fprintln(w)
		t.Inc()
	})
	t.Done()
	return err
}

const usageMessage = `` +
	`SparseySparse converts dense format abundance tables to sparse format.

Usage:
sprspr < INPUT_FILE > OUTPUT_FILE

Reading standard input...`
