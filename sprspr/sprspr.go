// Converts dense format tables to sparse format tables.
package main

import (
	"fmt"
	"os"

	"github.com/fluhus/frackyfrac/parser"
	"github.com/fluhus/gostuff/ptimer"
)

func main() {
	fmt.Fprintln(os.Stderr, usageMessage)
	t := ptimer.New()
	err := parser.ParseAbundance(os.Stdin, 2, func(m map[string]float64) {
		first := true
		for k, v := range m {
			if first {
				first = false
			} else {
				fmt.Print("\t")
			}
			fmt.Printf("%s:%g", k, v)
		}
		fmt.Println()
		t.Inc()
	})
	t.Done()
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR:", err)
		os.Exit(2)
	}
}

const usageMessage = `` +
	`SparseySparse converts dense format abundance tables to sparse format.

Usage:
sprspr < INPUT_FILE > OUTPUT_FILE

Reading standard input...`
