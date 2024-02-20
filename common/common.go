// Package common provides utilities that are common to the binaries in this
// repository.
package common

import (
	"fmt"
	"iter"
	"os"
)

// ExitIfError prints the given error and exits if err is not nil.
// Should generally be used from main functions.
func ExitIfError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR:", err)
		os.Exit(2)
	}
}

// IterPairs returns an iterator over pairs of elements in s.
func IterPairs[T any](s []T) iter.Seq2[[2]T, error] {
	return func(yield func([2]T, error) bool) {
		for i := range s {
			for j := range i {
				if !yield([2]T{s[i], s[j]}, nil) {
					return
				}
			}
		}
	}
}
