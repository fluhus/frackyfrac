// Package common provides utilities that are common to the binaries in this
// repository.
package common

import (
	"fmt"
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
