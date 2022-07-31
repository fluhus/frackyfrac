// Package parser parses frcfrc's input tables.
package parser

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"regexp"
	"strconv"

	"github.com/fluhus/gostuff/ppln"
)

// Splits input rows into individual values.
var splitter = regexp.MustCompile(`\S+`)

// ParseAbundance parses the input abundance table. Returns a map for each row
// where keys are species and values are abundances.
func ParseAbundance(r io.Reader, ngoroutines int,
	f func(map[string]float64)) error {
	var names []string
	err := ppln.Serial(ngoroutines,
		func(push func(string), stop func() bool) error {
			sc := bufio.NewScanner(r)
			sc.Buffer(nil, 1<<25)
			for sc.Scan() {
				if stop() {
					break
				}
				if names == nil {
					parts := splitter.FindAllString(sc.Text(), -1)
					if len(parts) == 0 {
						return fmt.Errorf("row #1 has 0 values")
					}
					names = parts
					continue
				}
				push(sc.Text())
			}
			return sc.Err()
		},
		func(a string, _, _ int) (parseResult, error) {
			r := parseRow(a, names)
			return r, r.err
		},
		func(a parseResult) error {
			f(a.m)
			return nil
		})
	if err != nil {
		return err
	}
	return nil
}

func parseRow(row string, names []string) parseResult {
	parts := splitter.FindAllString(row, -1)
	if len(parts) != len(names) {
		return parseResult{nil, fmt.Errorf("has %d values, expected %d",
			len(parts), len(names))}
	}
	m := map[string]float64{}
	for i := range parts {
		f, err := strconv.ParseFloat(parts[i], 64)
		if err != nil {
			return parseResult{nil, fmt.Errorf("value #%d: %v", i+1, err)}
		}
		if math.IsNaN(f) || math.IsInf(f, 0) || f < 0 {
			return parseResult{nil, fmt.Errorf("value #%d: bad value: %f",
				i+1, f)}
		}
		if f == 0 {
			continue
		}
		m[names[i]] = f
	}
	return parseResult{m, nil}
}

// ParseSparseAbundance parses the input sparse abundance table. Returns a map
// for each row where keys are species and values are abundances.
func ParseSparseAbundance(r io.Reader, ngoroutines int,
	f func(map[string]float64)) error {
	err := ppln.Serial(ngoroutines,
		func(push func(string), stop func() bool) error {
			sc := bufio.NewScanner(r)
			sc.Buffer(nil, 1<<25)
			for sc.Scan() {
				if stop() {
					break
				}
				push(sc.Text())
			}
			return sc.Err()
		},
		func(a string, _, _ int) (parseResult, error) {
			r := parseSparseRow(a)
			return r, r.err
		},
		func(a parseResult) error {
			f(a.m)
			return nil
		})
	if err != nil {
		return err
	}
	return nil
}

func parseSparseRow(row string) parseResult {
	parts := splitter.FindAllString(row, -1)
	m := make(map[string]float64, len(parts)*11/10)
	for i := range parts {
		species, val, err := splitSparse(parts[i])
		if err != nil {
			return parseResult{nil, fmt.Errorf("value #%d: %v", i+1, err)}
		}
		if species == "" {
			return parseResult{nil,
				fmt.Errorf("value #%d: empty species name", i+1)}
		}
		f, err := strconv.ParseFloat(val, 64)
		if err != nil {
			return parseResult{nil, fmt.Errorf("value #%d: %v", i+1, err)}
		}
		if math.IsNaN(f) || math.IsInf(f, 0) || f < 0 {
			return parseResult{nil,
				fmt.Errorf("value #%d: bad value: %f", i+1, f)}
		}
		if f == 0 {
			return parseResult{nil,
				fmt.Errorf("value #%d: zeros are not allowed in sparse format",
					i+1)}
		}
		m[species] = f
	}
	return parseResult{m, nil}
}

func splitSparse(s string) (string, string, error) {
	last := -1
	for i, c := range s {
		if c == ':' {
			last = i
		}
	}
	if last == -1 {
		return "", "", fmt.Errorf("no colon in %q", s)
	}
	return s[:last], s[last+1:], nil
}

// The result of parsing a single row.
type parseResult struct {
	m   map[string]float64
	err error
}
