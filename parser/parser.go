// Package parser parses frcfrc's input tables.
package parser

import (
	"bufio"
	"fmt"
	"io"
	"iter"
	"math"
	"regexp"
	"strconv"

	"github.com/fluhus/gostuff/ppln/v2"
)

// Splits input rows into individual values.
var splitter = regexp.MustCompile(`\S+`)

// ParseAbundance parses the input abundance table. Returns a map for each row
// where keys are species and values are abundances.
func ParseAbundance(r io.Reader, ngoroutines int,
	f func(map[string]float64)) error {
	var names []string
	err := ppln.Serial(
		ngoroutines,
		func(yield func(string, error) bool) {
			for row, err := range iterRows(r) {
				if err != nil {
					yield("", err)
					return
				}
				if names == nil {
					parts := splitter.FindAllString(row, -1)
					if len(parts) == 0 {
						yield("", fmt.Errorf("row #1 has 0 values"))
						return
					}
					names = parts
					continue
				}
				if !yield(row, err) {
					return
				}
			}
		},
		func(a string, _, _ int) (map[string]float64, error) {
			return parseRow(a, names)
		},
		func(m map[string]float64) error {
			f(m)
			return nil
		})
	if err != nil {
		return err
	}
	return nil
}

func parseRow(row string, names []string) (map[string]float64, error) {
	parts := splitter.FindAllString(row, -1)
	if len(parts) != len(names) {
		return nil, fmt.Errorf("has %d values, expected %d",
			len(parts), len(names))
	}
	m := map[string]float64{}
	for i := range parts {
		f, err := strconv.ParseFloat(parts[i], 64)
		if err != nil {
			return nil, fmt.Errorf("value #%d: %v", i+1, err)
		}
		if math.IsNaN(f) || math.IsInf(f, 0) || f < 0 {
			return nil, fmt.Errorf("value #%d: bad value: %f",
				i+1, f)
		}
		if f == 0 {
			continue
		}
		m[names[i]] = f
	}
	return m, nil
}

// ParseSparseAbundance parses the input sparse abundance table. Returns a map
// for each row where keys are species and values are abundances.
func ParseSparseAbundance(r io.Reader, ngoroutines int,
	f func(map[string]float64)) error {
	err := ppln.Serial(ngoroutines,
		iterRows(r),
		func(a string, _, _ int) (map[string]float64, error) {
			return parseSparseRow(a)
		},
		func(m map[string]float64) error {
			f(m)
			return nil
		})
	if err != nil {
		return err
	}
	return nil
}

func parseSparseRow(row string) (map[string]float64, error) {
	parts := splitter.FindAllString(row, -1)
	m := make(map[string]float64, len(parts)*11/10)
	for i := range parts {
		species, val, err := splitSparse(parts[i])
		if err != nil {
			return nil, fmt.Errorf("value #%d: %v", i+1, err)
		}
		if species == "" {
			return nil, fmt.Errorf("value #%d: empty species name", i+1)
		}
		f, err := strconv.ParseFloat(val, 64)
		if err != nil {
			return nil, fmt.Errorf("value #%d: %v", i+1, err)
		}
		if math.IsNaN(f) || math.IsInf(f, 0) || f < 0 {
			return nil, fmt.Errorf("value #%d: bad value: %f", i+1, f)
		}
		if f == 0 {
			return nil, fmt.Errorf(
				"value #%d: zeros are not allowed in sparse format", i+1)
		}
		m[species] = f
	}
	return m, nil
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

func iterRows(r io.Reader) iter.Seq2[string, error] {
	return func(yield func(string, error) bool) {
		sc := bufio.NewScanner(r)
		sc.Buffer(nil, 1<<25)
		for sc.Scan() {
			if !yield(sc.Text(), nil) {
				return
			}
		}
		if err := sc.Err(); err != nil {
			yield("", err)
		}
	}
}
