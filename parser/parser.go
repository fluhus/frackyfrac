// Package parser parses frcfrc's input tables.
package parser

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"regexp"
	"strconv"

	"github.com/fluhus/frackyfrac/ppln"
)

// Splits input rows into individual values.
var splitter = regexp.MustCompile(`\S+`)

// ParseAbundance parses the input abundance table. Returns a map for each row
// where keys are species and values are abundances.
func ParseAbundance(r io.Reader, ngoroutines int) ([]map[string]float64, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(nil, 1<<25)
	var names []string
	var result []map[string]float64
	var err error
	ppln.Serial(ngoroutines,
		func(push func(interface{}), s ppln.Stopper) {
			for sc.Scan() {
				if s.Stopped() {
					break
				}
				if names == nil {
					parts := splitter.FindAllString(sc.Text(), -1)
					if len(parts) == 0 {
						err = fmt.Errorf("row #1 has 0 values")
					}
					names = parts
					continue
				}
				push(sc.Text())
			}
		},
		func(a interface{}, s ppln.Stopper) interface{} {
			return parseRow(a.(string), names)
		},
		func(a interface{}, s ppln.Stopper) {
			aa := a.(parseResult)
			if aa.err != nil && err == nil { // First error.
				s.Stop()
				err = aa.err
			}
			result = append(result, aa.m)
		})
	if err != nil {
		return nil, err
	}
	if sc.Err() != nil {
		return nil, sc.Err()
	}
	return result, nil
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
func ParseSparseAbundance(r io.Reader, ngoroutines int) ([]map[string]float64, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(nil, 1<<25)
	var result []map[string]float64
	var err error
	ppln.Serial(ngoroutines,
		func(push func(interface{}), s ppln.Stopper) {
			for sc.Scan() {
				if s.Stopped() {
					break
				}
				push(sc.Text())
			}
		},
		func(a interface{}, s ppln.Stopper) interface{} {
			return parseSparseRow(a.(string))
		},
		func(a interface{}, s ppln.Stopper) {
			aa := a.(parseResult)
			if aa.err != nil && err == nil { // Take first error.
				s.Stop()
				err = aa.err
			}
			result = append(result, aa.m)
		})
	if err != nil {
		return nil, err
	}
	if sc.Err() != nil {
		return nil, sc.Err()
	}
	return result, nil
}

func parseSparseRow(row string) parseResult {
	parts := splitter.FindAllString(row, -1)
	m := map[string]float64{}
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