package main

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"regexp"
	"sort"
	"strconv"

	"github.com/fluhus/biostuff/formats/newick"
	"github.com/fluhus/frackyfrac/ppln"
)

// TODO(amit): Clean this up.
const divideByUnion = true

// Splits input rows into individual values.
var splitter = regexp.MustCompile(`\S+`)

// Parses the input abundance table. Returns a map for each row where keys are
// species and values are abundances.
func parseAbundance(r io.Reader) ([]map[string]float64, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(nil, 1<<25)
	var names []string
	var result []map[string]float64
	for sc.Scan() {
		if names == nil {
			parts := splitter.FindAllString(sc.Text(), -1)
			if len(parts) == 0 {
				return nil, fmt.Errorf("row #1 has 0 values")
			}
			names = parts
			continue
		}
		m, err := parseRow(sc.Text(), names)
		if err != nil {
			return nil, fmt.Errorf("row #%d: %v", len(result)+2, err)
		}
		result = append(result, m)
	}
	if sc.Err() != nil {
		return nil, sc.Err()
	}
	return result, nil
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
			return nil, fmt.Errorf("value #%d: bad value: %f", i+1, f)
		}
		if f == 0 {
			continue
		}
		m[names[i]] = f
	}
	return m, nil
}

// Parses the input sparse abundance table. Returns a map for each row where
// keys are species and values are abundances.
func parseSparseAbundance(r io.Reader) ([]map[string]float64, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(nil, 1<<25)
	var result []map[string]float64
	type bla struct {
		m   map[string]float64
		err error
	}
	var err error
	ppln.Serial(*nt,
		func(c chan<- interface{}) {
			for sc.Scan() {
				c <- sc.Text()
			}
		},
		func(a interface{}) interface{} {
			m, err := parseSparseRow(a.(string))
			return bla{m, err}
		},
		func(a interface{}) {
			aa := a.(bla)
			if aa.err != nil {
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

func parseSparseRow(row string) (map[string]float64, error) {
	parts := splitter.FindAllString(row, -1)
	m := map[string]float64{}
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
			return nil, fmt.Errorf("value #%d: "+
				"zeros are not allowed in sparse format", i+1)
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

// Converts an abundance map to a list of flat nodes. Returns the sum of
// abundances under the given tree.
func abundanceToFlatNodes(abnd map[string]float64, tree *newick.Node,
	enum map[*newick.Node]int, result *[]flatNode) float64 {
	sum := 0.0
	for _, c := range tree.Children {
		sum += abundanceToFlatNodes(abnd, c, enum, result)
	}
	if a := abnd[tree.Name]; a > 0 {
		sum += a
	}
	if sum > 0 {
		*result = append(*result, flatNode{enum[tree], tree.Distance, sum})
	}
	return sum
}

// Sorts and divides abundances by their sum.
func normalizeFlatNodes(nodes []flatNode) {
	sort.Slice(nodes, func(i, j int) bool {
		return nodes[i].id < nodes[j].id
	})
	sum := 0.0
	for i := range nodes {
		sum += nodes[i].abnd
	}
	for i := range nodes {
		nodes[i].abnd /= sum
	}
}

// Returns the sum of distances in a tree and its subtrees.
func treeSum(tree *newick.Node) float64 {
	sum := 0.0
	tree.PreOrder(func(n *newick.Node) bool {
		sum += n.Distance
		return true
	})
	return sum
}

// Returns all the unique names in the tree.
func treeNames(tree *newick.Node) map[string]struct{} {
	m := map[string]struct{}{}
	tree.PreOrder(func(n *newick.Node) bool {
		m[n.Name] = struct{}{}
		return true
	})
	return m
}

// Validates that all species names in the given abundances list are in the
// tree.
func validateSpecies(abnd []map[string]float64, tree *newick.Node) error {
	species := treeNames(tree)
	for i, m := range abnd {
		for name, val := range m {
			if _, ok := species[name]; !ok {
				return fmt.Errorf(
					"sample #%d has value %v for species %q "+
						"which is not in the tree",
					i+1, val, name)
			}
		}
	}
	return nil
}

// Returns the unifrac distances between the given abundances, in flat pyramid
// order.
func unifrac(abnd []map[string]float64, tree *newick.Node, weighted bool,
) []float64 {
	sets := make([][]flatNode, 0, len(abnd))
	enum := enumerateNodes(tree)
	ppln.Serial(*nt,
		func(c chan<- interface{}) {
			for i := range abnd {
				c <- i
			}
		}, func(a interface{}) interface{} {
			var set []flatNode
			abundanceToFlatNodes(abnd[a.(int)], tree, enum, &set)
			normalizeFlatNodes(set)
			return set
		}, func(a interface{}) {
			sets = append(sets, a.([]flatNode))
		})
	return unifracDists(sets, tree, weighted)
}

// Assigns an arbitrary unique number to each node in the tree.
func enumerateNodes(tree *newick.Node) map[*newick.Node]int {
	m := map[*newick.Node]int{}
	tree.PreOrder(func(n *newick.Node) bool {
		m[n] = len(m)
		return true
	})
	return m
}

// Represents a node in a tree. Used for comparisons using slices rather than
// tree objects.
type flatNode struct {
	id   int     // Node unique ID.
	dist float64 // Distance from parent.
	abnd float64 // Sum of abundances under this node.
}

// Returns unweighted UniFrac between the two samples, not divided by the
// tree's sum.
func unifracDistUnweighted(a, b []flatNode) float64 {
	result := 0.0
	common := 0.0
	i, j := 0, 0
	for i < len(a) && j < len(b) {
		if a[i].id < b[j].id {
			result += a[i].dist
			i++
			continue
		}
		if a[i].id > b[j].id {
			result += b[j].dist
			j++
			continue
		}
		if a[i].dist != b[j].dist {
			panic(fmt.Sprintf("mismatching distances at (%d,%d): %f, %f",
				i, j, a[i].dist, b[j].dist))
		}
		common += a[i].dist
		i++
		j++
	}
	for _, x := range a[i:] {
		result += x.dist
	}
	for _, x := range b[j:] {
		result += x.dist
	}
	if divideByUnion {
		result /= (result + common)
	}
	return result
}

// Returns weighted UniFrac between the two samples.
func unifracDistWeighted(a, b []flatNode) float64 {
	numer := 0.0
	denom := 0.0
	i, j := 0, 0
	for i < len(a) && j < len(b) {
		if a[i].id < b[j].id {
			numer += a[i].dist * a[i].abnd
			denom += a[i].dist * a[i].abnd
			i++
			continue
		}
		if a[i].id > b[j].id {
			numer += b[j].dist * b[j].abnd
			denom += b[j].dist * b[j].abnd
			j++
			continue
		}
		if a[i].dist != b[j].dist {
			panic(fmt.Sprintf("mismatching distances at (%d,%d): %f, %f",
				i, j, a[i].dist, b[j].dist))
		}
		numer += a[i].dist * math.Abs(a[i].abnd-b[j].abnd)
		denom += a[i].dist * (a[i].abnd + b[j].abnd)
		i++
		j++
	}
	for _, x := range a[i:] {
		numer += x.dist * x.abnd
		denom += x.dist * x.abnd
	}
	for _, x := range b[j:] {
		numer += x.dist * x.abnd
		denom += x.dist * x.abnd
	}
	return numer / denom
}

// Returns the UniFrac distances between the given samples, in flat pyramid
// order.
func unifracDists(x [][]flatNode, tree *newick.Node, weighted bool) []float64 {
	sum := 1.0
	if !divideByUnion {
		sum = treeSum(tree)
	}

	type task struct {
		a, b []flatNode
	}
	result := make([]float64, 0, len(x)*(len(x)-1)/2)

	ppln.Serial(*nt,
		func(push chan<- interface{}) {
			for i, a := range x {
				for _, b := range x[i+1:] {
					push <- task{a, b}
				}
			}
		},
		func(a interface{}) interface{} {
			aa := a.(task)
			if weighted {
				return unifracDistWeighted(aa.a, aa.b)
			} else {
				return unifracDistUnweighted(aa.a, aa.b) / sum
			}
		}, func(a interface{}) {
			result = append(result, a.(float64))
		})

	return result
}
