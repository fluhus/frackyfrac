package main

import (
	"fmt"
	"iter"
	"math"
	"os"
	"runtime"
	"sort"

	"github.com/fluhus/biostuff/formats/newick"
	"github.com/fluhus/frackyfrac/common"
	"github.com/fluhus/gostuff/ppln"
)

const (
	// Check for abundance only if a node is a leaf.
	flatNodeOptimization = true

	// Print debug information, for development.
	debugPrints = false
)

func init() {
	if flatNodeOptimization && debugPrints {
		fmt.Fprintln(os.Stderr, "*** Flat node optimization ***")
	}
}

// Converts an abundance map to a list of flat nodes. Returns the sum of
// abundances under the given tree.
func abundanceToFlatNodes(abnd map[string]float64, tree *newick.Node,
	enum map[*newick.Node]int, result *[]flatNode) float64 {
	sum := 0.0
	for _, c := range tree.Children {
		sum += abundanceToFlatNodes(abnd, c, enum, result)
	}
	if flatNodeOptimization {
		if len(tree.Children) == 0 {
			if a := abnd[tree.Name]; a > 0 {
				sum += a
			}
		}
	} else {
		if a := abnd[tree.Name]; a > 0 {
			sum += a
		}
	}
	if sum > 0 {
		*result = append(*result, flatNode{enum[tree], sum})
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

// Returns all the unique names in the tree.
func treeNames(tree *newick.Node) map[string]struct{} {
	m := map[string]struct{}{}
	for n := range tree.PreOrder() {
		m[n.Name] = struct{}{}
	}
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
) iter.Seq[float64] {
	sets := make([][]flatNode, 0, len(abnd))
	enum := enumerateNodes(tree)
	fmt.Fprintln(os.Stderr, "Converting abundances")
	ppln.Serial[map[string]float64, []flatNode](
		*nt,
		ppln.SliceInput(abnd),
		func(a map[string]float64, _, _ int) ([]flatNode, error) {
			var set []flatNode
			abundanceToFlatNodes(a, tree, enum, &set)
			if !*nnorm {
				normalizeFlatNodes(set)
			}
			return set, nil
		},
		func(a []flatNode) error {
			sets = append(sets, a)
			return nil
		})
	treeDists := make([]float64, len(enum))
	for k, v := range enum {
		treeDists[v] = k.Distance
	}
	runtime.GC()
	fmt.Fprintln(os.Stderr, "Calculating distances")
	return unifracDists(sets, treeDists, weighted)
}

// Assigns an arbitrary unique number to each node in the tree.
func enumerateNodes(tree *newick.Node) map[*newick.Node]int {
	m := map[*newick.Node]int{}
	for n := range tree.PreOrder() {
		m[n] = len(m)
	}
	return m
}

// Represents a node in a tree. Used for comparisons using slices rather than
// tree objects.
type flatNode struct {
	id   int     // Node unique ID.
	abnd float64 // Sum of abundances under this node.
}

// Returns unweighted UniFrac between the two samples, not divided by the
// tree's sum.
func unifracDistUnweighted(a, b []flatNode, treeDists []float64) float64 {
	result := 0.0
	common := 0.0
	i, j := 0, 0
	for i < len(a) && j < len(b) {
		if a[i].id < b[j].id {
			result += treeDists[a[i].id]
			i++
			continue
		}
		if a[i].id > b[j].id {
			result += treeDists[b[j].id]
			j++
			continue
		}
		common += treeDists[a[i].id]
		i++
		j++
	}
	for _, x := range a[i:] {
		result += treeDists[x.id]
	}
	for _, x := range b[j:] {
		result += treeDists[x.id]
	}
	result /= (result + common)
	return result
}

// Returns weighted UniFrac between the two samples.
func unifracDistWeighted(a, b []flatNode, treeDists []float64) float64 {
	numer := 0.0
	denom := 0.0
	i, j := 0, 0
	for i < len(a) && j < len(b) {
		if a[i].id < b[j].id {
			numer += treeDists[a[i].id] * a[i].abnd
			denom += treeDists[a[i].id] * a[i].abnd
			i++
			continue
		}
		if a[i].id > b[j].id {
			numer += treeDists[b[j].id] * b[j].abnd
			denom += treeDists[b[j].id] * b[j].abnd
			j++
			continue
		}
		numer += treeDists[a[i].id] * math.Abs(a[i].abnd-b[j].abnd)
		denom += treeDists[a[i].id] * (a[i].abnd + b[j].abnd)
		i++
		j++
	}
	for _, x := range a[i:] {
		numer += treeDists[x.id] * x.abnd
		denom += treeDists[x.id] * x.abnd
	}
	for _, x := range b[j:] {
		numer += treeDists[x.id] * x.abnd
		denom += treeDists[x.id] * x.abnd
	}
	return numer / denom
}

// Returns the UniFrac distances between the given samples, in flat pyramid
// order.
func unifracDists(nodes [][]flatNode, treeDists []float64, weighted bool,
) iter.Seq[float64] {
	return func(yield func(float64) bool) {
		ppln.Serial(*nt,
			common.IterPairs(nodes),
			func(a [2][]flatNode, _, _ int) (float64, error) {
				aa := a
				if weighted {
					return unifracDistWeighted(aa[0], aa[1], treeDists), nil
				} else {
					return unifracDistUnweighted(aa[0], aa[1], treeDists), nil
				}
			}, func(a float64) error {
				if !yield(a) {
					return fmt.Errorf("")
				}
				return nil
			})
	}
}
