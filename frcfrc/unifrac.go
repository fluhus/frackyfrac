package main

import (
	"fmt"
	"math"
	"sort"

	"github.com/fluhus/biostuff/formats/newick"
	"github.com/fluhus/gostuff/ppln"
)

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
	forEach func(float64) bool) {
	sets := make([][]flatNode, 0, len(abnd))
	enum := enumerateNodes(tree)
	ppln.Serial(*nt,
		func(push func(int), s ppln.Stopper) {
			for i := range abnd {
				push(i)
			}
		}, func(a int, _, _ int, s ppln.Stopper) []flatNode {
			var set []flatNode
			abundanceToFlatNodes(abnd[a], tree, enum, &set)
			normalizeFlatNodes(set)
			return set
		}, func(a []flatNode, s ppln.Stopper) {
			sets = append(sets, a)
		})
	unifracDists(sets, tree, weighted, forEach)
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
	result /= (result + common)
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
func unifracDists(nodes [][]flatNode, tree *newick.Node, weighted bool,
	forEach func(float64) bool) {
	type task struct {
		a, b []flatNode
	}

	ppln.Serial(*nt,
		func(push func(task), s ppln.Stopper) {
			for i, a := range nodes {
				for _, b := range nodes[:i] {
					push(task{a, b})
				}
			}
		},
		func(a task, _, _ int, s ppln.Stopper) float64 {
			aa := a
			if weighted {
				return unifracDistWeighted(aa.a, aa.b)
			} else {
				return unifracDistUnweighted(aa.a, aa.b)
			}
		}, func(a float64, s ppln.Stopper) {
			if !forEach(a) {
				s.Stop()
			}
		})
}
