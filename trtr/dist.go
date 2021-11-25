package main

import (
	"fmt"
	"math"
	"os"

	"github.com/fluhus/biostuff/formats/newick"
	"github.com/fluhus/gostuff/clustering"
)

// Returns the Jaccard dissimilarity between 2 sketches.
func dist(a, b []uint64) float64 {
	i, j := 0, 0
	common := 0
	for i < len(a) && j < len(b) {
		if a[i] < b[j] {
			i++
			continue
		}
		if a[i] > b[j] {
			j++
			continue
		}
		common++
		i++
		j++
	}
	return 1 - float64(common)/float64(len(a)+len(b)-common)
}

// Returns the entropy of dividing the given distances to buckets.
func entropy(x []float64, buckets int) float64 {
	counts := map[int]int{}
	fbuckets := float64(buckets)
	for _, xx := range x {
		if xx < 0 || xx > 1 {
			panic(fmt.Sprintf("value out of [0,1]: %f", xx))
		}
		counts[int(xx*fbuckets)]++
	}
	result := 0.0
	for _, count := range counts {
		p := float64(count) / float64(len(x))
		result -= p * math.Log2(p)
	}
	return result
}

// Integer square-root.
func isqrt(i int) int {
	return int(math.Round(math.Sqrt(float64(i))))
}

// A tree node with depth instead of length.
type deepNode struct {
	name     string
	depth    float64
	children []*deepNode
}

// Creates a tree from the given sketches with names as the leaf names.
func makeTree(sketches [][]uint64, names []string) *newick.Node {
	if len(sketches) != len(names) {
		panic(fmt.Sprintf("mismatching lengths: %d, %d",
			len(sketches), len(names)))
	}
	hcl := clustering.Agglo(len(sketches), clustering.AggloMax,
		func(i, j int) float64 {
			return dist(sketches[i], sketches[j])
		})
	var nodes []*deepNode
	var distances []float64
	for _, name := range names {
		nodes = append(nodes, &deepNode{name: name})
	}
	for i := 0; i < hcl.Len(); i++ {
		step := hcl.Step(i)
		node1, node2 := nodes[step.C1], nodes[step.C2]
		depth := step.D / 2 // Divided by the 2 children.
		// TODO(amit): Consider dividing new depth relative to each child node's
		// depth.
		parent := &deepNode{children: []*deepNode{node1, node2}, depth: depth}
		nodes[step.C2] = parent
		distances = append(distances, step.D)
	}
	fmt.Fprintln(os.Stderr, "Entropy:", entropy(distances, isqrt(len(distances))))
	return nodes[len(nodes)-1].toNewickNode()
}

// Converts a node with depth to a node with distance.
func (n *deepNode) toNewickNode() *newick.Node {
	node := &newick.Node{
		Name:     n.name,
		Children: nil,
	}
	if len(n.children) > 0 {
		node.Children = make([]*newick.Node, len(n.children))
	}
	for i, c := range n.children {
		node.Children[i] = c.toNewickNode()
		node.Children[i].Distance = n.depth - c.depth
	}
	return node
}
