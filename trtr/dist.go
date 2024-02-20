package main

import (
	"fmt"
	"math"
	"os"

	"github.com/fluhus/biostuff/formats/newick"
	"github.com/fluhus/frackyfrac/common"
	"github.com/fluhus/gostuff/clustering"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/minhash"
	"github.com/fluhus/gostuff/ppln/v2"
	"golang.org/x/exp/maps"
)

const (
	// In development. If true, prints out the entropy of the distances.
	// Considering whether this information is valueable.
	reportEntropy = false
)

// Returns the entropy of the given distances.
func entropy(x []float64) float64 {
	counts := map[float64]int{}
	for _, xx := range x {
		counts[math.Round(xx*100)]++
	}
	return gnum.Entropy(maps.Values(counts))
}

// A tree node with depth instead of length.
type deepNode struct {
	name     string
	depth    float64
	children []*deepNode
}

// Creates a tree from the given sketches with names as the leaf names.
func makeTree(sketches []*minhash.MinHash[uint64], names []string) *newick.Node {
	if len(sketches) != len(names) {
		panic(fmt.Sprintf("mismatching lengths: %d, %d",
			len(sketches), len(names)))
	}
	var distances []float64
	ppln.Serial[[2]*minhash.MinHash[uint64], float64](
		*nt,
		common.IterPairs(sketches),
		func(a [2]*minhash.MinHash[uint64], i, g int) (float64, error) {
			return jaccardToMash(a[0].Jaccard(a[1])), nil
		},
		func(a float64) error {
			distances = append(distances, a)
			return nil
		},
	)
	hcl := clustering.Agglo(len(sketches), clustering.AggloAverage,
		func(i, j int) float64 {
			if i == j {
				return 0
			}
			return distances[ijToN(i, j)]
		})
	fmt.Printf("Ditances: [%.2f,%.2f] mean=%.2f+-%.2f\n",
		gnum.Min(distances), gnum.Max(distances), gnum.Mean(distances),
		gnum.Std(distances))
	if reportEntropy {
		fmt.Fprintf(os.Stderr, "Entropy=%.2f\n", entropy(distances))
	}
	var nodes []*deepNode
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
	}
	return nodes[len(nodes)-1].toNewickNode()
}

// Converts coordinates in the distance matrix to an index in the linear vector.
func ijToN(i, j int) int {
	if i == j {
		panic(fmt.Sprintf("i=j is not allowed (i=j=%v)", i))
	}
	if i < j {
		i, j = j, i
	}
	return i*(i-1)/2 + j
}

// Converts a Jaccard similarity score to Mash distance.
func jaccardToMash(jac float64) float64 {
	if jac == 0 {
		return 1
	}
	return min(-math.Log(2*jac/(1+jac))/float64(*k), 1)
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
