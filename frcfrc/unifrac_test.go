package main

import (
	"fmt"
	"reflect"
	"strings"
	"testing"

	"github.com/fluhus/biostuff/formats/newick"
)

func TestUniFrac_simple(t *testing.T) {
	treeText := "(s2:3,s1:1,s3:5);"
	tree, err := parseTree(treeText)
	if err != nil {
		t.Fatal("failed to parse tree:", err)
	}
	abnd := []map[string]float64{
		{"s1": 1, "s2": 1},
		{"s3": 1, "s2": 1},
	}
	want := []float64{6.0 / 9.0}
	var got []float64
	for f := range unifrac(abnd, tree, false) {
		got = append(got, f)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unifrac(%v, %q, false)=%v, want %v",
			abnd, treeText, got, want)
	}
}

func TestUniFrac_complex(t *testing.T) {
	treeText := "((s1:1,s2:3,s3:5):3,(s4:2,s5:2,s6:2):4,(s7:3,s8:2,s9:1):5);"
	tree, err := parseTree(treeText)
	if err != nil {
		t.Fatal("failed to parse tree:", err)
	}
	abnd := []map[string]float64{
		{"s1": 1, "s2": 1, "s5": 1, "s9": 1},
		{"s3": 1, "s4": 1, "s5": 1, "s6": 1},
		{"s7": 1, "s9": 1},
	}
	want := []float64{19.0 / 28.0, 16.0 / 22.0, 1.0}
	var got []float64
	for f := range unifrac(abnd, tree, false) {
		got = append(got, f)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unifrac(%v, %q, false)=%v, want %v",
			abnd, treeText, got, want)
	}
}

func TestUniFrac_weighted(t *testing.T) {
	treeText := "((s1:1,s2:3):2,(s3:2,s4:5):1);"
	tree, err := parseTree(treeText)
	if err != nil {
		t.Fatal("failed to parse tree:", err)
	}
	abnd := []map[string]float64{
		{"s1": 4, "s2": 1},
		{"s3": 3, "s2": 2},
	}
	want := []float64{22.0 / 36.0}
	var got []float64
	for f := range unifrac(abnd, tree, true) {
		got = append(got, f)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unifrac(%v, %q, true)=%v, want %v",
			abnd, treeText, got, want)
	}
}

func parseTree(s string) (*newick.Node, error) {
	for tr, err := range newick.Reader(strings.NewReader(s)) {
		return tr, err
	}
	return nil, fmt.Errorf("no tree provided")
}
