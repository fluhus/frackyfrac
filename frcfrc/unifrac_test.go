package main

import (
	"reflect"
	"strings"
	"testing"

	"github.com/fluhus/biostuff/formats/newick"
)

func TestParseAbundance(t *testing.T) {
	input := "   aa  bbbb    \n1\t2\n 3  \t  4 \t\n"
	want := []map[string]float64{
		{"aa": 1, "bbbb": 2},
		{"aa": 3, "bbbb": 4},
	}
	got, err := parseAbundance(strings.NewReader(input))
	if err != nil {
		t.Fatalf("parseAbundance(%q) failed: %v", input, err)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("parseAbundance(%q)=%v, want %v", input, got, want)
	}
}

func TestParseAbundanceSparse(t *testing.T) {
	input := "a:11 b:222  \n  b:32 c:7\n\nd:1\tc:4\ta:10\n"
	want := []map[string]float64{
		{"a": 11, "b": 222},
		{"b": 32, "c": 7},
		{},
		{"d": 1, "c": 4, "a": 10},
	}
	got, err := parseSparseAbundance(strings.NewReader(input))
	if err != nil {
		t.Fatalf("parseAbundanceSparse(%q) failed: %v", input, err)
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("parseAbundanceSparse(%q)=%v, want %v", input, got, want)
	}
}
func TestSplitSparse(t *testing.T) {
	tests := []struct {
		in, want1, want2 string
	}{
		{"a:b", "a", "b"},
		{"c:d:e::f", "c:d:e:", "f"},
		{":", "", ""},
		{"a:", "a", ""},
		{":b", "", "b"},
	}
	for _, test := range tests {
		a, b, err := splitSparse(test.in)
		if err != nil {
			t.Errorf("splitSparse(%q) failed: %v", test.in, err)
			continue
		}
		if a != test.want1 || b != test.want2 {
			t.Errorf("splitSparse(%q)=%q,%q want %q,%q",
				test.in, a, b, test.want1, test.want2)
		}
	}
}

func TestSplitSparse_bad(t *testing.T) {
	tests := []string{"", "a", "aaa"}
	for _, test := range tests {
		a, b, err := splitSparse(test)
		if err == nil {
			t.Errorf("splitSparse(%q)=%q,%q want error",
				test, a, b)
		}
	}
}

func TestUniFrac_simple(t *testing.T) {
	treeText := "(s2:3,s1:1,s3:5);"
	tree, err := newick.NewReader(strings.NewReader(treeText)).Read()
	if err != nil {
		t.Fatal("failed to parse tree:", err)
	}
	abnd := []map[string]float64{
		{"s1": 1, "s2": 1},
		{"s3": 1, "s2": 1},
	}
	want := []float64{6.0 / 9.0}
	got := unifrac(abnd, tree, false)
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unifrac(%v, %q, false)=%v, want %v",
			abnd, treeText, got, want)
	}
}

func TestUniFrac_complex(t *testing.T) {
	treeText := "((s1:1,s2:3,s3:5):3,(s4:2,s5:2,s6:2):4,(s7:3,s8:2,s9:1):5);"
	tree, err := newick.NewReader(strings.NewReader(treeText)).Read()
	if err != nil {
		t.Fatal("failed to parse tree:", err)
	}
	abnd := []map[string]float64{
		{"s1": 1, "s2": 1, "s5": 1, "s9": 1},
		{"s3": 1, "s4": 1, "s5": 1, "s6": 1},
		{"s7": 1, "s9": 1},
	}
	want := []float64{19.0 / 28.0, 16.0 / 22.0, 1.0}
	got := unifrac(abnd, tree, false)
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unifrac(%v, %q, false)=%v, want %v",
			abnd, treeText, got, want)
	}
}

func TestUniFrac_weighted(t *testing.T) {
	treeText := "((s1:1,s2:3):2,(s3:2,s4:5):1);"
	tree, err := newick.NewReader(strings.NewReader(treeText)).Read()
	if err != nil {
		t.Fatal("failed to parse tree:", err)
	}
	abnd := []map[string]float64{
		{"s1": 4, "s2": 1},
		{"s3": 3, "s2": 2},
	}
	want := []float64{22.0 / 36.0}
	got := unifrac(abnd, tree, true)
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("unifrac(%v, %q, true)=%v, want %v",
			abnd, treeText, got, want)
	}
}
