package parser

import (
	"reflect"
	"strings"
	"testing"
)

func TestParseAbundance(t *testing.T) {
	input := "   aa  bbbb    \n1\t2\n 3  \t  4 \t\n"
	want := []map[string]float64{
		{"aa": 1, "bbbb": 2},
		{"aa": 3, "bbbb": 4},
	}
	var got []map[string]float64
	err := ParseAbundance(strings.NewReader(input), 1,
		func(m map[string]float64) {
			got = append(got, m)
		})
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
	var got []map[string]float64
	err := ParseSparseAbundance(strings.NewReader(input), 1,
		func(m map[string]float64) {
			got = append(got, m)
		})
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
