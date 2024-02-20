package common

import (
	"slices"
	"testing"
)

func TestIterPairs(t *testing.T) {
	s := []int{1, 2, 4, 8}
	want := [][2]int{{2, 1}, {4, 1}, {4, 2}, {8, 1}, {8, 2}, {8, 4}}
	var got [][2]int
	for i, _ := range IterPairs(s) {
		got = append(got, i)
	}
	if !slices.Equal(got, want) {
		t.Fatalf("pairsInput(%v)=%v, want %v", s, got, want)
	}
}
