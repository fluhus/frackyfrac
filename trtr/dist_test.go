package main

import (
	"reflect"
	"strings"
	"testing"

	"github.com/fluhus/biostuff/formats/newick"
)

func TestDist(t *testing.T) {
	tests := []struct {
		a, b []uint64
		want float64
	}{
		{[]uint64{1, 2, 3}, []uint64{1, 2, 3}, 0},
		{[]uint64{1, 2, 3}, []uint64{2, 3, 4}, 0.5},
		{[]uint64{2, 3, 4}, []uint64{1, 2, 3}, 0.5},
		{[]uint64{1, 2, 3, 4, 5}, []uint64{1, 3, 5}, 0.4},
	}
	for _, test := range tests {
		if got := dist(test.a, test.b); got != test.want {
			t.Errorf("dist(%d,%d)=%f, want %f",
				test.a, test.b, got, test.want)
		}
	}
}

func TestMakeTree(t *testing.T) {
	sketches := [][]uint64{
		{1, 2, 3},
		{5, 6, 7},
		{5, 6, 7, 8},
		{2, 3, 4},
	}
	names := []string{"a", "b", "c", "d"}
	wantString := "((b:0.125,c:0.125):0.375,(a:0.25,d:0.25):0.25);"
	want, err := newick.NewReader(strings.NewReader(wantString)).Read()
	if err != nil {
		t.Fatal("failed to parse want string")
	}
	got := makeTree(sketches, names)
	if !reflect.DeepEqual(got, want) {
		gotText, _ := got.MarshalText()
		wantText, _ := want.MarshalText()
		t.Fatalf("makeTree(...)=%s, want %s", gotText, wantText)
	}
}
