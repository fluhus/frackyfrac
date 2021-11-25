package sketch

import (
	"reflect"
	"sort"
	"testing"
)

func TestSketch(t *testing.T) {
	tests := []struct {
		n     int
		input []uint64
		want  []uint64
	}{
		{
			3,
			[]uint64{1, 2, 2, 2, 2, 1, 1, 3, 3, 3, 1, 2, 3, 1, 3, 3, 2},
			[]uint64{1, 2, 3},
		},
		{
			3,
			[]uint64{1, 2, 3, 4, 5, 6, 7, 8, 9},
			[]uint64{1, 2, 3},
		},
		{
			3,
			[]uint64{9, 8, 7, 6, 5, 4, 3, 2, 1},
			[]uint64{1, 2, 3},
		},
		{
			5,
			[]uint64{40, 19, 55, 10, 32, 1, 100, 5, 99, 16, 16},
			[]uint64{1, 5, 10, 16, 19},
		},
	}
	for _, test := range tests {
		skch := New(test.n)
		for _, k := range test.input {
			skch.Add(k)
		}
		got := skch.h.h
		sort.Slice(got, func(i, j int) bool {
			return got[i] < got[j]
		})
		if !reflect.DeepEqual(skch.h.h, test.want) {
			t.Errorf("New(%d).Add(%v)=%v, want %v",
				test.n, test.input, got, test.want)
		}
	}
}

func TestJSON(t *testing.T) {
	input := New(5)
	input.Add(1)
	input.Add(4)
	input.Add(9)
	input.Add(16)
	input.Add(25)
	jsn, err := input.MarshalJSON()
	if err != nil {
		t.Fatalf("Sketch(1,2,3,4,5).MarshalJSON() failed: %v", err)
	}
	got := New(2)
	err = got.UnmarshalJSON(jsn)
	if err != nil {
		t.Fatalf("UnmarshalJSON(%q) failed: %v", jsn, err)
	}
	if !reflect.DeepEqual(got, input) {
		t.Fatalf("UnmarshalJSON(%q)=%v, want %v", jsn, got, input)
	}
}
