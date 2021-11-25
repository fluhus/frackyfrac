package sketch

import (
	"container/heap"
	"encoding/json"
	"fmt"
)

// A Sketch is a min-hash collection. Retains the k lowest unique values out of all
// the values that were added to it.
type Sketch struct {
	h sketchHeap
	k int
}

// New returns an empty sketch that stores k values.
func New(k int) *Sketch {
	if k < 1 {
		panic(fmt.Sprintf("invalid n: %d, should be positive", k))
	}
	return &Sketch{sketchHeap{nil, map[uint64]struct{}{}}, k}
}

// Add tries to add x to the sketch. x is added only if it does not already exist
// in the sketch, and there are less than k elements lesser than x.
// Returns true if x was added and false if not.
func (s *Sketch) Add(x uint64) bool {
	if len(s.h.h) == s.k && x >= s.h.h[0] {
		// x is too large.
		return false
	}
	if _, ok := s.h.m[x]; ok {
		return false
	}
	if len(s.h.h) == s.k {
		heap.Pop(&s.h)
	}
	heap.Push(&s.h, x)
	return true
}

// View returns the underlying slice of values.
func (s *Sketch) View() []uint64 {
	return s.h.h
}

// MarshalJSON implements the json.Marshaler interface.
func (s *Sketch) MarshalJSON() ([]byte, error) {
	return json.Marshal(struct {
		K int      `json:"k"`
		H []uint64 `json:"h"`
	}{s.k, s.h.h})
}

// UnmarshalJSON implements the json.Unmarshaler interface.
func (s *Sketch) UnmarshalJSON(b []byte) error {
	var raw struct {
		K int
		H []uint64
	}
	if err := json.Unmarshal(b, &raw); err != nil {
		return err
	}
	s.k = raw.K
	s.h.h = raw.H
	s.h.m = make(map[uint64]struct{}, raw.K)
	for _, i := range raw.H {
		s.h.m[i] = struct{}{}
	}
	heap.Init(&s.h)
	return nil
}

// A heap that panics if an element is inserted twice.
type sketchHeap struct {
	h []uint64
	m map[uint64]struct{}
}

// Implement heap.Interface.

func (s *sketchHeap) Len() int {
	return len(s.h)
}

func (s *sketchHeap) Less(i, j int) bool {
	return s.h[i] > s.h[j] // Max heap.
}

func (s *sketchHeap) Swap(i, j int) {
	s.h[i], s.h[j] = s.h[j], s.h[i]
}

func (s *sketchHeap) Push(x interface{}) {
	ix := x.(uint64)
	if _, ok := s.m[ix]; ok {
		panic(fmt.Sprintf("element %d already exists in the sketch map", ix))
	}
	s.m[ix] = struct{}{}
	s.h = append(s.h, ix)
}

func (s *sketchHeap) Pop() interface{} {
	x := s.h[len(s.h)-1]
	if _, ok := s.m[x]; !ok {
		panic(fmt.Sprintf("element %d does not exist in the sketch map", x))
	}
	delete(s.m, x)
	s.h = s.h[:len(s.h)-1]
	return x
}
