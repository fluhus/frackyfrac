package main

import "testing"

func TestIJToN(t *testing.T) {
	tests := []struct {
		i, j, n int
	}{
		{1, 0, 0}, {2, 0, 1}, {2, 1, 2}, {3, 0, 3}, {3, 1, 4}, {3, 2, 5},
	}
	for _, test := range tests {
		if n := ijToN(test.i, test.j); n != test.n {
			t.Errorf("ijToN(%v,%v)=%v, want %v",
				test.i, test.j, n, test.n)
		}
		if n := ijToN(test.j, test.i); n != test.n {
			t.Errorf("ijToN(%v,%v)=%v, want %v",
				test.j, test.i, n, test.n)
		}
	}
}

func TestIJToN_equal(t *testing.T) {
	defer func() { recover() }()
	n := ijToN(2, 2)
	t.Fatalf("ijToN(2,2)=%v, want panic", n)
}
