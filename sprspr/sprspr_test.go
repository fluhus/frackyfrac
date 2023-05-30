package main

import (
	"bytes"
	"strings"
	"testing"

	"golang.org/x/exp/slices"
)

func TestToSparse(t *testing.T) {
	tests := []struct {
		input string
		want  string
	}{
		{"s1\n1.3", "s1:1.3"},
		{"s1\ts2\n0\t4\n3\t0\n", "s2:4\ns1:3"},
		{"s1\ts2\ts3\n4\t3\t2\n5\t0\t8\n0\t0\t10",
			"s1:4\ts2:3\ts3:2\ns1:5\ts3:8\ns3:10"},
	}
	for _, test := range tests {
		buf := bytes.NewBuffer(nil)
		if err := toSparse(strings.NewReader(test.input), buf); err != nil {
			t.Errorf("toSparse(%q) failed: %v", test.input, err)
		}
		got := strings.TrimSuffix(buf.String(), "\n")
		rows := strings.Split(got, "\n")
		for i := range rows {
			parts := strings.Split(rows[i], "\t")
			slices.Sort(parts)
			rows[i] = strings.Join(parts, "\t")
		}
		got = strings.Join(rows, "\n")
		if got != test.want {
			t.Errorf("toSparse(%q)=%q, want %q", test.input, got, test.want)
		}
	}
}
