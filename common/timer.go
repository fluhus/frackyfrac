package common

import (
	"fmt"
	"os"
	"strings"
	"time"
)

// A Timer measures time during iterative processes and prints the progress on
// exponential checkpoints.
type Timer struct {
	N int
	t time.Time
	s string
}

// Indexes of checkpoints.
var checkpoints = map[int]struct{}{}

// Initializes the checkpoints variable.
func init() {
	exp := 1
	for i := 0; i < 10; i++ {
		checkpoints[exp] = struct{}{}
		checkpoints[exp*2] = struct{}{}
		checkpoints[exp*5] = struct{}{}
		exp *= 10
	}
}

// Checks if i is in the checkpoints map.
func isCheckpoint(i int) bool {
	_, ok := checkpoints[i]
	return ok
}

// Prints the progress.
func (t *Timer) print() {
	since := time.Since(t.t)
	msg := strings.ReplaceAll(t.s, "*", fmt.Sprint(t.N))
	fmt.Fprintf(os.Stderr, "\r%s (%s) %s", fmtDuration(since),
		fmtDuration(since/time.Duration(t.N)), msg)
}

// Formats a duration in constant-width format.
func fmtDuration(d time.Duration) string {
	return fmt.Sprintf("%02d:%02d:%02d.%06d",
		d/time.Hour,
		d%time.Hour/time.Minute,
		d%time.Minute/time.Second,
		d%time.Second/time.Microsecond,
	)
}

// NewTimerMessasge returns a new timer that prints msg on checkpoints.
// A '*' character in msg will be replaced with the current count.
func NewTimerMessasge(msg string) *Timer {
	return &Timer{0, time.Now(), msg}
}

// NewTimer returns a new timer without a specific message.
func NewTimer() *Timer {
	return NewTimerMessasge("*")
}

// Inc increments t's counter and prints progress if reached a checkpoint.
func (t *Timer) Inc() {
	t.N++
	if isCheckpoint(t.N) {
		t.print()
	}
}

// Done prints progress as if a checkpoint was reached.
func (t *Timer) Done() {
	if t.N == 0 {
		t.N = 1
		t.print()
		t.N = 0
	} else {
		t.print()
	}
	fmt.Fprintln(os.Stderr)
}
