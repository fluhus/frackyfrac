package ppln

import (
	"fmt"
	"math/rand"
	"testing"
	"time"
)

func TestSerial(t *testing.T) {
	for _, nt := range []int{1, 2, 4, 8} {
		t.Run(fmt.Sprint(nt), func(t *testing.T) {
			n := nt * 100
			var result []int
			Serial(nt, func(c chan<- interface{}) {
				for i := 0; i < n; i++ {
					time.Sleep(time.Millisecond * time.Duration(rand.Intn(5)))
					c <- i
				}
			}, func(i interface{}) interface{} {
				ii := i.(int)
				return ii * ii
			}, func(i interface{}) {
				result = append(result, i.(int))
			})
			for i := range result {
				if result[i] != i*i {
					t.Errorf("result[%d]=%d, want %d", i, result[i], i*i)
				}
			}
		})
	}
}
