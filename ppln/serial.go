// Package ppln provides a generic multi-goroutine pipeline.
package ppln

import (
	"container/heap"
	"fmt"
	"sync"
)

const (
	// Size of pipeline channel buffers per goroutine.
	chanLenCoef = 1000
)

// Serial starts a multi-goroutine transformation pipeline that maintains the
// order of the inputs.
//
// Pusher sends the inputs into the given channel. It should not close it.
// Transform takes an input, processes it and returns the result. Puller acts on
// a single result.
//
// Puller will be called on the results by the same order of pusher's input.
func Serial(ngoroutines int,
	pusher func(push func(interface{}), s Stopper),
	mapper func(a interface{}, s Stopper) interface{},
	puller func(a interface{}, s Stopper)) {
	if ngoroutines < 1 {
		panic(fmt.Sprintf("bad number of goroutines: %d", ngoroutines))
	}

	stopper := make(Stopper)

	// An optimization for a single thread.
	if ngoroutines == 1 {
		pusher(func(a interface{}) {
			puller(mapper(a, stopper), stopper)
		}, stopper)
		return
	}

	push := make(chan serialItem, ngoroutines*chanLenCoef)
	pull := make(chan serialItem, ngoroutines*chanLenCoef)
	wait := &sync.WaitGroup{}
	wait.Add(ngoroutines)

	go func() {
		i := 0
		pusher(func(a interface{}) {
			push <- serialItem{i, a}
			i++
		}, stopper)
		close(push)
	}()
	for i := 0; i < ngoroutines; i++ {
		go func() {
			for item := range push {
				if stopper.Stopped() {
					break
				}
				pull <- serialItem{item.i, mapper(item.data, stopper)}
			}
			for range push { // Drain channel.
			}
			wait.Done()
		}()
	}
	go func() {
		items := &serialHeap{}
		for item := range pull {
			if stopper.Stopped() {
				break
			}
			items.put(item)
			for items.ok() {
				puller(items.pop(), stopper)
			}
		}
		for range pull { // Drain channel.
		}
		wait.Done()
	}()

	wait.Wait() // Wait for workers.
	wait.Add(1)
	close(pull)
	wait.Wait() // Wait for pull.
}

// A Stopper is used in pipelines to communicate that work should stop.
type Stopper chan struct{}

// Stop sets Stopped to true.
func (s Stopper) Stop() {
	close(s)
}

// Stopped returns whether Stop was called.
func (s Stopper) Stopped() bool {
	select {
	case <-s:
		return true
	default:
		return false
	}
}

// General data with a serial number.
type serialItem struct {
	i    int
	data interface{}
}

// A heap of serial items. Sorts by serial number.
type serialHeap struct {
	next int
	data []serialItem
}

// Checks whether the minimal element in the heap is the next in the series.
func (s *serialHeap) ok() bool {
	return len(s.data) > 0 && s.data[0].i == s.next
}

// Removes and returns the minimal element in the heap. Panics if the element
// is not the next in the series.
func (s *serialHeap) pop() interface{} {
	if !s.ok() {
		panic("get when not ok")
	}
	s.next++
	a := heap.Pop(s)
	return a.(serialItem).data
}

// Adds an item to the heap.
func (s *serialHeap) put(item serialItem) {
	if item.i < s.next {
		panic(fmt.Sprintf("put(%d) when next is %d", item.i, s.next))
	}
	heap.Push(s, item)
}

// Implementation of heap.Interface.

func (s *serialHeap) Len() int {
	return len(s.data)
}

func (s *serialHeap) Less(i, j int) bool {
	return s.data[i].i < s.data[j].i
}

func (s *serialHeap) Swap(i, j int) {
	s.data[i], s.data[j] = s.data[j], s.data[i]
}

func (s *serialHeap) Push(a interface{}) {
	s.data = append(s.data, a.(serialItem))
}

func (s *serialHeap) Pop() interface{} {
	a := s.data[len(s.data)-1]
	s.data = s.data[:len(s.data)-1]
	return a
}
