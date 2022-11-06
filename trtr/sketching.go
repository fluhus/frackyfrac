package main

import (
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/minhash"
	"github.com/spaolacci/murmur3"
)

// Creates a kmer sketch for the given fasta file.
func sketchFile(fin, fout string) error {
	f, err := aio.Open(fin)
	if err != nil {
		return err
	}
	mh := minhash.New[uint64](int(*n))
	h := murmur3.New64()
	err = (iterKmers(f, int(*k), func(kmer []byte) {
		h.Reset()
		h.Write(kmer)
		mh.Push(h.Sum64())
	}))
	f.Close()
	if err != nil {
		return err
	}
	mh.Sort()
	err = jio.Save(fout, mh)
	if err != nil {
		return err
	}
	return nil
}

// Loads the sketches saved in the given file list.
func loadSketches(files []string) ([]*minhash.MinHash[uint64], error) {
	result := make([]*minhash.MinHash[uint64], 0, len(files))
	for _, file := range files {
		x := minhash.New[uint64](1)
		err := jio.Load(file, &x)
		if err != nil {
			return nil, err
		}
		result = append(result, x)
	}
	return result, nil
}
