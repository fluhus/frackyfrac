package main

import (
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/minhash"
	"github.com/spaolacci/murmur3"
)

// Creates a kmer sketch for the given fasta file.
func sketchFile(fin, fout string) error {
	mh := minhash.New[uint64](int(*n))
	h := murmur3.New64()
	for kmer, err := range iterKmers(fin, int(*k)) {
		if err != nil {
			return err
		}
		h.Reset()
		h.Write(kmer)
		mh.Push(h.Sum64())
	}
	mh.Sort()
	if err := jio.Save(fout, mh); err != nil {
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
