package main

import (
	"sort"

	"github.com/fluhus/frackyfrac/sketch"
	"github.com/fluhus/gostuff/gzipf"
	"github.com/fluhus/gostuff/jsonf"
	"github.com/spaolacci/murmur3"
)

// Creates a kmer sketch for the given fasta file.
func sketchFile(fin, fout string) error {
	f, err := gzipf.Open(fin)
	if err != nil {
		return err
	}
	skch := sketch.New(int(*n))
	h := murmur3.New64()
	err = (iterKmers(f, int(*k), func(kmer []byte) {
		h.Reset()
		h.Write(kmer)
		skch.Add(h.Sum64())
	}))
	f.Close()
	if err != nil {
		return err
	}
	hashes := skch.View()
	sort.Slice(hashes, func(i, j int) bool {
		return hashes[i] < hashes[j]
	})
	err = jsonf.Save(fout, hashes)
	if err != nil {
		return err
	}
	return nil
}

// Loads the sketches saved in the given file list.
func loadSketches(files []string) ([][]uint64, error) {
	result := make([][]uint64, 0, len(files))
	// cl := util.NewCooler()
	for _, file := range files {
		x := make([]uint64, 0, *n)
		err := jsonf.Load(file, &x)
		if err != nil {
			return nil, err
		}
		result = append(result, x)
	}
	return result, nil
}
