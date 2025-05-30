package main

import (
	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/mash/v2"
	"github.com/fluhus/gostuff/jio"
)

// Creates a kmer sketch for the given fasta file.
func sketchFile(fas []*fasta.Fasta, fout string) error {
	mh := mash.New(int(*n), int(*k)).AddSeq(iterSeqs(fas))
	if err := jio.Write(fout, mh); err != nil {
		return err
	}
	return nil
}

// Loads the sketches saved in the given file list.
func loadSketches(files []string) ([]*mash.Masher, error) {
	result := make([]*mash.Masher, 0, len(files))
	for _, file := range files {
		x := &mash.Masher{}
		err := jio.Read(file, &x)
		if err != nil {
			return nil, err
		}
		result = append(result, x)
	}
	return result, nil
}
