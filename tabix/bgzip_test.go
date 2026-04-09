package tabix

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"io"
	"strings"
	"testing"
)

func TestBGZipWriter_SmallWrite(t *testing.T) {
	var buf bytes.Buffer
	w := NewBGZipWriter(&buf)

	data := []byte("hello bgzip\n")
	if _, err := w.Write(data); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Should be decompressible with standard gzip.
	gr, err := gzip.NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	got, err := io.ReadAll(gr)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != string(data) {
		t.Errorf("got %q, want %q", got, data)
	}
}

func TestBGZipWriter_LargeWrite(t *testing.T) {
	var buf bytes.Buffer
	w := NewBGZipWriter(&buf)

	// Write more than one block (>64KB).
	line := "ACGTACGTACGTACGTACGTACGTACGTACGT\n"
	total := 0
	for total < bgzfMaxBlockSize*2+1000 {
		n, err := w.Write([]byte(line))
		if err != nil {
			t.Fatal(err)
		}
		total += n
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Decompress with standard gzip (handles concatenated members).
	result := &bytes.Buffer{}
	r := bytes.NewReader(buf.Bytes())
	for {
		gr, err := gzip.NewReader(r)
		if err != nil {
			break
		}
		if _, err := io.Copy(result, gr); err != nil {
			t.Fatal(err)
		}
		gr.Close()
	}

	if result.Len() != total {
		t.Errorf("decompressed size = %d, want %d", result.Len(), total)
	}
}

func TestBGZipWriter_BlockStructure(t *testing.T) {
	var buf bytes.Buffer
	w := NewBGZipWriter(&buf)

	// Write exactly one full block + some extra.
	data := make([]byte, bgzfMaxBlockSize+100)
	for i := range data {
		data[i] = byte('A' + i%26)
	}
	if _, err := w.Write(data); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Parse the output to verify block structure.
	raw := buf.Bytes()
	blocks := 0
	pos := 0
	for pos < len(raw) {
		if pos+bgzfHeaderSize > len(raw) {
			t.Fatalf("truncated block at offset %d", pos)
		}
		// Check gzip magic.
		if raw[pos] != 0x1f || raw[pos+1] != 0x8b {
			t.Fatalf("bad magic at offset %d: %x %x", pos, raw[pos], raw[pos+1])
		}
		// Check FEXTRA flag.
		if raw[pos+3]&0x04 == 0 {
			t.Fatalf("FEXTRA not set at offset %d", pos)
		}
		// Check BC subfield.
		if raw[pos+12] != 'B' || raw[pos+13] != 'C' {
			t.Fatalf("BC subfield not found at offset %d", pos)
		}
		// Read BSIZE.
		bsize := int(binary.LittleEndian.Uint16(raw[pos+16:pos+18]))
		blockLen := bsize + 1
		if pos+blockLen > len(raw) {
			t.Fatalf("block extends past end: offset=%d, blockLen=%d, total=%d", pos, blockLen, len(raw))
		}
		pos += blockLen
		blocks++
	}

	// Should have 3 blocks: full block + partial block + EOF block.
	if blocks != 3 {
		t.Errorf("got %d blocks, want 3", blocks)
	}
}

func TestBGZipWriter_EOFBlock(t *testing.T) {
	var buf bytes.Buffer
	w := NewBGZipWriter(&buf)
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Empty file should just be the EOF block.
	if !bytes.Equal(buf.Bytes(), bgzfEOFBlock) {
		t.Errorf("empty close should produce EOF block only, got %d bytes", buf.Len())
	}
}

func TestBGZipWriter_GzipCompatibility(t *testing.T) {
	var buf bytes.Buffer
	w := NewBGZipWriter(&buf)

	// Write a variety of data.
	lines := []string{
		"chr1\t100\t200\tread1\t0\t+\n",
		"chr1\t150\t250\tread2\t0\t-\n",
		"chr2\t0\t1000\tread3\t0\t+\n",
	}
	expected := strings.Join(lines, "")
	for _, line := range lines {
		if _, err := w.Write([]byte(line)); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Standard gzip should decompress it fully.
	gr, err := gzip.NewReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatal(err)
	}
	got, err := io.ReadAll(gr)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != expected {
		t.Errorf("got %q, want %q", got, expected)
	}
}
