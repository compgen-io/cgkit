package tabix

import (
	"compress/flate"
	"encoding/binary"
	"hash/crc32"
	"io"
	"os"
)

const (
	bgzfMaxBlockSize = 65536 // 64KB max uncompressed block size
	bgzfHeaderSize   = 18   // gzip header (10) + XLEN (2) + extra field (6)
	bgzfFooterSize   = 8    // CRC32 (4) + ISIZE (4)
	bgzfOverhead     = bgzfHeaderSize + bgzfFooterSize // 26 bytes
)

// bgzfEOFBlock is the standard 28-byte empty BGZF EOF marker block.
var bgzfEOFBlock = []byte{
	0x1f, 0x8b, 0x08, 0x04, // gzip magic, deflate, FEXTRA
	0x00, 0x00, 0x00, 0x00, // mtime
	0x00, 0xff, // xfl, OS=unknown
	0x06, 0x00, // XLEN=6
	0x42, 0x43, // SI1='B', SI2='C'
	0x02, 0x00, // SLEN=2
	0x1b, 0x00, // BSIZE=27 (block size - 1)
	0x03, 0x00, // empty deflate block
	0x00, 0x00, 0x00, 0x00, // CRC32=0
	0x00, 0x00, 0x00, 0x00, // ISIZE=0
}

// BGZipWriter writes bgzip (blocked gzip) compressed data. Each block is an
// independent gzip member with a BC extra field containing the block size,
// enabling random access via virtual offsets. The output is compatible with
// standard gzip decompression and can be indexed with tabix.
type BGZipWriter struct {
	w      io.Writer
	f      *os.File // non-nil if opened by NewBGZipFile
	buf    []byte   // uncompressed data buffer
	n      int      // bytes in buf
	level  int      // deflate compression level
	closed bool
}

// NewBGZipWriter creates a new BGZipWriter that writes to w. An optional
// compression level may be provided (flate.NoCompression through
// flate.BestCompression); if omitted, flate.DefaultCompression is used.
func NewBGZipWriter(w io.Writer, level ...int) *BGZipWriter {
	lvl := flate.DefaultCompression
	if len(level) > 0 {
		lvl = level[0]
	}
	return &BGZipWriter{
		w:     w,
		buf:   make([]byte, bgzfMaxBlockSize),
		level: lvl,
	}
}

// NewBGZipFile creates a new BGZipWriter that writes to the named file.
// The file is created (or truncated) and will be closed when the writer
// is closed. An optional compression level may be provided; if omitted,
// flate.DefaultCompression is used.
func NewBGZipFile(filename string, level ...int) (*BGZipWriter, error) {
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	w := NewBGZipWriter(f, level...)
	w.f = f
	return w, nil
}

// Write buffers data and flushes complete blocks. Each full 64KB buffer is
// compressed and written as an independent bgzip block.
func (b *BGZipWriter) Write(p []byte) (int, error) {
	written := 0
	for len(p) > 0 {
		space := bgzfMaxBlockSize - b.n
		if space > len(p) {
			space = len(p)
		}
		copy(b.buf[b.n:], p[:space])
		b.n += space
		p = p[space:]
		written += space

		if b.n == bgzfMaxBlockSize {
			if err := b.flushBlock(); err != nil {
				return written, err
			}
		}
	}
	return written, nil
}

// Flush writes any buffered data as a bgzip block.
func (b *BGZipWriter) Flush() error {
	if b.n > 0 {
		return b.flushBlock()
	}
	return nil
}

// Close flushes remaining data, writes the EOF block, and closes the
// underlying file if the writer was created with NewBGZipFile.
func (b *BGZipWriter) Close() error {
	if b.closed {
		return nil
	}
	b.closed = true
	if err := b.Flush(); err != nil {
		return err
	}
	if _, err := b.w.Write(bgzfEOFBlock); err != nil {
		return err
	}
	if b.f != nil {
		return b.f.Close()
	}
	return nil
}

// flushBlock compresses the current buffer as a single bgzip block.
func (b *BGZipWriter) flushBlock() error {
	data := b.buf[:b.n]
	b.n = 0

	// Deflate compress into a temporary buffer.
	var compressed []byte
	{
		// Use a growable buffer via a simple writer.
		var cw byteSliceWriter
		fw, err := flate.NewWriter(&cw, b.level)
		if err != nil {
			return err
		}
		if _, err := fw.Write(data); err != nil {
			return err
		}
		if err := fw.Close(); err != nil {
			return err
		}
		compressed = cw.buf
	}

	// Calculate block size: header + compressed + footer.
	bsize := uint16(bgzfOverhead + len(compressed) - 1)
	checksum := crc32.ChecksumIEEE(data)
	isize := uint32(len(data))

	// Write gzip header with BGZF extra field.
	var header [bgzfHeaderSize]byte
	header[0] = 0x1f // gzip magic
	header[1] = 0x8b
	header[2] = 0x08 // deflate
	header[3] = 0x04 // FEXTRA
	// bytes 4-7: mtime = 0
	header[8] = 0x00 // xfl
	header[9] = 0xff // OS = unknown
	binary.LittleEndian.PutUint16(header[10:12], 6) // XLEN
	header[12] = 'B' // SI1
	header[13] = 'C' // SI2
	binary.LittleEndian.PutUint16(header[14:16], 2) // SLEN
	binary.LittleEndian.PutUint16(header[16:18], bsize)

	if _, err := b.w.Write(header[:]); err != nil {
		return err
	}
	if _, err := b.w.Write(compressed); err != nil {
		return err
	}

	// Write footer: CRC32 + ISIZE.
	var footer [bgzfFooterSize]byte
	binary.LittleEndian.PutUint32(footer[0:4], checksum)
	binary.LittleEndian.PutUint32(footer[4:8], isize)
	_, err := b.w.Write(footer[:])
	return err
}

// byteSliceWriter is a minimal writer that appends to a byte slice.
type byteSliceWriter struct {
	buf []byte
}

func (w *byteSliceWriter) Write(p []byte) (int, error) {
	w.buf = append(w.buf, p...)
	return len(p), nil
}
