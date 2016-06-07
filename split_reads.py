#!/bin/env python
import argparse
import math
import operator
import pysam
import re


def make_windows(length, window, slide):
    """
    For a given length, return an iterator for intervals of length `window` with
    a slide of `slide`.

    >>> list(make_windows(8, 4, 0))
    [(0, 4), (4, 8)]
    >>> list(make_windows(8, 5, 0))
    [(0, 5), (5, 8)]
    >>> list(make_windows(8, 8, 0))
    [(0, 8)]
    >>> list(make_windows(8, 4, 2))
    [(0, 4), (2, 6), (4, 8)]
    >>> list(make_windows(8, 5, 2))
    [(0, 5), (2, 7), (4, 8)]
    >>> list(make_windows(7, 8, 0))
    [(0, 7)]
    """
    if slide == 0:
        windows = xrange(0, length, window)
    else:
        windows = xrange(0, length, slide)

    for start in windows:
        yield (start, min(start + window, length))

        # At most, only output one window at the end of the sequence.
        if length <= start + window:
            break


def fragment_sequence(sequence, window, slide=0):
    """Fragment a given sequence to the requested window length without a slide.

    >>> fragment_sequence("ACTGACTG", 4, 0)
    ['ACTG', 'ACTG']
    >>> fragment_sequence("ACTGACTG", 5, 0)
    ['ACTGA', 'CTG']
    >>> fragment_sequence("ACTGACTG", 8, 0)
    ['ACTGACTG']

    Fragment a given sequence to the requested window length with a slide.

    >>> fragment_sequence("ACTGACTG", 4, 2)
    ['ACTG', 'TGAC', 'ACTG']
    >>> fragment_sequence("ACTGACTG", 5, 2)
    ['ACTGA', 'TGACT', 'ACTG']

    Remove gap bases from input sequence and return the longest non-gap
    fragment. Don't return any sequence if the entire input is gap bases.

    >>> fragment_sequence("NNNNNNNN", 4, 2)
    []
    >>> fragment_sequence("ACTGNNNN", 4, 2)
    ['ACTG']
    >>> fragment_sequence("ACTGNNTA", 4, 2)
    ['ACTG']
    >>> fragment_sequence("ACNNACTA", 4, 2)
    ['ACTA']
    """
    # Check sequence for gap bases and keep the longest of the non-gap pieces in
    # the sequence.
    sequences = []
    sequence_pieces = [(piece, len(piece)) for piece in re.split("N+", sequence) if len(piece) > 0]

    if len(sequence_pieces) > 0:
        sorted_sequence_pieces = sorted(sequence_pieces, key=operator.itemgetter(1), reverse=True)
        sequence = sorted_sequence_pieces[0][0]
        sequence_length = sorted_sequence_pieces[0][1]

        if sequence_length > window:
            # Split sequence into two or more reads of the given length.
            window_ranges = make_windows(sequence_length, window, slide)
            sequences = [sequence[start:end] for start, end in window_ranges]
        else:
            # Output the sequence as is.
            sequences = [sequence]

    return sequences


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split reads from FASTQ or BAM into fragments of a given size and output fragments as FASTQ records")
    parser.add_argument("input", help="FASTQ or BAM with sequences to fragment")
    parser.add_argument("window", type=int, help="length to fragment each sequence to")
    parser.add_argument("--slide", type=int, default=0, help="length to slide the given window size across the input sequences")
    parser.add_argument("--full_length_only", action="store_true", help="omit sequences that are shorter than the requested window size")
    args = parser.parse_args()

    if args.input.endswith(".bam"):
        input_file = pysam.AlignmentFile(args.input, check_header=False, check_sq=False)
        is_bam = True
    else:
        input_file = pysam.FastqFile(args.input)
        is_bam = False

    for record in input_file:
        if is_bam:
            record_name = "%s_%s" % (record.qname, (record.is_read1 and "1" or "2"))
            sequence = record.seq
            quality = record.qual
        else:
            record_name = record.name.replace("/", "_")
            sequence = record.sequence
            quality = record.quality

        sequences = fragment_sequence(sequence, args.window, args.slide)
        qualities = fragment_sequence(quality, args.window, args.slide)

        for i in xrange(len(sequences)):
            if not args.full_length_only or len(sequences[i]) == args.window:
                print "@%s_%s" % (record_name, i)
                print sequences[i]
                print "+"
                print qualities[i]
