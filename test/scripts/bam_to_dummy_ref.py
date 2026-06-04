#!/usr/bin/env python3
"""Create a dummy FASTA reference from a BAM file header.

By default the reference sequences are filled with Ns, but they keep the
correct names, lengths and order from the @SQ lines in the BAM header.
With --random, random A/C/G/T bases are emitted instead.
"""

import argparse
import random
import sys

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a dummy FASTA reference whose contigs match the @SQ "
            "entries of a BAM file. Useful for tools that need a reference "
            "but you do not care about the actual sequence content."
        )
    )
    parser.add_argument("bam", help="Input BAM file (index not required)")
    parser.add_argument("fasta", help="Output FASTA file")
    parser.add_argument(
        "--random",
        action="store_true",
        help="Fill sequences with random A/C/G/T bases instead of Ns.",
    )
    parser.add_argument(
        "--line-width",
        type=int,
        default=80,
        help="Number of bases per FASTA line (default: 80).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (only used with --random).",
    )
    return parser.parse_args()


def get_references(bam_path):
    """Return (names, lengths) preserving the order from the @SQ header."""
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except (OSError, ValueError) as exc:
        sys.exit(f"Error: could not open BAM file '{bam_path}': {exc}")

    names = bam.references
    lengths = bam.lengths
    bam.close()

    if names is None or len(names) == 0:
        sys.exit("Error: no @SQ records found in the BAM header.")

    return list(names), list(lengths)


def write_dummy_fasta(fasta_path, names, lengths, random_bases, line_width):
    """Write the dummy FASTA. line_width must be >= 1."""
    if line_width < 1:
        sys.exit("Error: --line-width must be >= 1.")

    alphabet = "ACGT" if random_bases else "N"
    full_n_line = "N" * line_width  # pre-build for speed

    with open(fasta_path, "w") as out:
        for name, length in zip(names, lengths):
            if name == "":
                sys.exit("Error: empty contig name in BAM header.")
            out.write(f">{name}\n")

            if length <= 0:
                continue

            if random_bases:
                # Generate in chunks to avoid allocating a huge string
                remaining = length
                while remaining > 0:
                    chunk = min(remaining, line_width)
                    out.write("".join(random.choices(alphabet, k=chunk)) + "\n")
                    remaining -= chunk
            else:
                # Ns: just repeat a pre-built line
                n_full = length // line_width
                n_rem = length % line_width
                for _ in range(n_full):
                    out.write(full_n_line + "\n")
                if n_rem:
                    out.write("N" * n_rem + "\n")


def main():
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    names, lengths = get_references(args.bam)
    write_dummy_fasta(
        fasta_path=args.fasta,
        names=names,
        lengths=lengths,
        random_bases=args.random,
        line_width=args.line_width,
    )

    total = sum(lengths)
    print(
        f"Wrote {len(names)} contig(s), {total:,} bp total, to '{args.fasta}'."
    )


if __name__ == "__main__":
    main()
