#!/usr/bin/env python 
import pysam
import random
import argparse
import sys

ERROR_TYPES = ["pnext", "tlen", "flags", "rnext", "mc_tag"]

def corrupt_sam(input_bam, output_bam, error_types, rate=0.1, seed=None):
    if seed is not None:
        random.seed(seed)

    with pysam.AlignmentFile(input_bam, "r") as fin, \
         pysam.AlignmentFile(output_bam, "w", header=fin.header) as fout:

        refs = list(fin.header.references)
        corrupted = 0

        for read in fin:
            if not read.is_paired or random.random() >= rate:
                fout.write(read)
                continue

            err = random.choice(error_types)

            if err == "pnext" and read.next_reference_start:
                read.next_reference_start += random.randint(50, 500)

            elif err == "tlen":
                read.template_length = 0

            elif err == "flags":
                read.flag ^= 0x20  # flip mate-reverse

            elif err == "rnext":
                wrong = [r for r in refs if r != read.next_reference_name]
                if wrong:
                    ref = random.choice(wrong)
                    read.next_reference_name = ref
                    read.next_reference_id = fin.header.get_tid(ref)

            elif err == "mc_tag":
                tags = {k: v for k, v in read.get_tags() if k != "MC"}
                read.set_tags(list(tags.items()))

            fout.write(read)
            corrupted += 1

    return corrupted


def main():
    parser = argparse.ArgumentParser(
        description="Inject fixmate-relevant errors into a BAM/SAM file."
    )
    parser.add_argument("input", help="Input BAM/SAM (name-sorted recommended)")
    parser.add_argument("output", help="Output BAM/SAM")
    parser.add_argument(
        "-e", "--errors",
        nargs="+",
        choices=ERROR_TYPES,
        default=ERROR_TYPES,
        metavar="ERR",
        help=f"Error types to inject (default: all). Choices: {', '.join(ERROR_TYPES)}",
    )
    parser.add_argument(
        "-r", "--rate",
        type=float,
        default=0.1,
        help="Fraction of paired reads to corrupt (default: 0.1)",
    )
    parser.add_argument(
        "-s", "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility",
    )
    args = parser.parse_args()

    if not 0 < args.rate <= 1:
        parser.error("--rate must be between 0 (exclusive) and 1")

    n = corrupt_sam(args.input, args.output, args.errors, args.rate, args.seed)
    print(f"Done: {n} reads corrupted → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
