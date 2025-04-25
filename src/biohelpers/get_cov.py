import argparse
import sys

import pysam
import pysamstats


def main():
    parser = argparse.ArgumentParser(
        description="Calculate genomic regions with mean coverage ≥ threshold",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example: python get_cov.py -s input.bam -m 10 -o output.bed",
    )
    parser.add_argument(
        "-s",
        "--sam",
        required=True,
        help="Input SAM/BAM file path (must be sorted and indexed)",
    )
    parser.add_argument(
        "-m",
        "--min_coverage",
        type=int,
        default=0,
        help="Minimum mean coverage threshold (default: 0)",
    )
    parser.add_argument("-o", "--output", help="Output file path (default: stdout)")
    args = parser.parse_args()

    try:
        # Setup output handler
        output_file = args.output
        fh = open(output_file, "w") if output_file else sys.stdout
        fh.write("chr\tstart\tend\tmean_cov\n")

        with pysam.AlignmentFile(args.sam, "rb") as bam:
            if not bam.references:
                raise ValueError("BAM file lacks reference sequences")

            # Process each chromosome
            for chrom in bam.references:
                chrom_len = bam.get_reference_length(chrom)
                current_start = None
                current_sum = 0
                current_count = 0

                # Generate coverage statistics
                for record in pysamstats.stat_coverage(
                    bam, chrom=chrom, start=1, end=chrom_len, pad=True
                ):
                    pos = record["pos"]
                    cov = record["reads_all"]

                    if cov >= args.min_coverage:
                        if current_start is None:
                            current_start = pos
                        current_sum += cov
                        current_count += 1
                    else:
                        if current_start is not None:
                            end = pos - 1
                            mean_cov = current_sum / current_count
                            fh.write(
                                f"{chrom}\t{current_start}\t{end}\t{mean_cov:.2f}\n"
                            )
                            current_start = None
                            current_sum = 0
                            current_count = 0

                # Handle remaining region at chromosome end
                if current_start is not None:
                    end = chrom_len
                    mean_cov = current_sum / current_count
                    fh.write(f"{chrom}\t{current_start}\t{end}\t{mean_cov:.2f}\n")

        # Close file handle if output to file
        if output_file:
            fh.close()

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
