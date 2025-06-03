#!/usr/bin/env python3
"""
run_fastp.py - Batch processing tool for paired-end sequencing data quality control using fastp

Author: Your Name
Date: 2025-06-03
"""

import argparse
import glob
import logging
import os
import subprocess
import sys
from pathlib import Path


def setup_logging(log_file=None, verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s - %(levelname)s - %(message)s"

    if log_file:
        logging.basicConfig(
            level=level,
            format=format_str,
            handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
        )
    else:
        logging.basicConfig(level=level, format=format_str)


def get_sample_names(input_dir, pattern="*.R1.raw.fastq.gz"):
    """Get sample names list from input directory"""
    r1_files = glob.glob(os.path.join(input_dir, pattern))
    samples = []

    for r1_file in r1_files:
        basename = os.path.basename(r1_file)
        # Remove .R1.raw.fastq.gz suffix to get sample name
        sample_name = basename.replace(".R1.raw.fastq.gz", "")
        samples.append(sample_name)

    return sorted(samples)


def check_input_files(input_dir, sample_name):
    """Check if input files exist"""
    r1_file = os.path.join(input_dir, f"{sample_name}.R1.raw.fastq.gz")
    r2_file = os.path.join(input_dir, f"{sample_name}.R2.raw.fastq.gz")

    if os.path.exists(r1_file) and os.path.exists(r2_file):
        return r1_file, r2_file
    else:
        return None, None


def run_fastp_single(sample_name, input_dir, output_dir, report_dir, args):
    """Run fastp for single sample"""
    logging.info(f"Processing sample: {sample_name}")

    # Check input files
    input_r1, input_r2 = check_input_files(input_dir, sample_name)
    if not input_r1 or not input_r2:
        logging.error(f"Missing input files for sample: {sample_name}")
        return False

    # Define output files
    output_r1 = os.path.join(output_dir, f"{sample_name}.R1.clean.fastq.gz")
    output_r2 = os.path.join(output_dir, f"{sample_name}.R2.clean.fastq.gz")
    html_report = os.path.join(report_dir, f"{sample_name}.fastp.html")
    json_report = os.path.join(report_dir, f"{sample_name}.fastp.json")
    log_file = os.path.join(report_dir, f"{sample_name}.fastp.log")

    logging.info(f"  Input files: {input_r1}, {input_r2}")
    logging.info(f"  Output files: {output_r1}, {output_r2}")

    # Build fastp command
    cmd = [
        "fastp",
        "-i",
        input_r1,
        "-I",
        input_r2,
        "-o",
        output_r1,
        "-O",
        output_r2,
        "-h",
        html_report,
        "-j",
        json_report,
        "--thread",
        str(args.threads),
        "--qualified_quality_phred",
        str(args.quality_threshold),
        "--unqualified_percent_limit",
        str(args.unqualified_percent),
        "--n_base_limit",
        str(args.n_base_limit),
        "--length_required",
        str(args.min_length),
        "--cut_window_size",
        str(args.cut_window_size),
        "--cut_mean_quality",
        str(args.cut_mean_quality),
    ]

    # Add optional parameters
    if args.detect_adapter:
        cmd.append("--detect_adapter_for_pe")
    if args.correction:
        cmd.append("--correction")
    if args.cut_front:
        cmd.append("--cut_front")
    if args.cut_tail:
        cmd.append("--cut_tail")
    if args.fix_mgi_id:
        cmd.append("--fix_mgi_id")
    if args.verbose:
        cmd.append("--verbose")

    # Run fastp
    try:
        with open(log_file, "w") as log_f:
            result = subprocess.run(
                cmd, stdout=log_f, stderr=subprocess.STDOUT, check=False
            )

        if result.returncode == 0:
            logging.info(f"  ✓ {sample_name} processing completed")
            return True
        else:
            logging.warning(
                f"  ✗ {sample_name} processing failed, trying error-tolerant mode..."
            )
            return run_fastp_retry(
                sample_name,
                input_r1,
                input_r2,
                output_r1,
                output_r2,
                html_report,
                json_report,
                report_dir,
                args,
            )

    except Exception as e:
        logging.error(f"  ✗ {sample_name} processing error: {e}")
        return False


def run_fastp_retry(
    sample_name,
    input_r1,
    input_r2,
    output_r1,
    output_r2,
    html_report,
    json_report,
    report_dir,
    args,
):
    """Retry fastp with error-tolerant mode"""
    retry_log = os.path.join(report_dir, f"{sample_name}.fastp_retry.log")

    # Build error-tolerant mode command
    cmd = [
        "fastp",
        "-i",
        input_r1,
        "-I",
        input_r2,
        "-o",
        output_r1,
        "-O",
        output_r2,
        "-h",
        html_report,
        "-j",
        json_report,
        "--thread",
        str(max(1, args.threads // 2)),  # Reduce thread count
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--fix_mgi_id",
        "--verbose",
    ]

    try:
        with open(retry_log, "w") as log_f:
            result = subprocess.run(
                cmd, stdout=log_f, stderr=subprocess.STDOUT, check=False
            )

        if result.returncode == 0:
            logging.info(f"  ✓ {sample_name} error-tolerant mode processing completed")
            return True
        else:
            logging.error(f"  ✗ {sample_name} completely failed, please check raw data")
            return False

    except Exception as e:
        logging.error(f"  ✗ {sample_name} error-tolerant mode processing error: {e}")
        return False


def run_fastp(args):
    """Main processing function"""
    # Setup logging
    log_file = os.path.join(args.output_dir, "run_fastp.log") if args.log else None
    setup_logging(log_file, args.verbose)

    logging.info("Starting fastp batch processing")
    logging.info(f"Input directory: {args.input_dir}")
    logging.info(f"Output directory: {args.output_dir}")

    # Create output directories
    output_dir = args.output_dir
    report_dir = os.path.join(output_dir, "fastp_reports")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(report_dir).mkdir(parents=True, exist_ok=True)

    # Get sample list
    samples = get_sample_names(args.input_dir, args.pattern)

    if not samples:
        logging.error(f"No files matching {args.pattern} found in {args.input_dir}")
        return False

    logging.info(f"Found {len(samples)} samples: {', '.join(samples)}")
    logging.info("=" * 50)

    # Process each sample
    success_count = 0
    total_count = len(samples)

    for i, sample in enumerate(samples, 1):
        logging.info(f"Progress: {i}/{total_count}")
        if run_fastp_single(sample, args.input_dir, output_dir, report_dir, args):
            success_count += 1
        logging.info("-" * 30)

    # Output summary
    logging.info("=" * 50)
    logging.info(f"Processing completed! Success: {success_count}/{total_count}")
    logging.info(f"Clean data saved in: {output_dir}")
    logging.info(f"QC reports saved in: {report_dir}")

    return success_count == total_count


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Batch processing paired-end sequencing data using fastp",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i 01.data -o 02.clean_data
  %(prog)s -i raw_data -o clean_data -t 8 -q 25 --verbose
  %(prog)s -i input_dir -o output_dir -p "*.R1.fastq.gz" --log
        """,
    )

    # Required arguments
    parser.add_argument(
        "-i",
        "--input-dir",
        required=True,
        help="Input directory containing raw fastq.gz files",
    )

    parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory for clean data"
    )

    # File pattern arguments
    parser.add_argument(
        "-p", "--pattern", default="*.R1.raw.fastq.gz", help="R1 file matching pattern"
    )

    # fastp QC parameters
    parser.add_argument(
        "-t", "--threads", type=int, default=16, help="Number of threads to use"
    )

    parser.add_argument(
        "-q",
        "--quality-threshold",
        type=int,
        default=20,
        help="Quality value threshold",
    )

    parser.add_argument(
        "-u",
        "--unqualified-percent",
        type=int,
        default=40,
        help="Percentage of low quality bases allowed",
    )

    parser.add_argument(
        "-n", "--n-base-limit", type=int, default=10, help="Number of N bases allowed"
    )

    parser.add_argument(
        "-l",
        "--min-length",
        type=int,
        default=50,
        help="Minimum sequence length required",
    )

    parser.add_argument(
        "--cut-window-size",
        type=int,
        default=4,
        help="Sliding window size for trimming",
    )

    parser.add_argument(
        "--cut-mean-quality",
        type=int,
        default=20,
        help="Mean quality threshold for sliding window",
    )

    # Boolean parameters
    parser.add_argument(
        "--detect-adapter",
        action="store_true",
        default=True,
        help="Auto detect adapter for PE data",
    )

    parser.add_argument(
        "--correction", action="store_true", default=True, help="Enable base correction"
    )

    parser.add_argument(
        "--cut-front", action="store_true", default=True, help="Trim front of sequences"
    )

    parser.add_argument(
        "--cut-tail", action="store_true", default=True, help="Trim tail of sequences"
    )

    parser.add_argument(
        "--fix-mgi-id",
        action="store_true",
        default=True,
        help="Fix MGI sequencer ID issues",
    )

    # Logging parameters
    parser.add_argument("--log", action="store_true", help="Save log to file")

    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    # Parse arguments
    args = parser.parse_args()

    # Check input directory
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist")
        sys.exit(1)

    # Run main program
    success = run_fastp(args)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
