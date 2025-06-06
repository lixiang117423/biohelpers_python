#!/usr/bin/env python3
"""
RNA-seq Analysis Pipeline
Complete workflow including HISAT2 indexing, alignment, StringTie quantification and expression matrix merging
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


def run_command(cmd, description=""):
    """Execute shell command and handle errors"""
    print(f"Running: {description}")
    print(f"Command: {cmd}")

    try:
        result = subprocess.run(
            cmd, shell=True, check=True, capture_output=True, text=True
        )
        print(f"✓ {description} completed")
        return result
    except subprocess.CalledProcessError as e:
        print(f"✗ {description} failed")
        print(f"Error message: {e.stderr}")
        sys.exit(1)


def build_hisat2_index(genome_path, threads):
    """Build HISAT2 genome index"""
    genome_dir = os.path.dirname(genome_path)
    genome_name = os.path.splitext(os.path.basename(genome_path))[0]
    index_prefix = os.path.join(genome_dir, f"{genome_name}.hisat2.index")

    # Check if index already exists
    if os.path.exists(f"{index_prefix}.1.ht2"):
        print(f"HISAT2 index already exists: {index_prefix}")
        return index_prefix

    cmd = f"hisat2-build -p {threads} {genome_path} {index_prefix}"
    run_command(cmd, "Building HISAT2 index")
    return index_prefix


def run_hisat2_mapping(index_prefix, fastq1, fastq2, output_bam, threads):
    """Run HISAT2 alignment and sorting"""
    # Ensure output directory exists
    output_dir = os.path.dirname(output_bam)
    os.makedirs(output_dir, exist_ok=True)

    cmd = (
        f"hisat2 -x {index_prefix} -1 {fastq1} -2 {fastq2} -p {threads} | "
        f"samtools sort -@ {threads} -O BAM -o {output_bam} -"
    )

    run_command(cmd, f"HISAT2 alignment and sorting -> {output_bam}")


def run_stringtie(bam_file, gtf_file, output_gtf, threads):
    """Run StringTie quantification"""
    # Ensure output directory exists
    output_dir = os.path.dirname(output_gtf)
    os.makedirs(output_dir, exist_ok=True)

    cmd = f"stringtie -p {threads} -G {gtf_file} -o {output_gtf} -e {bam_file}"
    run_command(cmd, f"StringTie quantification -> {output_gtf}")


def extract_gtf_values(gtf_file, sample_name, output_file):
    """
    Extract gene ID, FPKM and TPM values from GTF file
    and add sample name to the last column
    """
    # Regular expression to match key-value pairs (e.g., gene_id "YLG01G000001")
    pattern = re.compile(r'(\w+)\s+"([^"]+)"')

    with open(gtf_file, "r") as infile, open(output_file, "w") as outfile:
        # 1. First write the header
        outfile.write("gene_id\ttranscript_id\tcov\tFPKM\tTPM\tsample\n")

        for line in infile:
            # Skip comment lines and lines without FPKM
            if line.startswith("#") or "FPKM" not in line:
                continue

            # Split line (keep complete attributes in column 9)
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            # Parse attributes column
            attributes = {}
            for match in pattern.finditer(cols[8]):
                key, value = match.groups()
                attributes[key] = value

            # Get required values (ensure all keys exist)
            if all(key in attributes for key in ["gene_id", "FPKM", "TPM"]):
                gene_id = attributes["gene_id"]
                transcript_id = attributes["transcript_id"]
                fpkm = attributes["FPKM"]
                tpm = attributes["TPM"]
                cov = attributes["cov"]

                # Write to output file (tab-separated)
                outfile.write(
                    f"{gene_id}\t{transcript_id}\t{cov}\t{fpkm}\t{tpm}\t{sample_name}\n"
                )


def process_single_sample(args, sample_info, index_prefix):
    """Process complete workflow for a single sample"""
    sample_name = sample_info["name"]
    fastq1 = sample_info["fastq1"]
    fastq2 = sample_info["fastq2"]

    print(f"\n{'=' * 60}")
    print(f"Processing sample: {sample_name}")
    print(f"{'=' * 60}")

    # Set file paths
    bam_file = os.path.join(args.output, f"{sample_name}.sorted.bam")
    stringtie_output = os.path.join(
        args.output, "stringtie_output", f"{sample_name}.gtf"
    )
    fpkm_output = os.path.join(args.output, "fpkm_output", f"{sample_name}.fpkm.txt")

    # 1. HISAT2 alignment
    if os.path.exists(bam_file):
        print(f"✓ BAM file already exists, skipping alignment: {bam_file}")
    else:
        run_hisat2_mapping(index_prefix, fastq1, fastq2, bam_file, args.threads)

    # 2. StringTie quantification
    if os.path.exists(stringtie_output):
        print(f"✓ GTF file already exists, skipping StringTie: {stringtie_output}")
    else:
        # run_hisat2_mapping(index_prefix, fastq1, fastq2, bam_file, args.threads)
        run_stringtie(bam_file, args.gtf, stringtie_output, args.threads)

    # 3. Extract FPKM values
    os.makedirs(os.path.dirname(fpkm_output), exist_ok=True)
    # extract_fpkm(stringtie_output, sample_name, fpkm_output)
    extract_gtf_values(stringtie_output, sample_name, fpkm_output)

    # print(f"Successfully extracted values to {output_file}")

    # 4. Remove BAM file if requested
    if args.remove.lower() in ["yes", "y"]:
        if os.path.exists(bam_file):
            os.remove(bam_file)
            print(f"✓ Removed BAM file: {bam_file}")
    else:
        print(f"✓ BAM file retained: {bam_file}")

    return fpkm_output


# def handle_duplicates(df, method="sum"):
#     """
#     Handle duplicate gene_id entries
#     method: 'sum' - sum values, 'mean' - average values, 'first' - take first, 'max' - take maximum
#     """
#     if method == "sum":
#         # For same gene_id and sample, sum FPKM and TPM values
#         df_grouped = (
#             df.groupby(["gene_id", "sample"])
#             .agg({"cov": "sum", "FPKM": "sum", "TPM": "sum"})
#             .reset_index()
#         )
#     elif method == "mean":
#         # For same gene_id and sample, average FPKM and TPM values
#         df_grouped = (
#             df.groupby(["gene_id", "sample"])
#             .agg({"cov": "mean", "FPKM": "mean", "TPM": "mean"})
#             .reset_index()
#         )
#     elif method == "first":
#         # For same gene_id and sample, take first value
#         df_grouped = df.groupby(["gene_id", "sample"]).first().reset_index()
#     elif method == "max":
#         # For same gene_id and sample, take maximum value
#         df_grouped = (
#             df.groupby(["gene_id", "sample"])
#             .agg({"cov": "max", "FPKM": "max", "TPM": "max"})
#             .reset_index()
#         )
#     else:
#         raise ValueError("method must one of 'sum', 'mean', 'first', or 'max'.")

#     return df_grouped


def merge_expression_matrix(fpkm_files, output_dir):
    """Merge expression matrix from all samples and generate three output files"""
    # print(f"\n{'=' * 60}")
    # print("Merging expression matrix")
    # print(f"{'=' * 60}")

    all_data = []

    # Read all FPKM files
    for fpkm_file in fpkm_files:
        if os.path.exists(fpkm_file):
            print(f"Read file: {fpkm_file}")
            try:
                # Read file, assuming six columns: gene_id, transcript_id, cov, FPKM, TPM, sample
                df = pd.read_csv(
                    fpkm_file,
                    sep="\t",
                    header=None,
                    names=["gene_id", "transcript_id", "cov", "FPKM", "TPM", "sample"],
                )

                # Check data types and convert
                df["cov"] = pd.to_numeric(df["cov"], errors="coerce")
                df["FPKM"] = pd.to_numeric(df["FPKM"], errors="coerce")
                df["TPM"] = pd.to_numeric(df["TPM"], errors="coerce")

                # Remove NaN values
                df = df.dropna()

                # # Check for duplicate gene_id entries
                # duplicates = df.duplicated(subset=["gene_id"], keep=False)
                # if duplicates.any():
                #     print(f"  Warning: Found {duplicates.sum()} duplicate gene_id entries")
                #     print(f"  Using sum method to handle duplicates")
                #     df = handle_duplicates(df, method="sum")

                all_data.append(df)
                # print(f"  Successfully read {len(df)} rows of data")

            except Exception as e:
                print(f"  Error reading file {fpkm_file}: {e}")
                continue
        else:
            print(f"File does not exist: {fpkm_file}")

    if not all_data:
        print("✗ No valid FPKM data found")
        return

    print(f"\nFound {len(all_data)} valid data files")

    # Merge all data
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"Total {len(combined_df)} rows after merging")

    # Get sample name list
    sample_names = sorted(combined_df["sample"].unique())
    print(f"Sample names: {', '.join(sample_names)}")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # 1. Generate all.fpkm.tpm.txt (direct merge of all data)
    print("\nGenerating all.fpkm.tpm.txt...")
    all_output_file = os.path.join(output_dir, "all.fpkm.tpm.txt")
    # Reorder columns: gene_id, transcript_id, cov, FPKM, TPM, sample
    combined_df_ordered = combined_df[
        ["gene_id", "transcript_id", "cov", "FPKM", "TPM", "sample"]
    ]
    combined_df_ordered.to_csv(all_output_file, sep="\t", index=False, header=True)
    print(f"✓ Save completed: {all_output_file} ({len(combined_df_ordered)} rows)")

    # # 2. Generate fpkm.txt (FPKM matrix)
    # print("\nGenerating fpkm.txt...")
    # fpkm_output_file = os.path.join(output_dir, "fpkm.txt")

    # try:
    #     # Create FPKM matrix
    #     fpkm_data = combined_df[["gene_id", "FPKM", "sample"]].copy()

    #     # Handle possible duplicates
    #     fpkm_data_grouped = (
    #         fpkm_data.groupby(["gene_id", "sample"])["FPKM"].sum().reset_index()
    #     )

    #     # Create pivot table
    #     fpkm_matrix = fpkm_data_grouped.pivot(
    #         index="gene_id", columns="sample", values="FPKM"
    #     )
    #     fpkm_matrix = fpkm_matrix.fillna(0)

    #     # Reset index to make gene_id the first column
    #     fpkm_matrix.reset_index(inplace=True)

    #     # Save file (without column names)
    #     fpkm_matrix.to_csv(fpkm_output_file, sep="\t", index=False, header=False)
    #     print(
    #         f"✓ Save completed: {fpkm_output_file} ({fpkm_matrix.shape[0]} genes × {fpkm_matrix.shape[1] - 1} samples)"
    #     )

    # except Exception as e:
    #     print(f"Error generating FPKM matrix: {e}")

    # # 3. Generate tpm.txt (TPM matrix)
    # print("\nGenerating tpm.txt...")
    # tpm_output_file = os.path.join(output_dir, "tpm.txt")

    # try:
    #     # Create TPM matrix
    #     tpm_data = combined_df[["gene_id", "TPM", "sample"]].copy()

    #     # Handle possible duplicates
    #     tpm_data_grouped = (
    #         tpm_data.groupby(["gene_id", "sample"])["TPM"].sum().reset_index()
    #     )

    #     # Create pivot table
    #     tpm_matrix = tpm_data_grouped.pivot(
    #         index="gene_id", columns="sample", values="TPM"
    #     )
    #     tpm_matrix = tpm_matrix.fillna(0)

    #     # Reset index to make gene_id the first column
    #     tpm_matrix.reset_index(inplace=True)

    #     # Save file (without column names)
    #     tpm_matrix.to_csv(tpm_output_file, sep="\t", index=False, header=False)
    #     print(
    #         f"✓ Save completed: {tpm_output_file} ({tpm_matrix.shape[0]} genes × {tpm_matrix.shape[1] - 1} samples)"
    #     )

    # except Exception as e:
    #     print(f"Error generating TPM matrix: {e}")

    print(f"\n{'=' * 60}")
    print("Expression matrix merging completed!")
    print(f"{'=' * 60}")
    print(f"Output files:")
    print(f"  - all.fpkm.tpm.txt: Direct merge of all data")
    # print(f"  - fpkm.txt: FPKM expression matrix")
    # print(f"  - tpm.txt: TPM expression matrix")


def parse_fastq_pattern(pattern):
    """Parse fastq file pattern to extract prefix, suffix and read identifiers"""
    if "*" not in pattern:
        raise ValueError("File pattern must contain * as sample name placeholder")

    # Split pattern string
    parts = pattern.split("*")
    if len(parts) != 2:
        raise ValueError("File pattern can only contain one * placeholder")

    prefix = parts[0]  # Part before *
    suffix = parts[1]  # Part after *

    # Detect read indicators (R1, R2 or 1, 2)
    read_indicators = ["R1", "R2", "_1", "_2", ".1", ".2"]
    read1_indicator = None
    read2_indicator = None

    for indicator in read_indicators:
        if indicator in suffix:
            if indicator.endswith("1") or indicator == "R1":
                read1_indicator = indicator
                if indicator == "R1":
                    read2_indicator = "R2"
                elif indicator == "_1":
                    read2_indicator = "_2"
                elif indicator == ".1":
                    read2_indicator = ".2"
            break

    if not read1_indicator:
        raise ValueError(
            f"Cannot identify read indicator (R1/R2, _1/_2, .1/.2) in pattern '{pattern}'"
        )

    return prefix, suffix, read1_indicator, read2_indicator


def parse_input_samples(input_path, fastq_pattern=None):
    """Parse input sample information"""
    samples = []

    if os.path.isdir(input_path):
        if fastq_pattern:
            # Use user-specified file pattern
            prefix, suffix, read1_indicator, read2_indicator = parse_fastq_pattern(
                fastq_pattern
            )

            # Build search pattern
            search_pattern = os.path.join(input_path, f"{prefix}*{suffix}")
            fastq_files = glob.glob(search_pattern)

            # Filter R1 files
            r1_files = [
                f for f in fastq_files if read1_indicator in os.path.basename(f)
            ]

            for fq1 in r1_files:
                # Build corresponding R2 file path
                fq2 = fq1.replace(read1_indicator, read2_indicator)

                if os.path.exists(fq2):
                    # Extract sample name
                    basename = os.path.basename(fq1)
                    sample_name = basename.replace(prefix, "").replace(suffix, "")

                    samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})
                else:
                    print(f"Warning: Cannot find corresponding R2 file: {fq2}")
        else:
            # Default mode: automatically search for common fastq file pairs
            patterns = [
                ("*_1.fq.gz", "*_2.fq.gz"),
                ("*_R1.fq.gz", "*_R2.fq.gz"),
                ("*.R1.fastq.gz", "*.R2.fastq.gz"),
                ("*_1.fastq.gz", "*_2.fastq.gz"),
                ("*.1.fq.gz", "*.2.fq.gz"),
                ("*_R1.fq", "*_R2.fq"),
                ("*_1.fastq", "*_2.fastq"),
            ]

            for pattern1, pattern2 in patterns:
                search_path1 = os.path.join(input_path, pattern1)
                fastq_files = glob.glob(search_path1)

                for fq1 in fastq_files:
                    # Build corresponding R2 file path
                    fq2 = fq1.replace(
                        pattern1.replace("*", ""), pattern2.replace("*", "")
                    )

                    if os.path.exists(fq2):
                        # Extract sample name
                        basename = os.path.basename(fq1)
                        # Remove fixed parts from pattern to get sample name
                        sample_name = basename
                        for to_remove in [
                            pattern1.replace("*", ""),
                            pattern2.replace("*", ""),
                        ]:
                            if to_remove in sample_name:
                                sample_name = sample_name.replace(to_remove, "")
                                break

                        # Ensure sample name is not empty
                        if not sample_name:
                            sample_name = os.path.splitext(
                                os.path.splitext(basename)[0]
                            )[0]

                        samples.append(
                            {"name": sample_name, "fastq1": fq1, "fastq2": fq2}
                        )

                # If files found, don't try other patterns
                if samples:
                    break
    else:
        # If it's a file, assume it's a sample information file
        # Format: sample_name\tfastq1_path\tfastq2_path
        with open(input_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    samples.append(
                        {"name": parts[0], "fastq1": parts[1], "fastq2": parts[2]}
                    )

    return samples


def main():
    parser = argparse.ArgumentParser(
        description="RNA-seq analysis pipeline: HISAT2 + StringTie"
    )

    # Required parameters
    parser.add_argument("-g", "--genome", required=True, help="Genome fasta file path")
    parser.add_argument(
        "-f", "--gtf", required=True, help="Gene annotation GTF file path"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input fastq file directory or sample information file",
    )
    parser.add_argument("-o", "--output", required=True, help="Output directory")

    # Optional parameters
    parser.add_argument(
        "-p",
        "--pattern",
        default=None,
        help='Fastq file naming pattern, e.g., "*.R1.fastq.gz" or "*_1.fq.gz", * represents sample name',
    )
    parser.add_argument(
        "-r",
        "--remove",
        default="no",
        choices=["yes", "y", "no", "n"],
        help="Remove BAM files after processing (default: no)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=8, help="Number of threads (default: 8)"
    )
    # parser.add_argument(
    #     "-k",
    #     "--keep_transcript",
    #     default="no",
    #     choices=["yes", "y", "no", "n"],
    #     help="Keep transcript or not (default: no)",
    # )

    args = parser.parse_args()

    # Check if required files exist
    if not os.path.exists(args.genome):
        print(f"Error: Genome file does not exist: {args.genome}")
        sys.exit(1)

    if not os.path.exists(args.gtf):
        print(f"Error: GTF file does not exist: {args.gtf}")
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    print("RNA-seq analysis pipeline started")
    print(f"Genome file: {args.genome}")
    print(f"Input path: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"GTF file: {args.gtf}")
    print(f"Threads: {args.threads}")
    print(f"Remove BAM files: {args.remove}")
    if args.pattern:
        print(f"File pattern: {args.pattern}")

    # 1. Build HISAT2 index
    print(f"\n{'=' * 60}")
    print("Step 1: Building HISAT2 index")
    print(f"{'=' * 60}")
    index_prefix = build_hisat2_index(args.genome, args.threads)

    # 2. Parse input samples
    print(f"\n{'=' * 60}")
    print("Step 2: Parsing input samples")
    print(f"{'=' * 60}")
    samples = parse_input_samples(args.input, args.pattern)

    if not samples:
        print("Error: No valid sample files found")
        sys.exit(1)

    print(f"Found {len(samples)} samples:")
    for sample in samples:
        print(f"  - {sample['name']}: {sample['fastq1']}, {sample['fastq2']}")

    # 3. Process all samples
    print(f"\n{'=' * 60}")
    print("Step 3: Processing all samples")
    print(f"{'=' * 60}")
    fpkm_files = []

    for sample_info in samples:
        fpkm_file = process_single_sample(args, sample_info, index_prefix)
        fpkm_files.append(fpkm_file)

    # 4. Merge expression matrix
    print(f"\n{'=' * 60}")
    print("Step 4: Merging expression matrix")
    print(f"{'=' * 60}")
    merge_expression_matrix(fpkm_files, args.output)

    print(f"\n{'=' * 60}")
    print("Analysis completed!")
    print(f"{'=' * 60}")
    print(f"Output files in: {args.output}")
    print("  - all.fpkm.tpm.txt: FPKM and TPM matrix for all samples.")
    print("  - fpkm.txt: FPKM only matrix for all samples.")
    print("  - tpm.txt: TPM only matrix for all samples.")


if __name__ == "__main__":
    main()
