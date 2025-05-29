#!/usr/bin/env python3
"""
VCF file chromosome SNP and INDEL statistics script with per-sample analysis
"""

import argparse
import sys
from collections import defaultdict

def classify_variant(ref, alt):
    """
    Classify variant type based on REF and ALT sequences
    
    Args:
        ref: Reference sequence
        alt: Alternative sequence
    
    Returns:
        str: 'SNP', 'INDEL', or 'MNV'
    """
    # Handle multiple alt alleles, take the first one
    if ',' in alt:
        alt = alt.split(',')[0]
    
    # SNP: same length and single base
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    # INDEL: different lengths
    elif len(ref) != len(alt):
        return 'INDEL'
    # MNV (multi-nucleotide variant): same length but more than one base
    else:
        return 'MNV'

def is_variant_in_sample(genotype):
    """
    Check if a sample has the variant based on genotype
    
    Args:
        genotype: Genotype string (e.g., "0/1", "1/1", "0|1")
    
    Returns:
        bool: True if sample has variant
    """
    if not genotype or genotype in ['.', './.', '.|.']:
        return False
    
    # Extract genotype calls
    gt = genotype.split(':')[0]  # Take only the GT field
    
    # Handle phased (|) and unphased (/) genotypes
    if '/' in gt:
        alleles = gt.split('/')
    elif '|' in gt:
        alleles = gt.split('|')
    else:
        return False
    
    # Check if any allele is non-reference (not 0)
    for allele in alleles:
        if allele != '0' and allele != '.':
            return True
    
    return False

def parse_vcf(vcf_file):
    """
    Parse VCF file and count variants per chromosome for total and per sample
    
    Args:
        vcf_file: VCF file path
    
    Returns:
        tuple: (total_stats, sample_stats, chromosome_lengths, sample_names)
    """
    total_stats = defaultdict(lambda: {'SNP': 0, 'INDEL': 0, 'MNV': 0, 'TOTAL': 0})
    sample_stats = defaultdict(lambda: defaultdict(lambda: {'SNP': 0, 'INDEL': 0, 'MNV': 0, 'TOTAL': 0}))
    chromosome_lengths = {}
    sample_names = []
    
    try:
        with open(vcf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Process header lines, extract chromosome length information
                if line.startswith('##contig='):
                    try:
                        contig_info = line[line.find('<')+1:line.find('>')]
                        parts = contig_info.split(',')
                        chrom_id = None
                        chrom_length = None
                        
                        for part in parts:
                            if part.startswith('ID='):
                                chrom_id = part.split('=')[1]
                            elif part.startswith('length='):
                                chrom_length = int(part.split('=')[1])
                        
                        if chrom_id and chrom_length:
                            chromosome_lengths[chrom_id] = chrom_length
                    except Exception as e:
                        print(f"Warning: Failed to parse contig info at line {line_num}: {e}", file=sys.stderr)
                
                # Process column header line to get sample names
                elif line.startswith('#CHROM'):
                    fields = line.split('\t')
                    if len(fields) > 9:
                        sample_names = fields[9:]  # Sample names start from column 10
                    continue
                
                # Skip other comment lines
                elif line.startswith('#'):
                    continue
                
                # Skip empty lines
                if not line:
                    continue
                
                # Parse VCF data line
                fields = line.split('\t')
                if len(fields) < 5:
                    print(f"Warning: Malformed line {line_num}, skipping", file=sys.stderr)
                    continue
                
                chrom = fields[0]  # Chromosome
                pos = fields[1]    # Position
                ref = fields[3]    # Reference sequence
                alt = fields[4]    # Alternative sequence
                
                # Skip invalid variants
                if alt == '.' or ref == '.':
                    continue
                
                # Classify variant type
                variant_type = classify_variant(ref, alt)
                
                # Count total variants
                total_stats[chrom][variant_type] += 1
                total_stats[chrom]['TOTAL'] += 1
                
                # Count per-sample variants
                if len(fields) > 9 and sample_names:
                    for i, sample_name in enumerate(sample_names):
                        sample_idx = i + 9
                        if sample_idx < len(fields):
                            genotype = fields[sample_idx]
                            if is_variant_in_sample(genotype):
                                sample_stats[sample_name][chrom][variant_type] += 1
                                sample_stats[sample_name][chrom]['TOTAL'] += 1
                
                # If no length info from header, infer minimum length from variant positions
                if chrom not in chromosome_lengths:
                    try:
                        pos_int = int(pos)
                        if chrom not in chromosome_lengths or pos_int > chromosome_lengths.get(chrom, 0):
                            chromosome_lengths[chrom] = pos_int
                    except ValueError:
                        pass
    
    except FileNotFoundError:
        print(f"Error: File not found: {vcf_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Problem reading file - {e}", file=sys.stderr)
        sys.exit(1)
    
    return total_stats, sample_stats, chromosome_lengths, sample_names

def print_total_statistics_clean(stats, lengths):
    """
    Print total statistics in clean tab-separated format
    
    Args:
        stats: Total chromosome statistics dictionary
        lengths: Chromosome lengths dictionary
    """
    print("# TOTAL VARIANT STATISTICS")
    print("Chromosome\tLength\tSNP\tINDEL\tMNV\tTotal")
    
    # Sort chromosomes by name
    sorted_chroms = sorted(stats.keys(), key=lambda x: (
        0 if x.isdigit() else 1,  # Numeric chromosomes first
        int(x) if x.isdigit() else x  # Numeric by value, others alphabetically
    ))
    
    total_snp = 0
    total_indel = 0
    total_mnv = 0
    total_variants = 0
    total_length = 0
    
    for chrom in sorted_chroms:
        snp_count = stats[chrom]['SNP']
        indel_count = stats[chrom]['INDEL']
        mnv_count = stats[chrom]['MNV']
        total_count = stats[chrom]['TOTAL']
        
        # Get chromosome length
        chrom_length = lengths.get(chrom, 'NA')
        if chrom_length != 'NA':
            total_length += chrom_length
        
        print(f"{chrom}\t{chrom_length}\t{snp_count}\t{indel_count}\t{mnv_count}\t{total_count}")
        
        total_snp += snp_count
        total_indel += indel_count
        total_mnv += mnv_count
        total_variants += total_count
    
    # Print totals
    total_length_val = total_length if total_length > 0 else 'NA'
    print(f"TOTAL\t{total_length_val}\t{total_snp}\t{total_indel}\t{total_mnv}\t{total_variants}")
    
    # Additional statistics as comments
    if total_length > 0:
        print(f"# Total genome length: {total_length} bp")
        print(f"# Variant density: {total_variants/total_length*1000:.4f} variants/kb")
        print(f"# SNP density: {total_snp/total_length*1000:.4f} SNPs/kb")
        print(f"# INDEL density: {total_indel/total_length*1000:.4f} INDELs/kb")

def print_sample_statistics_by_chromosome_clean(sample_stats, lengths, sample_names, total_stats):
    """
    Print per-chromosome statistics with all samples in clean format
    
    Args:
        sample_stats: Per-sample chromosome statistics dictionary
        lengths: Chromosome lengths dictionary
        sample_names: List of sample names
        total_stats: Total chromosome statistics dictionary
    """
    if not sample_names:
        print("\n# No sample information found in VCF file")
        return
    
    print("\n# PER-CHROMOSOME VARIANT STATISTICS")
    
    # Header
    header = "Chromosome\tLength\tType\t" + "\t".join(sample_names) + "\tSample_Total\tVCF_Total"
    print(header)
    
    # Get all chromosomes
    all_chroms = set()
    all_chroms.update(total_stats.keys())
    for sample_data in sample_stats.values():
        all_chroms.update(sample_data.keys())
    
    # Sort chromosomes by name
    sorted_chroms = sorted(all_chroms, key=lambda x: (
        0 if x.isdigit() else 1,
        int(x) if x.isdigit() else x
    ))
    
    for chrom in sorted_chroms:
        chrom_length = lengths.get(chrom, 'NA')
        
        # Print data for each variant type
        for variant_type in ['SNP', 'INDEL', 'MNV', 'TOTAL']:
            sample_counts = []
            sample_total = 0
            
            # Get counts for each sample
            for sample_name in sample_names:
                sample_data = sample_stats[sample_name]
                if chrom in sample_data:
                    count = sample_data[chrom][variant_type]
                else:
                    count = 0
                sample_counts.append(str(count))
                sample_total += count
            
            # Get VCF total for this chromosome and variant type
            if chrom in total_stats:
                vcf_total = total_stats[chrom][variant_type]
            else:
                vcf_total = 0
            
            # Print the line
            line_parts = [chrom, str(chrom_length), variant_type] + sample_counts + [str(sample_total), str(vcf_total)]
            print("\t".join(line_parts))
    
    # Print density information as comments
    print("\n# VARIANT DENSITY INFORMATION")
    print("Chromosome\tLength\tSample_Density_per_kb\tVCF_Density_per_kb")
    
    for chrom in sorted_chroms:
        chrom_length = lengths.get(chrom, 'NA')
        
        if chrom_length != 'NA':
            # Calculate sample total for this chromosome
            sample_total = 0
            for sample_name in sample_names:
                sample_data = sample_stats[sample_name]
                if chrom in sample_data:
                    sample_total += sample_data[chrom]['TOTAL']
            
            # Calculate VCF total for this chromosome
            vcf_total = total_stats[chrom]['TOTAL'] if chrom in total_stats else 0
            
            # Calculate densities
            sample_density = sample_total / chrom_length * 1000
            vcf_density = vcf_total / chrom_length * 1000
            
            print(f"{chrom}\t{chrom_length}\t{sample_density:.4f}\t{vcf_density:.4f}")
        else:
            print(f"{chrom}\tNA\tNA\tNA")

def print_sample_summary_clean(sample_stats, sample_names):
    """
    Print summary table of all samples in clean format
    
    Args:
        sample_stats: Per-sample chromosome statistics dictionary
        sample_names: List of sample names
    """
    if not sample_names:
        return
    
    print("\n# SAMPLE SUMMARY")
    print("Sample\tSNP\tINDEL\tMNV\tTotal")
    
    for sample_name in sample_names:
        sample_data = sample_stats[sample_name]
        
        total_snp = sum(chrom_data['SNP'] for chrom_data in sample_data.values())
        total_indel = sum(chrom_data['INDEL'] for chrom_data in sample_data.values())
        total_mnv = sum(chrom_data['MNV'] for chrom_data in sample_data.values())
        total_variants = sum(chrom_data['TOTAL'] for chrom_data in sample_data.values())
        
        print(f"{sample_name}\t{total_snp}\t{total_indel}\t{total_mnv}\t{total_variants}")

def main():
    parser = argparse.ArgumentParser(
        description='Count SNPs and INDELs per chromosome in VCF files (total and per-sample)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.vcf
  %(prog)s sample.vcf -o output.txt
  %(prog)s variants.vcf --output results.txt --sample-details

Notes:
  - SNP: Single nucleotide polymorphism (REF and ALT both length 1)
  - INDEL: Insertion/deletion variant (REF and ALT different lengths)
  - MNV: Multi-nucleotide variant (REF and ALT same length but >1)
  - Per-sample analysis requires sample columns in VCF (columns 10+)
  - Only variants present in each sample (non-reference genotypes) are counted
  - Compressed VCF files (.vcf.gz) need to be decompressed first
        """
    )
    
    parser.add_argument('vcf_file', 
                       help='Input VCF file path')
    parser.add_argument('-o', '--output', 
                       help='Output file path (optional, prints to stdout if not specified)')
    parser.add_argument('-v', '--verbose', 
                       action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--sample-details', 
                       action='store_true',
                       help='Show detailed per-sample statistics for each chromosome')
    parser.add_argument('--total-only', 
                       action='store_true',
                       help='Show only total statistics, skip per-sample analysis')
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Analyzing VCF file: {args.vcf_file}")
        print("Processing...")
    
    # Parse VCF file
    total_stats, sample_stats, lengths, sample_names = parse_vcf(args.vcf_file)
    
    if not total_stats:
        print("Warning: No valid variant records found", file=sys.stderr)
        return
    
    # Output results
    if args.output:
        # Redirect output to file
        with open(args.output, 'w') as f:
            import contextlib
            with contextlib.redirect_stdout(f):
                print_total_statistics_clean(total_stats, lengths)
                
                if not args.total_only and sample_names:
                    print_sample_summary_clean(sample_stats, sample_names)
                    
                    if args.sample_details:
                        print_sample_statistics_by_chromosome_clean(sample_stats, lengths, sample_names, total_stats)
        
        if args.verbose:
            print(f"Results saved to: {args.output}")
    else:
        print_total_statistics_clean(total_stats, lengths)
        
        if not args.total_only and sample_names:
            print_sample_summary_clean(sample_stats, sample_names)
            
            if args.sample_details:
                print_sample_statistics_by_chromosome_clean(sample_stats, lengths, sample_names, total_stats)

if __name__ == "__main__":
    main()