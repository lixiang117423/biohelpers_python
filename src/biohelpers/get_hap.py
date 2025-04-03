import argparse
import gzip
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Extract haplotype information from VCF files')
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF file path')
    parser.add_argument('-c', '--chr', required=True, help='Chromosome identifier')
    parser.add_argument('-p', '--position', type=int, required=True, help='Target SNP position')
    parser.add_argument('-s', '--start', type=int, default=0, help='Upstream window size')
    parser.add_argument('-e', '--end', type=int, default=0, help='Downstream window size')
    parser.add_argument('-t', '--type', choices=['sample', 'hap'], required=True, help='Output type')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    return parser.parse_args()

def validate_args(args):
    if args.start < 0 or args.end < 0:
        raise ValueError('Window size cannot be negative')
    # 移除文件扩展名验证以支持更多格式

def process_genotype(gt, ref, alts):
    if gt in ('./.', '.|.'):
        return ('./.', './.', 'Missing', 'No phasing info')
    
    alleles = []
    separator = '|' if '|' in gt else '/'
    phase_info = 'Phased' if separator == '|' else 'Unphased'
    
    for code in gt.split(separator):
        if code == '.':
            return (gt, './.', 'Missing', 'No phasing info')
        try:
            idx = int(code)
            alleles.append(ref if idx == 0 else alts[idx-1])
        except (IndexError, ValueError):
            return (gt, './.', 'Missing', 'No phasing info')
    
    if len(set(alleles)) == 1:
        biological_meaning = 'Homozygous Reference' if alleles[0] == ref else 'Homozygous Alternate'
    else:
        biological_meaning = 'Heterozygous'
    
    base_combination = '/'.join(alleles)
    return (gt, base_combination, biological_meaning)

def parse_vcf(vcf_path, target_chr, start_pos, end_pos):
    samples = []
    hap_counts = defaultdict(int)
    # 通过魔数检测文件类型
    with open(vcf_path, 'rb') as test_f:
        header = test_f.read(2)
    
    if header == b'\x1f\x8b':
        open_func = gzip.open
    else:
        open_func = open
    
    try:
        with open_func(vcf_path, 'rt' if isinstance(open_func, type(open)) else 'rb') as f:
            sample_data = []
            for line in f:
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[9:]
                    continue
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                chrom, pos = fields[0], int(fields[1])
                
                if chrom != target_chr or not (start_pos <= pos <= end_pos):
                    continue
                
                ref = fields[3]
                alts = fields[4].split(',')
                for sample, gt in zip(samples, fields[9:]):
                    processed_gt = process_genotype(gt.split(':')[0], ref, alts)
                    sample_data.append((sample, processed_gt, ref, alts, chrom, pos))
                    hap_counts[processed_gt] += 1
            return sample_data, hap_counts
    except IOError as e:
        raise ValueError(f'Failed to open VCF file: {str(e)}')
    
def write_output(output_type, data, output_path, hap_counts):
    with open(output_path, 'w') as f:
        if output_type == 'sample':
            f.write('Chr\tPosition\tREF\tALT\tSample\tGT\tAlleles\tFrequency\tBiological_Meaning\n')
            total = sum(hap_counts.values()) or 1
            for sample, genotype_info, ref, alts, chrom, pos in data:
                freq = hap_counts[genotype_info] / total
                line = [
                    chrom, pos, ref, ",".join(alts), sample, genotype_info[0], genotype_info[1], f"{freq:.2%}", genotype_info[2]
                ]
                f.write('\t'.join(map(str, line)) + '\n')
        else:
            f.write('hap\tallel\tnumber\n')
            for gt, count in data.items():
                f.write(f'{gt}\t{count}\n')

def main():
    args = parse_args()
    validate_args(args)
    
    start_window = args.position - args.start
    end_window = args.position + args.end
    
    sample_data, hap_counts = parse_vcf(args.vcf, args.chr, start_window, end_window)
    
    if args.type == 'sample':
        write_output(args.type, sample_data, args.output, hap_counts)
    else:
        write_output(args.type, hap_counts, args.output, hap_counts)

if __name__ == '__main__':
    main()