import argparse
import os
import gzip
import io
from collections import defaultdict
from tabulate import tabulate

def parse_args():
    parser = argparse.ArgumentParser(description='Extract haplotype information from VCF files')
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF file path')
    parser.add_argument('-c', '--chr', required=True, help='Chromosome identifier')
    parser.add_argument('-p', '--position', type=int, required=True, help='Target SNP position')
    parser.add_argument('-s', '--start', type=int, default=0, help='Upstream window size')
    parser.add_argument('-e', '--end', type=int, default=0, help='Downstream window size')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    return parser.parse_args()

def validate_args(args):
    if args.start < 0 or args.end < 0:
        raise ValueError('Window size cannot be negative')
    if not args.output:
        raise ValueError('Output path cannot be empty')
    if not os.path.exists(args.vcf):
        raise FileNotFoundError(f'VCF file not found: {args.vcf}')
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

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
    import pandas as pd
    
    try:
        # 自动检测压缩格式并读取
        # 获取样本名称
        with gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[9:]
                    if not samples:
                        raise ValueError('No samples found in VCF header')
                    break
        
        # 读取数据并保留所有样本列
        # 加载数据并处理压缩格式
        df = pd.read_csv(
            vcf_path,
            compression='gzip' if vcf_path.endswith('.gz') else None,
            comment='#',
            sep='\t',
            header=None,
            usecols=[0,1,3,4] + list(range(9, 9 + len(samples))),
            low_memory=False,
            dtype={'CHROM':'str', 'POS':'int32', 'REF':'str', 'ALT':'str'},
            names=['CHROM','POS','REF','ALT'] + samples
        )
        
        # 向量化筛选并转换数据格式
        filtered = df[
            (df.CHROM == target_chr) & 
            (df.POS.between(start_pos, end_pos))
        ].melt(
            id_vars=['CHROM','POS','REF','ALT'],
            value_vars=samples,
            var_name='Sample',
            value_name='GT'
        )

        # 批量处理基因型数据
        sample_data = []
        hap_counts = defaultdict(int)
        for _, row in filtered.iterrows():
            gt = row['GT'].split(':')[0]
            processed = process_genotype(gt, row['REF'], row['ALT'].split(','))
            sample_data.append((
                row['Sample'],
                processed,
                row['REF'],
                row['ALT'].split(','),
                row['CHROM'],
                row['POS']
            ))
            hap_counts[processed] += 1
        
        return sample_data, hap_counts
        
    except Exception as e:
        raise RuntimeError(f'VCF processing error: {str(e)}') from e

# Removed legacy line-by-line processing function
    
def write_output(data, output_path, hap_counts):
    total = sum(hap_counts.values()) or 1
    output_data = []
    with open(output_path, 'w') as f:
        f.write('Chr\tPosition\tREF\tALT\tSample\tGT\tAlleles\tFrequency\tBiological_Meaning\n')
        for sample, genotype_info, ref, alts, chrom, pos in data:
            freq = hap_counts[genotype_info] / total
            line = [
                chrom, pos, ref, ",".join(alts), sample, genotype_info[0], genotype_info[1], f"{freq:.2%}", genotype_info[2]
            ]
            f.write('\t'.join(map(str, line)) + '\n')
            output_data.append((chrom, pos, ref, ",".join(alts)) + genotype_info + (freq,))
    return output_data, hap_counts, total

def format_console_output(hap_counts, total):
    headers = ['Chr', 'Position', 'REF', 'ALT', 'GT', 'Alleles', 'Frequency', 'Biological_Meaning']
    stats = defaultdict(list)
    
    for genotype_info, count in hap_counts.items():
        freq = count / total * 100
        key = (genotype_info[0], genotype_info[1], genotype_info[3])
        stats[key].append((freq, count))
    
    lines = []
    for (gt, alleles, bio_meaning), values in stats.items():
        total_freq = sum(f for f, _ in values)
        lines.append({
            'GT': gt,
            'Alleles': alleles,
            'Frequency': f"{total_freq:.1%}",
            'Biological_Meaning': bio_meaning
        })
    
    # 生成表格输出
    table = tabulate(
        [(item['GT'], item['Alleles'], item['Frequency'], item['Biological_Meaning']) for item in lines],
        headers=['Chr', 'Position', 'REF', 'ALT', 'GT', 'Alleles', 'Frequency', 'Biological_Meaning'],
        tablefmt='plain',
        numalign='center',
        stralign='center'
    )
    return table

def main():
    args = parse_args()
    validate_args(args)
    
    start_window = args.position - args.start
    end_window = args.position + args.end
    
    sample_data, hap_counts = parse_vcf(args.vcf, args.chr, start_window, end_window)
    
    output_data, hap_counts, total = write_output(sample_data, args.output, hap_counts)
    
    # 打印控制台统计信息
    print("\nVariant Statistics:")
    print(format_console_output(hap_counts, total))
    print(f"\nTotal samples processed: {total}")

if __name__ == '__main__':
    main()