import argparse
import gzip
import subprocess
import os
import tempfile
from collections import defaultdict

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
    # 移除文件扩展名验证以支持更多格式

def process_genotype(gt, ref, alts):
    if gt in ('./.', '.|.'):
        return ('./.', './.', 'Missing', 'No phasing info')
    
    alleles = []
    separator = '|' if '|' in gt else '/'
    # phase_info = 'Phased' if separator == '|' else 'Unphased'
    
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
    temp_dir = tempfile.mkdtemp()
    intermediate_vcf = os.path.normpath(os.path.join(temp_dir, 'filtered.vcf.gz'))

    try:
        if subprocess.run(['bcftools', '--version'], capture_output=True).returncode != 0:
            raise EnvironmentError('请先安装并配置bcftools环境')

        cmd = [
            'bcftools', 'view',
            '-r', f'{target_chr}:{start_pos}-{end_pos}',
            '-Oz', '-o', intermediate_vcf,
            os.path.normpath(vcf_path)
        ]
        subprocess.run(cmd, check=True)

        # 使用预处理后的文件继续处理
        samples = []
        hap_counts = defaultdict(int)
        sample_data = []
        with gzip.open(intermediate_vcf, 'rt') as f:
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
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f'bcftools执行失败: {e.stderr.decode()}') from e
    finally:
        if os.path.exists(intermediate_vcf):
            os.remove(intermediate_vcf)
        os.rmdir(temp_dir)
def write_output(data, output_path, hap_counts):
    with open(output_path, 'w') as f:
        f.write('Chr\tPosition\tREF\tALT\tSample\tGT\tAlleles\tFrequency\tBiological_Meaning\n')
        total = sum(hap_counts.values()) or 1
        for sample, genotype_info, ref, alts, chrom, pos in data:
            freq = hap_counts[genotype_info] / total
            line = [
                chrom, pos, ref, ",".join(alts), sample, genotype_info[0], genotype_info[1], f"{freq:.2%}", genotype_info[2]
            ]
            f.write('\t'.join(map(str, line)) + '\n')

def format_console_output(data):
    unique_stats = defaultdict(lambda: {'count':0, 'samples':set()})
    total_samples = len({entry[0] for entry in data})
    
    for entry in data:
        key = (
            entry[4],  # chrom
            entry[5],  # pos
            entry[2],  # ref
            ','.join(entry[3]),  # alts
            entry[1][0],  # GT
            entry[1][1]  # alleles
        )
        unique_stats[key]['count'] += 1
        unique_stats[key]['samples'].add(entry[0])
    
    print('\nConsolidated Haplotype Statistics:')
    print('Chr\tPosition\tREF\tALT\tGT\tAlleles\tUniqueSamples\tFrequency')
    for key in sorted(unique_stats.keys()):
        stats = unique_stats[key]
        freq = stats['count'] / len(data)
        print(f'{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{key[4]}\t{key[5]}\t{len(stats["samples"])}\t{freq:.2%}')


def main():
    args = parse_args()
    validate_args(args)
    start_window = args.position - args.start
    end_window = args.position + args.end
    sample_data, hap_counts = parse_vcf(args.vcf, args.chr, start_window, end_window)
    write_output(sample_data, args.output, hap_counts)
    format_console_output(sample_data)

if __name__ == '__main__':
    main()