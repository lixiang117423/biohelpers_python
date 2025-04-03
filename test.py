#!/usr/bin/env python3
import argparse
import pysam
from collections import defaultdict

def parse_genotype(gt, ref, alts):
    """解析基因型为等位碱基"""
    alleles = []
    for allele in gt:
        if allele is None:
            alleles.append("Missing")
        elif allele == 0:
            alleles.append(ref)
        elif allele <= len(alts):
            alleles.append(alts[allele-1])
        else:
            alleles.append("Missing")
    return alleles

def main():
    parser = argparse.ArgumentParser(description='Extract haplotypes from VCF')
    parser.add_argument('--vcf', '-v', required=True, help='Input VCF file path')
    parser.add_argument('--chr', '-c', required=True, help='Chromosome name')
    parser.add_argument('--position', '-p', type=int, required=True, help='Central SNP position')
    parser.add_argument('--start', '-s', type=int, default=0, help='Upstream window size')
    parser.add_argument('--end', '-e', type=int, default=0, help='Downstream window size')
    parser.add_argument('--type', '-t', choices=['sample', 'hap'], required=True, help='Output type')
    parser.add_argument('--output', '-o', required=True, help='Output file path')
    args = parser.parse_args()

    # 计算查询区域
    region_start = max(1, args.position - args.start)
    region_end = args.position + args.end

    # 读取VCF文件 
    vcf = pysam.VariantFile(args.vcf)
    samples = list(vcf.header.samples)
    hap_counter = defaultdict(int)
    sample_data = []

    try:
        # 获取目标区域的所有变异记录
        for record in vcf.fetch(args.chr, region_start-1, region_end):
            ref = record.ref
            alts = record.alts if record.alts else []
            
            # 遍历所有样本
            for sample in samples:
                gt = record.samples[sample].get('GT', (None, None))
                alleles = parse_genotype(gt, ref, alts)
                
                # 根据输出类型处理数据
                if args.type == "sample":
                    sample_allel = "/".join(alleles) if "Missing" not in alleles else "./."
                    sample_data.append((sample, sample_allel))
                else:
                    for allele in alleles:
                        hap_counter[allele] += 1
    except ValueError as e:
        print(f"Error: {str(e)}")
        exit(1)

    # 写入输出文件
    with open(args.output, 'w') as f:
        if args.type == "sample":
            f.write("sample\tallel\n")
            for sample, allel in sample_data:
                f.write(f"{sample}\t{allel}\n")
        else:
            sorted_hap = sorted(hap_counter.items(), key=lambda x: (-x[1], x[0]))
            f.write("hap\tallel\tnumber\n")
            for idx, (allel, count) in enumerate(sorted_hap, 1):
                f.write(f"{idx}\t{allel}\t{count}\n")

if __name__ == "__main__":
    main()
