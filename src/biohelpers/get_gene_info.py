#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
从 GFF3 文件中提取基因和转录本（mRNA）的整合信息。
对于每个转录本，输出一行，包含其自身的坐标和其父基因的坐标。
"""

import argparse
import csv
import sys
from collections import defaultdict

def parse_attributes(attr_string):
    """
    解析 GFF3 第九列的属性字符串。
    例如: "ID=gene1;Name=GeneA" -> {'ID': 'gene1', 'Name': 'GeneA'}
    """
    attributes = {}
    for part in attr_string.strip().split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            attributes[key.strip()] = value.strip()
    return attributes

def extract_gene_transcript_info(gff3_file, output_file):
    """
    主处理函数：
    1. 第一遍：扫描文件，将所有基因的信息存储在字典中。
    2. 第二遍：扫描文件，处理转录本，并结合已存储的基因信息写入输出文件。
    """
    # 定义感兴趣的特征类型
    transcript_types = {'mRNA', 'transcript'}
    
    # --- 第一遍：收集所有基因的信息 ---
    gene_data = {}
    try:
        with open(gff3_file, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) != 9:
                    continue
                
                feature_type = parts[2]
                if feature_type == 'gene':
                    attributes = parse_attributes(parts[8])
                    gene_id = attributes.get('ID')
                    if gene_id:
                        gene_data[gene_id] = {
                            'chr': parts[0],
                            'start': parts[3],
                            'end': parts[4],
                            'strand': parts[6]
                        }
    except FileNotFoundError:
        print(f"错误: 输入文件 '{gff3_file}' 未找到。", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"第一遍读取文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 第二遍：处理转录本并写入输出文件 ---
    try:
        with open(gff3_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')

            # 写入新的表头
            header = ['Gene_ID', 'Transcript_ID', 'Chromosome', 'Strand', 'Gene_Start', 'Gene_End', 'Transcript_Start', 'Transcript_End']
            writer.writerow(header)

            # 逐行读取GFF3文件
            for line in infile:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) != 9:
                    continue

                feature_type = parts[2]
                
                # 只处理转录本行
                if feature_type in transcript_types:
                    attributes = parse_attributes(parts[8])
                    transcript_id = attributes.get('ID')
                    gene_id = attributes.get('Parent')

                    if not transcript_id or not gene_id:
                        continue
                    
                    # 获取转录本自己的信息
                    transcript_start = parts[3]
                    transcript_end = parts[4]

                    # 从之前收集的基因数据中查找父基因信息
                    parent_gene_info = gene_data.get(gene_id)
                    
                    if parent_gene_info:
                        # 如果找到了父基因
                        chromosome = parent_gene_info['chr']
                        strand = parent_gene_info['strand']
                        gene_start = parent_gene_info['start']
                        gene_end = parent_gene_info['end']
                    else:
                        # 如果没找到父基因（孤儿转录本），则基因信息用'NA'填充
                        # 染色体和链信息可以从转录本行本身获取
                        print(f"警告: 未找到转录本 '{transcript_id}' 的父基因 '{gene_id}'。基因坐标将标记为 NA。", file=sys.stderr)
                        chromosome = parts[0]
                        strand = parts[6]
                        gene_start = 'NA'
                        gene_end = 'NA'
                    
                    output_row = [gene_id, transcript_id, chromosome, strand, gene_start, gene_end, transcript_start, transcript_end]
                    writer.writerow(output_row)

    except Exception as e:
        print(f"第二遍处理文件并写入时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"信息提取完成。结果已保存到 '{output_file}'。")


def main():
    """
    解析命令行参数并调用主处理函数。
    """
    parser = argparse.ArgumentParser(
        description="从GFF3文件中为每个转录本提取整合的基因和转录本信息。",
        epilog="示例: python extract_gff_info_v2.py -g input.gff3 -o gene_transcript_info.tsv"
    )
    
    parser.add_argument(
        '--gff3', '-g',
        required=True,
        type=str,
        help="输入的GFF3文件路径。"
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        type=str,
        help="输出的TSV文件路径。"
    )

    args = parser.parse_args()

    extract_gene_transcript_info(args.gff3, args.output)


if __name__ == '__main__':
    main()