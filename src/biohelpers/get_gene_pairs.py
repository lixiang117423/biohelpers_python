import argparse
from typing import Set, Dict, List, Tuple
import gffutils
from collections import defaultdict


def natural_chrom_key(chrom):
    try:
        num_part = ''.join(filter(str.isdigit, chrom))
        return int(num_part) if num_part else float('inf')
    except:
        return float('inf')


def build_gene_index(gff_path: str) -> Dict[str, List[Tuple]]:
    db = gffutils.create_db(gff_path, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge")
    chrom_index = defaultdict(list)
    gene_dict = {}

    # 提取基因和mRNA信息
    for feat in db.features_of_type(['gene', 'mRNA']):
        if feat.featuretype == 'gene':
            gene_id = feat.id
            chrom = feat.chrom
            start = feat.start
            end = feat.end
            strand = feat.strand
            gene_dict[gene_id] = {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcripts': []
    }
        elif feat.featuretype == 'mRNA':
            parent_gene = feat.attributes.get('Parent', [None])[0]
            if parent_gene in gene_dict:
                gene_dict[parent_gene]['transcripts'].append(feat.id)
                chrom = gene_dict[parent_gene]['chrom']
                chrom_index[chrom].append(parent_gene)

    # 自然排序染色体名称
    sorted_chroms = sorted(chrom_index.keys(), key=natural_chrom_key)

    # 按染色体自然顺序和基因起始位置排序
    sorted_index = {}
    for chrom in sorted_chroms:
        genes = chrom_index[chrom]
        sorted_genes = sorted(genes, key=lambda x: gene_dict[x]['start'])
        sorted_index[chrom] = {gene: i+1 for i, gene in enumerate(sorted_genes)}

    return sorted_index, gene_dict


def find_gene_pairs(target_genes: Set[str], index: Dict, distance: int) -> Set[Tuple]:
    pairs = set()
    for gene in target_genes:
        for chrom, genes in index.items():
            if gene in genes:
                idx = genes[gene]
                start = max(1, idx - distance)
                end = idx + distance
                
                # 获取窗口内的基因
                window_genes = [g for g, i in genes.items() if start <= i <= end]
                
                # 生成唯一配对
                for other in window_genes:
                    if other in target_genes and other != gene:
                        pair = tuple(sorted((gene, other)))
                        pairs.add(pair)
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Find gene pairs near target genes")
    parser.add_argument("--gff", "-g", required=True, help="Path to GFF file")
    parser.add_argument("--id", "-i", required=True, help="File containing target gene IDs")
    parser.add_argument("--type", "-t", choices=["gene", "mrna"], default="gene")
    parser.add_argument("--distance", "-d", type=int, default=3)
    parser.add_argument("--output", "-o", required=True, help="Output filename")
    args = parser.parse_args()

    # 生成染色体索引
    gene_index, gene_info = build_gene_index(args.gff)
    
    # 处理目标基因
    with open(args.id) as f:
        target_ids = set(line.strip() for line in f)
    
    # 转换mRNA到gene
    if args.type == "mrna":
        target_genes = set()
        mrna_to_gene = {mrna: gene for gene in gene_info for mrna in gene_info[gene]['transcripts']}
        for mrna in target_ids:
            if mrna in mrna_to_gene:
                target_genes.add(mrna_to_gene[mrna])
    else:
        target_genes = target_ids

    # 查找基因对
    pairs = find_gene_pairs(target_genes, gene_index, args.distance)

    # 写入结果
    with open(args.output, 'w') as f:
        f.write("Chr\tGene.1\tGene.2\tChr\tStart.1\tEnd.1\tStrand.1\tStart.2\tEnd.2\tStrand.2\n")
        # 按染色体自然顺序和起始位置排序
        def sort_key(pair):
            chrom = gene_info[pair[0]]['chrom']
            return (natural_chrom_key(chrom), gene_info[pair[0]]['start'])

        sorted_pairs = sorted(pairs, key=sort_key)
        for gene1, gene2 in sorted_pairs:
            chrom1 = gene_info[gene1]['chrom']
            s1 = gene_info[gene1]['start']
            e1 = gene_info[gene1]['end']
            st1 = gene_info[gene1]['strand']
            s2 = gene_info[gene2]['start']
            e2 = gene_info[gene2]['end']
            st2 = gene_info[gene2]['strand']
            f.write(f"{chrom1}\t{gene1}\t{gene2}\t{chrom1}\t{s1}\t{e1}\t{st1}\t{s2}\t{e2}\t{st2}\n")

if __name__ == "__main__":
    main()