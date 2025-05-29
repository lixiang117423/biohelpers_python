import argparse

import pandas as pd
from natsort import index_natsorted

# def natural_chrom_key(chrom):
#     try:
#         num_part = ''.join(filter(str.isdigit, chrom))
#         return int(num_part) if num_part else float('inf')
#     except:
#         return float('inf')
# def build_gene_index(gff_path: str) -> Dict[str, List[Tuple]]:
#     db = gffutils.create_db(gff_path, dbfn=":memory:", force=True, keep_order=True, merge_strategy="merge")
#     chrom_index = defaultdict(list)
#     gene_dict = {}
#     # 提取基因和mRNA信息
#     for feat in db.features_of_type(['gene', 'mRNA']):
#         if feat.featuretype == 'gene':
#             gene_id = feat.id
#             chrom = feat.chrom
#             start = feat.start
#             end = feat.end
#             strand = feat.strand
#             gene_dict[gene_id] = {
#         'chrom': chrom,
#         'start': start,
#         'end': end,
#         'strand': strand,
#         'transcripts': []
#     }
#         elif feat.featuretype == 'mRNA':
#             parent_gene = feat.attributes.get('Parent', [None])[0]
#             if parent_gene in gene_dict:
#                 gene_dict[parent_gene]['transcripts'].append(feat.id)
#                 chrom = gene_dict[parent_gene]['chrom']
#                 chrom_index[chrom].append(parent_gene)
#     # 自然排序染色体名称
#     sorted_chroms = sorted(chrom_index.keys(), key=natural_chrom_key)
#     # 按染色体自然顺序和基因起始位置排序
#     sorted_index = {}
#     for chrom in sorted_chroms:
#         genes = chrom_index[chrom]
#         sorted_genes = sorted(genes, key=lambda x: gene_dict[x]['start'])
#         sorted_index[chrom] = {gene: i+1 for i, gene in enumerate(sorted_genes)}
#     return sorted_index, gene_dict
# def find_gene_pairs(target_genes: Set[str], index: Dict, distance: int) -> Set[Tuple]:
#     pairs = set()
#     for gene in target_genes:
#         for chrom, genes in index.items():
#             if gene in genes:
#                 idx = genes[gene]
#                 start = max(1, idx - distance)
#                 end = idx + distance
#                 # 获取窗口内的基因
#                 window_genes = [g for g, i in genes.items() if start <= i <= end]
#                 # 生成唯一配对
#                 for other in window_genes:
#                     if other in target_genes and other != gene:
#                         pair = tuple(sorted((gene, other)))
#                         pairs.add(pair)
#     return pairs


def parse_gff_features(gff_path):
    """
    Enhanced GFF3 parser that extracts both gene and mRNA information while establishing their relationships
    Improvements:
    1. Simultaneously parses gene and mRNA information
    2. Establishes gene-mRNA hierarchical relationships
    3. Enhanced attribute parsing compatibility
    4. Added cross-feature validation
    """
    # Pre-allocate dual cache structure
    gene_records = []
    mrna_records = []
    gene_id_map = {}  # For gene existence validation

    with open(gff_path) as f:
        for line_num, line in enumerate(f, 1):
            # Skip comments and empty lines
            if line.startswith("#") or not line.strip():
                continue

            # Split fields and validate format
            parts = line.strip().split("\t")
            if len(parts) < 9:
                print(f"Line {line_num} format error, skipping: {line.strip()}")
                continue

            feature_type = parts[2].lower()
            if feature_type not in [
                "gene",
                "mrna",
            ]:  # Process key features according to GFF3 standard
                continue

            try:
                # Parse common fields
                seqid = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6] if parts[6] in ("+", "-", ".") else None

                # Enhanced attribute parsing (compatible with different formats)
                attributes = {}
                for pair in (
                    parts[8].replace("%20", " ").split(";")
                ):  # Handle URL encoding
                    pair = pair.strip()
                    if not pair:
                        continue
                    if "=" in pair:
                        key, value = pair.split("=", 1)
                        attributes[key.strip().lower()] = (
                            value.strip()
                        )  # Unified lowercase processing
                    else:
                        attributes[pair] = None

                # Feature type specific processing
                if feature_type == "gene":
                    # Gene ID parsing (compatible with NCBI/Ensembl different formats)
                    gene_id = (
                        attributes.get("id")
                        or attributes.get("geneid")
                        or f"GENE_{line_num}"
                    )
                    gene_id_map[gene_id] = True  # Register gene ID

                    gene_records.append(
                        {
                            "gene_id": gene_id,
                            "chr": seqid,
                            "gene_start": start,
                            "gene_end": end,
                            "strand": strand,
                            "gene_type": attributes.get("biotype")
                            or attributes.get("type"),
                            "gene_name": attributes.get("name")
                            or attributes.get("gene_name"),
                            # "attributes": str(attributes)
                        }
                    )

                elif feature_type == "mrna":
                    # Associate gene ID (handle different Parent formats)
                    parent_id = attributes.get("parent") or attributes.get("gene_id")
                    if not parent_id:
                        print(
                            f"Line {line_num} mRNA missing Parent attribute: {line.strip()}"
                        )
                        continue

                    # Standardize parent ID (handle ID:XXXX format)
                    parent_id = parent_id.split(":")[-1].replace(
                        "gene_", ""
                    )  # Compatible with NCBI's gene:XXX format

                    # Validate parent gene existence
                    if parent_id not in gene_id_map:
                        print(
                            f"Line {line_num} mRNA's parent gene doesn't exist: {parent_id}"
                        )
                        continue

                    mrna_records.append(
                        {
                            "gene_id": parent_id,
                            "mrna_id": attributes.get("id") or f"mRNA_{line_num}",
                            "mrna_start": start,
                            "mrna_end": end,
                            "exon_count": len([k for k in attributes if "exon" in k]),
                            "transcript_type": attributes.get("transcript_type")
                            or attributes.get("biotype"),
                            "product": attributes.get("product")
                            or attributes.get("description"),
                        }
                    )

            except Exception as e:
                print(f"Failed to parse line {line_num}: {line.strip()}")
                print(f"Error details: {str(e)}")
                continue

    # Create DataFrames and establish relationships
    gene_df = pd.DataFrame(gene_records)
    gene_df["index_temp"] = gene_df.groupby("chr").cumcount() + 1

    mrna_df = pd.DataFrame(mrna_records)

    # Left join merge (keep genes without mRNAs)
    merged_df = pd.merge(gene_df, mrna_df, on="gene_id", how="left")

    # Post-processing
    merged_df.sort_values(["chr", "gene_start", "mrna_start"], inplace=True)
    merged_df.reset_index(drop=True, inplace=True)
    # merged_df.sort_values(by=['chr', 'gene_start'], ascending=False, inplace=True)
    # merged_df['index'] = merged_df.groupby(['chr','gene_name']).ngroup() + 1

    return merged_df


def find_gene_pairs(
    gene_info: pd.DataFrame, gene_id_input: list, id_type: str, distnace: int
) -> pd.DataFrame:
    """
    Filter neighboring gene pairs based on chromosome index range of gene IDs

    Parameters:
        gene_info : DataFrame containing gene information, must include columns:
                   mrna_id, chr, index_temp, gene_id, gene_start, gene_end
        gene_id_input : List of gene IDs or file path

    Returns:
        DataFrame containing all qualified gene pairs
    """
    # Check if required columns exist
    required_cols = {
        "mrna_id",
        "chr",
        "strand",
        "index_temp",
        "gene_id",
        "gene_start",
        "gene_end",
    }
    if not required_cols.issubset(gene_info.columns):
        missing = required_cols - set(gene_info.columns)
        raise KeyError(f"Missing required columns: {missing}")

    # Initialize result container
    all_pairs = []

    # Process input type (support file path or ID list)
    if isinstance(gene_id_input, str):
        with open(gene_id_input) as f:
            gene_ids = [line.strip().split("\t")[0] for line in f]
    else:
        gene_ids = gene_id_input

    # Main processing logic
    for id_temp in gene_ids:
        # Find target gene information
        if id_type == "mrna":
            target_gene = gene_info[gene_info["mrna_id"] == id_temp]
        else:
            target_gene = gene_info[gene_info["gene_id"] == id_temp]

        if target_gene.empty:
            print(f"Warning: Gene {id_temp} not found")
            continue

        # Extract coordinates
        try:
            current_data = target_gene.iloc[0]
            chrom = current_data["chr"]
            idx = int(current_data["index_temp"])
        except (KeyError, ValueError) as e:
            print(f"Data format error: {e}")
            continue

        # Calculate index range
        search_range = (
            (gene_info["chr"] == chrom)
            & (
                gene_info["index_temp"]
                .astype(int)
                .between(idx - distnace, idx + distnace)
            )  # Include boundaries
            & (gene_info["index_temp"] != idx)
        )

        # Filter neighboring genes
        nearby_genes = gene_info[search_range]

        # Generate pairing information
        for _, neighbor in nearby_genes.iterrows():
            pair = {
                "chr": chrom,
                "gene_id_1": current_data["gene_id"],
                "gene_start_1": current_data["gene_start"],
                "gene_end_1": current_data["gene_end"],
                "strand_1": current_data["strand"],
                "gene_id_2": neighbor["gene_id"],
                "gene_start_2": neighbor["gene_start"],
                "gene_end_2": neighbor["gene_end"],
                "strand_2": neighbor["strand"],
            }
            all_pairs.append(pair)

    # Sort and merge col1 and col2 values into new column
    all_pairs = pd.DataFrame(all_pairs) if all_pairs else pd.DataFrame()
    all_pairs["sorted_cols"] = all_pairs.apply(
        lambda row: "_".join(sorted([row["gene_id_1"], row["gene_id_2"]])), axis=1
    )
    df_clean = all_pairs.drop_duplicates(subset=["chr", "sorted_cols"], keep="first")
    df_clean = df_clean.drop(columns=["sorted_cols"])
    df_clean = df_clean.sort_values(by=["chr", "gene_start_1"], ascending=[True, True])

    # 按自然顺序排序
    df_clean = df_clean.sort_values(
        "chr",
        key=lambda x: x.str.replace("chr", "").map(
            lambda s: (s.isdigit(), int(s) if s.isdigit() else s)
        ),
    )
    # 或使用natsort的索引排序
    df_clean = df_clean.iloc[index_natsorted(df_clean["chr"])]

    return df_clean


def main():
    parser = argparse.ArgumentParser(description="Find gene pairs near target genes")
    parser.add_argument("--gff", "-g", required=True, help="Path to GFF file")
    parser.add_argument(
        "--id", "-i", required=True, help="File containing target gene IDs"
    )
    parser.add_argument("--type", "-t", choices=["gene", "mrna"], default="gene")
    parser.add_argument(
        "--distance",
        "-d",
        type=int,
        default=3,
        help="The number of other genes between pairs of genes. The default value is 3.",
    )
    parser.add_argument("--output", "-o", required=True, help="Output filename")
    args = parser.parse_args()

    # 整理基因信息
    # gene_index, gene_info = build_gene_index(args.gff)
    gene_info = parse_gff_features(args.gff)

    # # 处理目标基因
    # with open(args.id) as f:
    #     target_ids = set(line.strip() for line in f)

    # # 转换mRNA到gene
    # if args.type == "mrna":
    #     target_genes = set()
    #     mrna_to_gene = {mrna: gene for gene in gene_info for mrna in gene_info[gene]['transcripts']}
    #     for mrna in target_ids:
    #         if mrna in mrna_to_gene:
    #             target_genes.add(mrna_to_gene[mrna])
    # else:
    #     target_genes = target_ids

    gene_df = pd.read_csv(args.id, header=None, names=["temp_id"])

    # 查找基因对
    # pairs = find_gene_pairs(target_genes, gene_index, args.distance)

    # 写入结果
    # with open(args.output, 'w') as f:
    #     f.write("Chr\tGene.1\tGene.2\tChr\tStart.1\tEnd.1\tStrand.1\tStart.2\tEnd.2\tStrand.2\n")
    #     # 按染色体自然顺序和起始位置排序
    #     def sort_key(pair):
    #         chrom = gene_info[pair[0]]['chrom']
    #         return (natural_chrom_key(chrom), gene_info[pair[0]]['start'])

    #     sorted_pairs = sorted(pairs, key=sort_key)
    #     for gene1, gene2 in sorted_pairs:
    #         chrom1 = gene_info[gene1]['chrom']
    #         s1 = gene_info[gene1]['start']
    #         e1 = gene_info[gene1]['end']
    #         st1 = gene_info[gene1]['strand']
    #         s2 = gene_info[gene2]['start']
    #         e2 = gene_info[gene2]['end']
    #         st2 = gene_info[gene2]['strand']
    #         f.write(f"{chrom1}\t{gene1}\t{gene2}\t{chrom1}\t{s1}\t{e1}\t{st1}\t{s2}\t{e2}\t{st2}\n")

    # 加载基因信息
    result_final = find_gene_pairs(
        gene_info=gene_info,
        gene_id_input=gene_df.temp_id,
        id_type=args.type,
        distnace=args.distance,
    )

    # # 保存结果
    result_final.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
