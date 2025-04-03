import argparse
import sys
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Process BLAST results and filter by E-value')
    parser.add_argument('-i', '--input', required=True, help='Input BLAST result file')
    parser.add_argument('-e', '--evalue', type=float, default=1e-5, help='E-value threshold (default: 1e-5)')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('-n', '--number', type=int, default=1, help='Number of top hits to retain per query (default: 1)')
    
    args = parser.parse_args()
    
    # Store best matches for each gene
    best_hits = defaultdict(list)
    
    try:
        with open(args.input, 'r') as f_in, open(args.output, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                
                evalue = float(fields[10])
                if evalue < args.evalue:
                    continue
                
                gene_id = fields[0]
                pident = float(fields[2])
                bitscore = float(fields[11])
                best_hits[gene_id].append((pident, bitscore, line))
            
            # Filter highest scoring records for each gene
            total_hits = 0
            header = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
            f_out.write(header)
            for gene, records in best_hits.items():
                if records:
                    sorted_records = sorted(records, key=lambda x: (-x[0], -x[1]))
                    # Validate field index validity
                    if len(fields) < 12:
                        print(f"Warning: Line {line.strip()} lacks required fields, skipped")
                        continue
                    for record in sorted_records[:args.number]:
                        f_out.write(record[2])
                        total_hits += 1
            
            print(f"Successfully processed {len(best_hits)} query sequences, retained {total_hits} records.")
        print(f"Results saved to: {args.output}")
    
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()