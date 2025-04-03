from collections import defaultdict
import os
import subprocess
import tempfile
from collections import defaultdict
from Bio import SeqIO

class GFF3Parser:
    def __init__(self, gff3_path):
        self.gff3_path = gff3_path
        self.transcripts = defaultdict(list)

    def parse(self):
        with open(self.gff3_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2]
                attributes = dict([item.split('=') for item in fields[8].split(';') if '=' in item])
                
                if feature_type == 'mRNA':
                    parent_gene = attributes.get('Parent', '').split(':')[-1]
                    self.transcripts[parent_gene].append({
                        'transcript_id': attributes['ID'],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'exons': []
                    })
                elif feature_type == 'exon':
                    parent_transcript = attributes.get('Parent', '').split(':')[-1]
                    for gene in self.transcripts.values():
                        for transcript in gene:
                            if transcript['transcript_id'] == parent_transcript:
                                transcript['exons'].append((int(fields[3]), int(fields[4])))

class TranscriptProcessor:
    @staticmethod
    def calculate_cds_length(transcript):
        return sum(end - start + 1 for start, end in transcript['exons'])

    @classmethod
    def get_longest_transcript(cls, gene_transcripts):
        return max(gene_transcripts, key=lambda x: cls.calculate_cds_length(x))

class CDSCalculator:
    def __init__(self):
        self.transcript_lengths = defaultdict(int)
        self.gene_metadata = defaultdict(dict)
        self.current_chromosome = None

    def calculate_from_gff(self, gff_path):
        gene_transcripts = defaultdict(list)
        with open(gff_path) as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    if line.startswith('##sequence-region'):
                        self.current_chromosome = line.split()[1]
                    continue
                fields = line.strip().split('\t')
                chrom = fields[0]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                
                if feature_type == 'gene':
                    gene_id = dict(item.split('=') for item in fields[8].split(';') if '=' in item)['ID']
                    self.gene_metadata[gene_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }
                
                if feature_type == 'mRNA':
                    attrs = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                    transcript_id = attrs['ID']
                    parent_gene = attrs.get('Parent', '').split(':')[-1]
                    gene_transcripts[parent_gene].append({
                        'id': transcript_id,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'chrom': chrom
                    })
                elif feature_type == 'CDS':
                    attrs = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                    transcript_id = attrs['Parent'].split(':')[-1]
                    length = end - start + 1
                    self.transcript_lengths[transcript_id] += length
        
        longest_transcripts = {}
        for gene_id, transcripts in gene_transcripts.items():
            longest_id = max(transcripts, key=lambda x: self.transcript_lengths.get(x['id'], 0))
            longest_transcripts[gene_id] = longest_id
        self.gene_transcripts = gene_transcripts
        return longest_transcripts

import argparse

def main():
    parser = argparse.ArgumentParser(description='Extract longest mRNA transcripts')
    parser.add_argument('-g', '--genome', required=True, help='Input genome FASTA file')
    parser.add_argument('-f', '--gff3', required=True, help='Input GFF3 annotation file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    args = parser.parse_args()
    
    parse_longest_mrna(args.genome, args.gff3, args.output)


def parse_longest_mrna(genome_path, gff3_path, output_path):
    # Calculate longest transcripts
    cds_calculator = CDSCalculator()
    longest_transcripts = cds_calculator.calculate_from_gff(gff3_path)
    
    # Generate gene info output path
    genome_dir = os.path.dirname(genome_path)
    genome_basename = os.path.basename(genome_path).split('.')[0]
    gene_info_path = os.path.join(genome_dir, f'{genome_basename}.gene.info.txt')
    
    # Write gene info file
    with open(gene_info_path, 'w') as info_file:
        info_file.write("mRNA_ID\tgene_ID\tmRNA_start\tmRNA_end\tgene_start\tgene_end\tstrand\tchr\n")
        for gene_id, transcript_info in longest_transcripts.items():
            gene_data = cds_calculator.gene_metadata.get(gene_id, {})
            info_file.write(f"{transcript_info['id']}\t{gene_id}\t"
                            f"{transcript_info['start']}\t{transcript_info['end']}\t"
                            f"{gene_data.get('start', '')}\t{gene_data.get('end', '')}\t"
                            f"{transcript_info['strand']}\t{transcript_info['chrom']}\n")
    
    # 输出统计信息
    # Output statistics
    print(f"Total genes processed: {gene_count}")
    print(f"Genes with multiple transcripts: {multi_isoform_genes}")
    print(f"Average transcript length: {avg_length:.2f}")
    
    # Temporary file path for transcript IDs
    # Preserve temporary files for debugging
    # Generate temporary ID file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, encoding='utf-8') as id_file:
        for gene_id, transcript_info in longest_transcripts.items():
            id_file.write(f"{transcript_info['id']}\n")
        temp_id_path = id_file.name
    
    # Debug log: output temp file path
    # 转录本ID的临时文件路径
    # Temporary file path for transcript IDs
    # print(f'Generated temp ID file path: {temp_id_path}')
    
    # Extract sequences using seqkit
    # Generate temporary protein sequence file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fa') as protein_file:
        # print(f'GFFREAD output path: {protein_file.name}')
        genome_norm = os.path.normpath(genome_path)
        gff3_norm = os.path.normpath(gff3_path)
        protein_norm = os.path.normpath(protein_file.name)
        gffread_cmd = f'gffread "{gff3_norm}" -g "{genome_norm}" -y "{protein_norm}"'
        subprocess.run(gffread_cmd, shell=True, check=True)
    
    # print(f'Temporary ID file: {temp_id_path}')
    # print(f'Protein file: {protein_file.name}')
    temp_id_path = os.path.normpath(temp_id_path)
    protein_path = os.path.normpath(protein_file.name)
    output_path = os.path.normpath(output_path)
    cmd = f'seqkit grep -f "{temp_id_path}" "{protein_path}" -o "{output_path}"'
    # print(f'Running seqkit to parse protein sequences of the longest  transcript: {cmd}')  # Debug log: output full command
    print("################################################################")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, encoding='utf-8')

    # subprocess.run(cmd, shell=True, capture_output=True, text=True, encoding='utf-8')
    
    # Debug log: output command execution results
    if not os.path.exists(protein_file.name):
        raise FileNotFoundError(f'Protein file {protein_file.name} not generated')
    if result.returncode != 0:
        print(f'Error exit code: {result.returncode}')
        print(f'Error output:\n{result.stderr}')
        raise RuntimeError(f'seqkit execution failed:\n{result.stderr}')
    else:
        print(f'Successfully extracted {len(longest_transcripts)} longest transcripts')
        print(f'Longest transcript protein sequences saved to: {output_path}')
        print(f'Gene and transcript information saved to: {gene_info_path}')
    
    # Clean up temporary files
    # 保留临时文件用于调试
    # Retain temporary files for debugging
    os.unlink(temp_id_path)
    os.unlink(protein_file.name)
    # print(f'Retained temporary files:\n{temp_id_path}\n{protein_file.name}')
    # print('Temporary files cleaned up')

if __name__ == '__main__':
    main()