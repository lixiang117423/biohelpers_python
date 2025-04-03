import argparse
import os
import subprocess
import sys
from pathlib import Path

def validate_fastq_files(folder):
    """Validate fastq.gz files in the directory"""
    fastq_files = list({f for f in Path(folder).glob('*.*') 
        if f.suffix in ('.fq', '.fastq') or 
        ''.join(f.suffixes[-2:]) in ('.fq.gz', '.fastq.gz')})
    if not fastq_files:
        raise FileNotFoundError(f"No fastq.gz files found in {folder}")
    return fastq_files

def pair_end_files(fastq_files, index):
    """Pair-end sequencing files matching"""
    pairs = {}
    for f in sorted(fastq_files):
        import re
        base_pattern = re.compile(r'(.+?)_[12](?:\..+)?$')
        base_match = base_pattern.match(f.name.split('.', 1)[0])
        if base_match:
            base = base_match.group(1)
        else:
            continue
        
        if '_1' in f.name:
            pairs.setdefault(base, {})['R1'] = f
        elif '_2' in f.name:
            pairs.setdefault(base, {})['R2'] = f
    
    unpaired = [str(f) for p in pairs.values() for f in p.values() if len(p) != 2]
    if unpaired:
        raise ValueError(f"Unpaired files detected: {', '.join(unpaired)}")
    
    return [(
        str(pair['R1']), 
        str(pair['R2']), 
        f"{base}.sorted.bam"
    ) for base, pair in pairs.items()]

def run_hisat2(index, threads, folder, output, method):
    """Execute HISAT2 alignment workflow"""
    os.makedirs(output, exist_ok=True)
    
    try:
        fastq_files = validate_fastq_files(folder)
        pairs = pair_end_files(fastq_files, index)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    for r1, r2, bam_name in pairs:
        sam_path = Path(output) / f"{Path(bam_name).stem}.sam"
        bam_path = Path(output) / bam_name
        
        # 合并hisat2比对和samtools排序
        pipe_cmd = f"hisat2 -x {index} -1 {r1} -2 {r2} -p {threads} | samtools sort -@ {threads} -O BAM -o {bam_path} -"
        
        if method == 'run':
            try:
                subprocess.run(pipe_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Command execution failed: {e}")
                sys.exit(1)
        # Merge hisat2 alignment and samtools sorting
        # Bugfix: Original code `elif method == 'save':` had syntax issues, need to add logic block
        pass  # Placeholder for future logic implementation
            script_path = Path("run_hisat2.sh")
            with open(script_path, 'w') as f:
                pass
            os.chmod(script_path, 0o755)
            for r1, r2, bam_name in pairs:
                sam_path = Path(output) / f"{Path(bam_name).stem}.sam"
                bam_path = Path(output) / bam_name
                
                pipe_cmd = f"hisat2 -x {index} -1 {r1} -2 {r2} -p {threads} | samtools sort -@ {threads} -O BAM -o {bam_path} -"
                
                with open(script_path, 'a') as f:
                    f.write(f"{pipe_cmd}\n")
        print(f"\nPlease run the following command to execute the alignment:\nbash {script_path}\n")

def main():
    parser = argparse.ArgumentParser(description='HISAT2 RNA-seq alignment pipeline')
    parser.add_argument('-x', '--index', required=True, help='Reference genome index path')
    parser.add_argument('-t', '--threads', type=int, default=os.cpu_count(), 
                       help='Number of threads (default: all cores)')
    parser.add_argument('-f', '--folder', required=True, 
                       help='Input directory containing FASTQ files')
    parser.add_argument('-m', '--method', choices=['run', 'save'], default='run',
                       help='Execution method: run immediately or save to script')
    parser.add_argument('-o', '--output', required=True, 
                       help='Output directory for BAM files')
    
    args = parser.parse_args()
    index_files = list(Path(args.index).parent.glob(f"{Path(args.index).name}*.ht2"))
    if not index_files:
        print(f"Error: No HISAT2 index files found with base name '{args.index}'. Please check the index path.")
        print(f"Expected files like: {args.index}.1.ht2, {args.index}.2.ht2, etc.")
        sys.exit(1)
    
    run_hisat2(args.index, args.threads, args.folder, args.output, args.method)

if __name__ == '__main__':
    main()