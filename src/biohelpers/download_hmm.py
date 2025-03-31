import argparse
import requests
import os
import sys
import gzip
import shutil
import tempfile

def main():
    parser = argparse.ArgumentParser(description='Download HMM profile from InterPro')
    parser.add_argument('-id', '--hmm_id', required=True, help='Pfam HMM ID (e.g. PF00010)')
    parser.add_argument('-o', '--output', required=True, help='Output directory path')
    
    args = parser.parse_args()
    
    # Construct API URL
    url = f'https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{args.hmm_id}?annotation=hmm'
    
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        # Create output directory if needed
        os.makedirs(args.output, exist_ok=True)
        output_file = os.path.join(args.output, f'{args.hmm_id}.hmm')
        
        # 使用临时文件保存压缩内容
        # Use temporary file to save compressed content
        with tempfile.NamedTemporaryFile(delete=False, suffix='.gz') as tmp_file:
            for chunk in response.iter_content(chunk_size=8192):
                tmp_file.write(chunk)
            tmp_path = tmp_file.name
        
        # 解压gzip文件到目标路径
        # Extract gzip file to target path
        with gzip.open(tmp_path, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # 清理临时文件
        # Cleanup temporary files
        os.unlink(tmp_path)
        
        print(f'Successfully downloaded HMM profile to {output_file}')
        
    except requests.exceptions.HTTPError as e:
        print(f'Error downloading HMM: {e}')
        if response.status_code == 404:
            print(f'HMM ID {args.hmm_id} not found')
        sys.exit(1)
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)

if __name__ == '__main__':
    main()