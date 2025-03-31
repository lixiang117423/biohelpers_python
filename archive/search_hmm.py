import argparse
import requests
import os
import sys
import gzip
import re
import shutil
from pathlib import Path
from datetime import datetime, timedelta

CACHE_DIR = os.path.join(Path.home(), '.biohelpers_cache')
DAT_FILE_URL = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
CACHE_EXPIRY_DAYS = 7

class HMMDatabase:
    def __init__(self):
        os.makedirs(CACHE_DIR, exist_ok=True)
        self.cache_file = os.path.join(CACHE_DIR, 'Pfam-A.hmm.dat')
        self.mapping = {}
        
    def needs_refresh(self):
        if not os.path.exists(self.cache_file):
            return True
        mod_time = datetime.fromtimestamp(os.path.getmtime(self.cache_file))
        return datetime.now() - mod_time > timedelta(days=CACHE_EXPIRY_DAYS)

    def download_database(self):
        os.makedirs(CACHE_DIR, exist_ok=True)
        print(f'Downloading Pfam database from {DAT_FILE_URL}')
        
        try:
            response = requests.get(DAT_FILE_URL, stream=True)
            response.raise_for_status()

            with gzip.open(response.raw, 'rt') as f_in:
                with open(self.cache_file, 'w', encoding='utf-8') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            print(f'Database cached at {self.cache_file}')
            print(f'Cache directory: {os.path.abspath(CACHE_DIR)}')
        except Exception as e:
            print(f'Failed to download database: {e}')
            sys.exit(1)

    def build_mapping(self):
        if not os.path.exists(self.cache_file):
            print(f'Error: Database file not found at {self.cache_file}')
            print('Please run the download first')
            sys.exit(1)
        
        self.total_entries = 0
        current_id = None
        current_desc = []
        print('Building description-ID mapping...')

        with open(self.cache_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#=GF AC'):
                    current_id = line.split()[-1].split('.')[0]  # Get PFxxxxx without version
                elif line.startswith('#=GF DE'):
                    current_desc.append(' '.join(line.split()[3:]).strip())
                elif line == '//' and current_id and current_desc:
                    full_desc = ' '.join(current_desc)
                    processed_desc = re.sub(r'\W+', '', full_desc).lower()
                    self.mapping[processed_desc] = (full_desc, current_id)
                    self.total_entries = len(self.mapping)
                    current_id = None
                    current_desc = []

    def search(self, keyword):
        processed_keyword = re.sub(r'[^\w-]+', '', keyword).lower().replace(' ', '')
        return {original_desc: hmm_id for processed_desc, (original_desc, hmm_id) in self.mapping.items() if processed_keyword in processed_desc}

def main():
    parser = argparse.ArgumentParser(description='Search Pfam HMM by description')
    parser.add_argument('-d', '--description', required=True, help='Search keyword for HMM description')
    
    args = parser.parse_args()
    
    db = HMMDatabase()
    
    if db.needs_refresh():
        db.download_database()
    
    db.build_mapping()
    results = db.search(args.description)
    
    if not results:
        print(f'No matching HMMs found in {self.total_entries} entries')
        print(f'Check cache file: {os.path.abspath(self.cache_file)}')
        sys.exit(1)
    
    output = [f'{hmm_id}\t{desc}' for desc, hmm_id in results.items()]
    
    print('\n'.join(output))

if __name__ == '__main__':
    main()