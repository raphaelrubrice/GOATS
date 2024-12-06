from pathlib import Path
import os
from datetime import datetime

base_dir = Path(__file__).parent

def get_most_recent_hits_file():
    
    """
    Returns the most recent hits_goats_query_*.txt file in the current directory
    
    """
    
    output_dirs = [d for d in base_dir.glob('*_outputs') if d.is_dir()]

    if not output_dirs:
        raise FileNotFoundError("No output directories found.")

    output_dirs.sort(key=lambda d: d.stat().st_mtime, reverse=True)
    most_recent_output_dir = output_dirs[0]

    hits_files = list(most_recent_output_dir.glob('hits_goats_query_*.txt'))

    if not hits_files:
        raise FileNotFoundError("No hits_goats_query_*.txt files found in the most recent output directory.")

    hits_files.sort(key=lambda f: f.stat().st_mtime, reverse=True)
    most_recent_hits_file = hits_files[0]

    return most_recent_hits_file

def main():
    try:
        most_recent_hits_file = get_most_recent_hits_file()
        print(f"Reading file: {most_recent_hits_file}")

        with most_recent_hits_file.open('r', encoding='utf-8') as file:
            content = file.read()
            print(content)

    except FileNotFoundError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()