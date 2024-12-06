from pathlib import Path
import os
from datetime import datetime
import re
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag  
import goatools.base
import subprocess
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



def extract_go_terms(filepath):
    go_terms = []
    pattern = r'GO:\d+'

    with filepath.open('r', encoding='utf-8') as file:
        content = file.read()
        go_terms = re.findall(pattern, content)

    if not go_terms:
        raise ValueError("No GO terms found in the provided file.")

    return set(go_terms)  

def plot_go_lineage_parent(go_terms):
    
    go_dag = GODag("go-basic.obo")  

    go_subdag = GoSubDag(go_terms, go_dag, prt=None)

    # Write the GO subdag to a DOT file for visualization
    dot_filename = "go_lineage_parent.dot"
    with open(dot_filename, 'w') as dot_file:
        
        
        dot_file.write("digraph GO {\n")
        
        
        for goid in go_subdag.rcntobj.go2parents:
            dot_file.write(f'    "{goid}" [label="{goid}"];\n')
            for parent in go_subdag.rcntobj.go2parents[goid]:
                dot_file.write(f'    "{parent}" -> "{goid}";\n')
        
        
        dot_file.write("}\n")
    
    # Convert the DOT file to a PNG image using graphviz (requires graphviz installed)
    os.system(f"dot -Tpng {dot_filename} -o go_lineage_parent.png")
    print(f"GO lineage plot saved as 'go_lineage_parent.png'")


def main():
    try:
        most_recent_hits_file = get_most_recent_hits_file()
        print(f"Reading file: {most_recent_hits_file}")

        
        go_terms = extract_go_terms(most_recent_hits_file)
        print(f"Extracted GO Terms: {go_terms}")

        
        plot_go_lineage_parent(go_terms)


    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()