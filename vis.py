from pathlib import Path
import os
from datetime import datetime
import re
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag  
import goatools.base
import subprocess
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from PIL import Image


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

def get_color_gradient(counts, min_count=None, max_count=None):
    """
    Generate a gradient color map based on counts.
    :param counts: A list of count values.
    :param min_count: Minimum value for the gradient (auto-calculated if None).
    :param max_count: Maximum value for the gradient (auto-calculated if None).
    :return: A tuple containing:
        - A dictionary of GO terms mapped to gradient colors.
        - The colormap used for the gradient.
        - The normalization object for the color scale.
    """
    min_count = min_count if min_count is not None else 0
    max_count = max_count if max_count is not None else max(counts.values())

    cmap = plt.get_cmap("Reds")  # Use a red gradient colormap
    norm = mcolors.Normalize(vmin=min_count, vmax=max_count)

    # Map counts to colors
    color_map = {}
    for go_term, count in counts.items():
        color = mcolors.to_hex(cmap(norm(count)))
        color_map[go_term] = color
    return color_map, cmap, norm





def extract_go_terms_with_counts(filepath):
    """
    Extract GO terms and their counts from the input file.
    Example file content: 'GO:0005515 | Frequency: 0.8889'
    :param filepath: Path to the input file.
    :return: A dictionary with GO terms as keys and counts as values.
    """
    go_counts = {}
    pattern = r'(GO:\d+)\s*\|\s*Frequency:\s*(\d*\.\d+)'

    with filepath.open('r', encoding='utf-8') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                go_term = match.group(1)
                count = float(match.group(2))
                go_counts[go_term] = count

    if not go_counts:
        raise ValueError("No GO terms with counts found in the file.")
    return go_counts





def plot_go_lineage_with_colors(go_terms, go_counts):
    """
    Generate a GO lineage plot with parent relationships and gradient colors based on counts.
    GO terms not in the list are displayed in grey.
    :param go_terms: Set of GO terms to visualize.
    :param go_counts: Dictionary with GO terms as keys and count values.
    """
    
    go_dag = GODag("go-basic.obo")  
    color_map, cmap, norm = get_color_gradient(go_counts)
    go_subdag = GoSubDag(go_terms, go_dag, prt=None)

    dot_filename = "go_lineage_colored.dot"
    image_filename = "go_lineage_colored.png"
    with open(dot_filename, 'w') as dot_file:
        dot_file.write("digraph GO {\n")
        for goid in go_subdag.go2obj:
            if goid in go_counts:
                color = color_map[goid]  # Color based on counts
                label = f"{goid}\n{go_counts[goid]:.2f}"
            else:
                color = "#D3D3D3"  # Grey for GO terms not in the list
                label = goid
            dot_file.write(f'    "{goid}" [label="{label}", style=filled, fillcolor="{color}"];')
        for goid in go_subdag.rcntobj.go2parents:
            for parent in go_subdag.rcntobj.go2parents[goid]:
                dot_file.write(f'    "{parent}" -> "{goid}";\n')
        dot_file.write("}\n")

    
    try:
        subprocess.run(['dot', '-Tpng', dot_filename, '-o', image_filename], check=True)
        print(f"GO lineage plot with colors saved as '{image_filename}'")
    except subprocess.CalledProcessError as e:
        print(f"Error creating GO lineage plot: {e}")

    # Combine the lineage plot and legend
    combined_filename = "go_lineage_with_legend.png"
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax, orientation='horizontal')
    cbar.set_label("GO Term Frequency")

    # Save the legend as a temporary image
    legend_filename = "legend_temp.png"
    plt.savefig(legend_filename)
    plt.close()

    # Combine the lineage image and legend
    lineage_img = Image.open(image_filename)
    legend_img = Image.open(legend_filename)
    total_height = lineage_img.height + legend_img.height
    combined_img = Image.new("RGBA", (lineage_img.width, total_height), (255, 255, 255, 255))
    combined_img.paste(lineage_img, (0, 0))
    combined_img.paste(legend_img, (0, lineage_img.height))
    combined_img.save(combined_filename)

    print(f"Combined lineage plot with legend saved as '{combined_filename}'")



def main():
    try:
        
        most_recent_hits_file = get_most_recent_hits_file()  
        go_counts = extract_go_terms_with_counts(most_recent_hits_file)
        go_terms = set(go_counts.keys())
        print(f"Extracted GO Terms: {go_terms}")
        plot_go_lineage_with_colors(go_terms, go_counts)

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()