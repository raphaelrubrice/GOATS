import numpy as np
import gzip
import os
from copy import deepcopy
from pronto import Ontology

def add_item(liste, item):
    if item not in liste:
        liste.append(item)
        return True
    return False

def parse_gaf(file_path, geneset, ontology, collections):
    gene_dico = {'ids':{coll:[] for coll in collections},
                    'str_ids':{coll:[] for coll in collections},
                    'terms':{coll:[] for coll in collections},
                    'predicate':{coll:[] for coll in collections},
                    'all_ids':[],
                    'all_collections':[]}
    info_dico = {gene: deepcopy(gene_dico) for gene in geneset}
    found_flags = {gene:False for gene in geneset}

    translate_collection = {'P':'BP', 'F':'MF', 'C':'CC'}
    description_memory = {}
    with gzip.open(file_path, 'rt') as f:  # Use 'open' instead if file is not gzipped
        for line in f:
            if line.startswith("!"):  # Skip header lines
                continue
            cols = line.strip().split('\t')
            all_names = cols[10].split('|') if not (cols[10] == "") else None
            if all_names == None:
                all_names = set(cols[2]) if not (cols[2] == "") else None
            if all_names != None:
                predicate = cols[3]
                go_term = cols[4]  
                go_number = int(go_term[3:])
                if go_term not in description_memory.keys():
                    go_description = ontology.get(go_term).name
                    description_memory[go_term] = go_description
                else:
                    go_description = description_memory[go_term]
                collection = translate_collection[cols[8]]
                for gene in geneset:
                    gene_dico = info_dico[gene]
                    if gene in all_names:
                        flag_int_id = add_item(gene_dico['ids'][collection], go_number)
                        flag_str_id = add_item(gene_dico['str_ids'][collection], go_term)

                        if flag_int_id and flag_str_id:
                            gene_dico['terms'][collection].append(go_description)
                            gene_dico['predicate'][collection].append(predicate)
                        found_flags[gene] = True
    for gene in geneset:
        if not found_flags[gene]:
            print(f'\nGene {gene} not found. Removing {gene}..')
            del info_dico[gene]
        else:
            for coll in collections:
                info_dico[gene]['all_ids'] = info_dico[gene]['all_ids'] + info_dico[gene]['ids'][coll]
                info_dico[gene]['all_collections'] = info_dico[gene]['all_collections'] + [coll]*len(info_dico[gene]['ids'][coll])
    return info_dico, list(info_dico.keys()), description_memory

def show_relations(infos, gene, format='str', query='predicate', show=True):
    out = {gene:[]} if format=='str' else {}
    if show:
        print(f'\n{gene} :')
    for section in ['MF', 'BP', 'CC']:
        if section in infos[gene]['ids'].keys():
            N = len(infos[gene]['ids'][section])
            for i in range(N):
                go_term = infos[gene]['ids'][section][i]
                pad = ''
                for z in range(7-len(str(go_term))):
                    pad += '0'
                go_relation = infos[gene][query][section][i]
            
                if format=='str':
                    added = {f"GO:{pad+str(go_term)}":go_relation}
                    out[gene].append(added)
                else:
                    out[go_term] = go_relation
                if show:
                    if format=='str':
                        mem = f"'GO:{pad+str(go_term)}'"
                        print(f'\n{mem} : {go_relation}')
                    else:
                        print(f'\n{go_term} : {go_relation}')
    return out

def find_top_terms(infos, geneset, threshold):
    N_geneset = len(list(geneset))
    assert N_geneset > 0, "Empty geneset."

    min_hit = int(np.floor(threshold * N_geneset))
    print(len([go_id for gene in infos.keys() for go_id in infos[gene]['all_ids']]))
    all_terms_geneset = set([go_id for gene in infos.keys() for go_id in infos[gene]['all_ids']])
    print(len(all_terms_geneset))
    kept = {}
    for GO in all_terms_geneset:
        count = 0
        for gene in infos.keys():
            if GO in infos[gene]['all_ids']:
                count += 1
        if count >= min_hit:
            kept[GO] = count / N_geneset
    return kept

def show_description(infos, go_id_list, show=True):
    out = {}
    for go_id in go_id_list:
        count = 0
        description_memory = []
        tuple_list = []
        for gene in infos.keys():
            if go_id in infos[gene]['all_ids']:
                count += 1
                go_id_index = infos[gene]['all_ids'].index(go_id)
                collection = infos[gene]['all_collections'][go_id_index]
                collection_go_index = infos[gene]['ids'][collection].index(go_id)
                gene_predicate = infos[gene]['predicate'][collection][collection_go_index]
                if count == 1:
                    go_description = infos[gene]['terms'][collection][collection_go_index]
                    description_memory.append(go_description)
                    if show:
                        print(f'\n{go_id} : {go_description}')
                if show:
                    print(f'\n{gene} -> {gene_predicate} -> {description_memory[0]}')
                tuple_list.append((gene, gene_predicate, description_memory[0]))
        out[go_id] = tuple_list
    return out

def load_config(config_file):
    """
    Loads configuration settings from a file into a dictionary.

    Parameters:
    config_file: str
        Path to the configuration file.

    Returns:
    dict
        A dictionary with configuration key-value pairs.
    """
    config = {}
    try:
        with open(config_file, 'r') as file:
            for line in file:
                # Ignore comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                # Split the line into key and value
                key, value = line.strip().split('=', 1)
                config[key.strip()] = value.strip()
    except FileNotFoundError:
        print(f"Error: Configuration file {config_file} not found.")
        exit(1)
    except Exception as e:
        print(f"Error while reading the configuration file: {e}")
        exit(1)
    
    return config

def create_info_dico_string(info_dico, genes, collection):
    """
    Create a string representation of the `info_dico` dictionary for saving to a text file.

    Args:
        info_dico (dict): Dictionary containing parsed GAF data for each gene.

    Returns:
        str: A formatted string representation of `info_dico`.
    """
    output_lines = []
    translate = {'MF':"  Molecular Function (MF):\n", "BP":"  Biological Process (BP):\n", "CC":"  Cellular Component (CC):\n"}
    for gene, data in info_dico.items():
        if gene in genes:
            for coll in collection:
                output_lines.append(f"Gene: {gene}\n")
                output_lines.append(translate[coll])
                for idx, go_id in enumerate(data['ids'][coll]):
                    term = data['terms'][coll][idx]
                    predicate = data['predicate'][coll][idx]
                    output_lines.append(f"    GO ID: {go_id} | Term: {term} | Predicate: {predicate}\n")

            # Add summary information
            output_lines.append("  All GO IDs:\n")
            output_lines.append(f"    {', '.join(map(str, data['all_ids']))}\n")
            output_lines.append("  All Collections:\n")
            output_lines.append(f"    {', '.join(data['all_collections'])}\n")
            output_lines.append("\n")

    # Join all lines into a single string
    return "".join(output_lines)

def create_top_terms_string(top_terms):
    """
    Create a string representation of the `top_terms` dictionary for saving to a text file.

    Args:
        top_terms (dict): Dictionary mapping GO terms to their relative frequency in a gene set.

    Returns:
        str: A formatted string representation of `top_terms`.
    """
    output_lines = []
    output_lines.append("Top GO Terms by Frequency in Gene Set:\n")
    output_lines.append("=======================================\n")
    
    # Format each term and its frequency
    for go_term, frequency in sorted(top_terms.items(), key=lambda x: -x[1]):  # Sort by frequency, descending
        output_lines.append(f"GO Term: GO:{str(go_term).zfill(7)} | Frequency: {frequency:.4f}\n")
    
    # Join all lines into a single string
    return "".join(output_lines)

def create_description_string(descriptions):
    """
    Create a string representation of the `descriptions` dictionary for saving to a text file.

    Args:
        descriptions (dict): Dictionary mapping GO terms to lists of tuples with gene, predicate, and description.

    Returns:
        str: A formatted string representation of `descriptions`.
    """
    output_lines = []
    output_lines.append("GO Term Descriptions and Gene Relationships:\n")
    output_lines.append("===========================================\n")

    for go_id, relationships in descriptions.items():
        output_lines.append(f"\nGO Term: GO:{str(go_id).zfill(7)}\n")
        output_lines.append("-------------------------------------------\n")

        if not relationships:
            output_lines.append("No associated genes found.\n")
        else:
            for gene, predicate, description in relationships:
                output_lines.append(f"Gene: {gene} | Predicate: {predicate} | Description: {description}\n")
    
    # Join all lines into a single string
    return "".join(output_lines)

def run_goats(config):
    path_geneset = config["GENESET"]
    geneset = []
    with open(path_geneset, 'r') as f:
        for line in f:
            gene_symbol = line.strip()
            if gene_symbol and 'str' in str(type(gene_symbol)):  # Only append if it's not an empty line
                geneset.append(gene_symbol)

    obo_file = config["OBO_FILE"]
    gaf_file = config["GAF_FILE"]
    collections = config["COLLECTIONS"].split(',')
    collections = [coll.replace("'", '') for coll in collections]
    threshold = float(config["THRESHOLD"])
    show_desc = True if int(config["SHOW_DESC"]) == 1 else False
    save = True if int(config["SAVE"]) == 1 else False
    folder = config["FOLDERNAME"]
    timestamp = config["TIMESTAMP"].replace("\"", "")
    if folder == "None":
        folder = f"{timestamp}_outputs"
    
    if save:
        try:
            os.mkdir(folder)
        except:
            folder = folder + "_" + timestamp
            os.mkdir(folder)
    os.rename(f"log_{timestamp}.txt", f"log_{folder}.txt")
    
    print("\nLoading Ontology..")
    ontology = Ontology(obo_file)

    print("\nRetrieving informations..")
    infos, genes_in_info, all_descriptions = parse_gaf(gaf_file, geneset, ontology, collections)
    if save:
        filename = folder + "/" + f"infos_goats_query_{timestamp}.txt"
        txt = create_info_dico_string(infos, genes_in_info, collections)
        with open(filename, 'w') as f:
            f.write(txt)
    print("\nFinding frequent GO terms..")
    hits = find_top_terms(infos, genes_in_info, threshold)
    if save:
        filename = folder + "/" + f"hits_goats_query_{timestamp}.txt"
        txt = create_top_terms_string(hits)
        with open(filename, 'w') as f:
            f.write(txt)
    print("\nMaking description object..")
    out = show_description(infos, hits.keys(), show=show_desc)
    if save:
        filename = folder + "/" + f"desc_goats_query_{timestamp}.txt"
        txt = create_description_string(out)
        with open(filename, 'w') as f:
            f.write(txt)
    print("\nQuery Completed.")
    if save:
        print(f"You can find outputs in ./{folder} folder")
    return infos, hits, out

if __name__ == "__main__":
    #Load the configuration file path from command-line arguments
    import sys
    if len(sys.argv) != 2:
        print("Usage: python goats.py <config_file>")
        exit(1)

    config_file = sys.argv[1]
    config = load_config(config_file)
    infos, hits, out = run_goats(config)

