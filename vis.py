import networkx as nx
import matplotlib.pyplot as plt


# {'query': 'GAPDH', '_id': '2597', '_score': 17.54872, 'go': {'BP': [{'evidence': 'ISS', 'gocategory': 'BP', 'id': 'GO:0000226', 'qualifier': 'involved_in', 'term': 'microtubule cytoskeleton organization'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0001819', 'pubmed': 22832495, 'qualifier': 'involved_in', 'term': 'positive regulation of cytokine production'},
#                                                                     {'evidence': 'IEA', 'gocategory': 'BP', 'id': 'GO:0006006', 'qualifier': 'involved_in', 'term': 'glucose metabolic process'},
#                                                                     {'evidence': 'IBA', 'gocategory': 'BP', 'id': 'GO:0006096', 'qualifier': 'involved_in', 'term': 'glycolytic process'},
#                                                                     {'evidence': 'IEA', 'gocategory': 'BP', 'id': 'GO:0006096', 'qualifier': 'involved_in', 'term': 'glycolytic process'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0010951', 'pubmed': 22832495, 'qualifier': 'involved_in', 'term': 'negative regulation of endopeptidase activity'},
#                                                                     {'evidence': 'TAS', 'gocategory': 'BP', 'id': 'GO:0016241', 'pubmed': 26626483, 'qualifier': 'involved_in', 'term': 'regulation of macroautophagy'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0017148', 'pubmed': 23071094, 'qualifier': 'involved_in', 'term': 'negative regulation of translation'},
#                                                                     {'evidence': 'IMP', 'gocategory': 'BP', 'id': 'GO:0017148', 'pubmed': 15479637, 'qualifier': 'involved_in', 'term': 'negative regulation of translation'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0031640', 'pubmed': 22832495, 'qualifier': 'involved_in', 'term': 'killing of cells of another organism'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0032481', 'pubmed': 27387501, 'qualifier': 'involved_in', 'term': 'positive regulation of type I interferon production'}, 
#                                                                     {'evidence': 'ISS', 'gocategory': 'BP', 'id': 'GO:0035606', 'qualifier': 'involved_in', 'term': 'peptidyl-cysteine S-trans-nitrosylation'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0043123', 'pubmed': 23332158, 'qualifier': 'involved_in', 'term': 'positive regulation of canonical NF-kappaB signal transduction'}, 
#                                                                     {'evidence': 'ISS', 'gocategory': 'BP', 'id': 'GO:0050821', 'qualifier': 'involved_in', 'term': 'protein stabilization'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0050832', 'pubmed': 22832495, 'qualifier': 'involved_in', 'term': 'defense response to fungus'},
#                                                                     {'evidence': 'ISS', 'gocategory': 'BP', 'id': 'GO:0051402', 'qualifier': 'involved_in', 'term': 'neuron apoptotic process'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0051873', 'pubmed': 22832495, 'qualifier': 'involved_in', 'term': 'killing by host of symbiont cells'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0061844', 'pubmed': 22832495, 'qualifier': 'involved_in', 'term': 'antimicrobial humoral immune response mediated by antimicrobial peptide'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'BP', 'id': 'GO:0071346', 'pubmed': 15479637, 'qualifier': 'involved_in', 'term': 'cellular response to type II interferon'}], 
#                                                              'CC': [{'evidence': 'HDA', 'gocategory': 'CC', 'id': 'GO:0005634', 'pubmed': 21630459, 'qualifier': 'located_in', 'term': 'nucleus'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0005634', 'pubmed': 28404743, 'qualifier': 'located_in', 'term': 'nucleus'}, 
#                                                                     {'evidence': 'ISS', 'gocategory': 'CC', 'id': 'GO:0005634', 'qualifier': 'located_in', 'term': 'nucleus'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0005737', 'pubmed': [11785981, 24507776], 'qualifier': 'located_in', 'term': 'cytoplasm'},
#                                                                     {'evidence': 'ISS', 'gocategory': 'CC', 'id': 'GO:0005737', 'qualifier': 'located_in', 'term': 'cytoplasm'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0005811', 'pubmed': 14741744, 'qualifier': 'located_in', 'term': 'lipid droplet'}, 
#                                                                     {'evidence': 'IBA', 'gocategory': 'CC', 'id': 'GO:0005829', 'qualifier': 'is_active_in', 'term': 'cytosol'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0005829', 'pubmed': [24101517, 28404743], 'qualifier': 'located_in', 'term': 'cytosol'}, 
#                                                                     {'evidence': 'ISS', 'gocategory': 'CC', 'id': 'GO:0005829', 'qualifier': 'located_in', 'term': 'cytosol'},
#                                                                     {'evidence': 'TAS', 'gocategory': 'CC', 'id': 'GO:0005829', 'qualifier': 'located_in', 'term': 'cytosol'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0005886', 'qualifier': 'located_in', 'term': 'plasma membrane'}, 
#                                                                     {'evidence': 'ISS', 'gocategory': 'CC', 'id': 'GO:0015630', 'qualifier': 'located_in', 'term': 'microtubule cytoskeleton'}, 
#                                                                     {'evidence': 'HDA', 'gocategory': 'CC', 'id': 'GO:0016020', 'pubmed': 19946888, 'qualifier': 'located_in', 'term': 'membrane'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0031965', 'qualifier': 'located_in', 'term': 'nuclear membrane'}, 
#                                                                     {'evidence': 'HDA', 'gocategory': 'CC', 'id': 'GO:0031982', 'pubmed': 19190083, 'qualifier': 'located_in', 'term': 'vesicle'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0043231', 'qualifier': 'located_in', 'term': 'intracellular membrane-bounded organelle'}, 
#                                                                     {'evidence': 'IEA', 'gocategory': 'CC', 'id': 'GO:0048471', 'qualifier': 'located_in', 'term': 'perinuclear region of cytoplasm'}, 
#                                                                     {'evidence': 'HDA', 'gocategory': 'CC', 'id': 'GO:0070062', 'pubmed': [11487543, 12519789, 19056867, 19199708, 20458337, 21362503, 23533145], 'qualifier': 'located_in', 'term': 'extracellular exosome'},
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:0097452', 'pubmed': [15479637, 23071094], 'qualifier': 'part_of', 'term': 'GAIT complex'}, 
#                                                                     {'evidence': 'IDA', 'gocategory': 'CC', 'id': 'GO:1990904', 'pubmed': 15479637, 'qualifier': 'part_of', 'term': 'ribonucleoprotein complex'}],
#                                                              'MF': [{'category': 'MF', 'evidence': 'IBA', 'id': 'GO:0004365', 'qualifier': 'enables', 'term': 'glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity'}, 
#                                                                     {'category': 'MF', 'evidence': 'ISS', 'id': 'GO:0004365', 'qualifier': 'enables', 'term': 'glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity'}, 
#                                                                     {'category': 'MF', 'evidence': 'NAS', 'id': 'GO:0004365', 'pubmed': [3170585, 7030790], 'qualifier': 'enables', 'term': 'glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity'}, 
#                                                                     {'category': 'MF', 'evidence': 'IPI', 'id': 'GO:0005515', 'pubmed': [11724794, 15628863, 16169070, 16799092, 17500595, 17540579, 20029029, 20392205, 20849852, 21044950, 23332158, 23348613, 23355646, 24658140, 25417112, 27387501, 28514442, 29028794, 29128334, 31980649, 32814053, 33961781], 'qualifier': 'enables', 'term': 'protein binding'},
#                                                                     {'category': 'MF', 'evidence': 'ISS', 'id': 'GO:0008017', 'qualifier': 'enables', 'term': 'microtubule binding'}, {'category': 'MF', 'evidence': 'IDA', 'id': 'GO:0019828', 'pubmed': 22832495, 'qualifier': 'enables', 'term': 'aspartic-type endopeptidase inhibitor activity'},
#                                                                     {'category': 'MF', 'evidence': 'ISS', 'id': 'GO:0035605', 'qualifier': 'enables', 'term': 'peptidyl-cysteine S-nitrosylase activity'}, {'category': 'MF', 'evidence': 'IPI', 'id': 'GO:0042802', 'pubmed': [20392205, 21988832], 'qualifier': 'enables', 'term': 'identical protein binding'}, 
#                                                                     {'category': 'MF', 'evidence': 'IEA', 'id': 'GO:0050661', 'qualifier': 'enables', 'term': 'NADP binding'}, 
#                                                                     {'category': 'MF', 'evidence': 'IEA', 'id': 'GO:0051287', 'qualifier': 'enables', 'term': 'NAD binding'}, 
#                                                                     {'category': 'MF', 'evidence': 'IPI', 'id': 'GO:0097718', 'pubmed': 11785981, 'qualifier': 'enables', 'term': 'disordered domain specific binding'}]}}


G = nx.DiGraph()

# Step 2: Add nodes (genes and GO terms) and their relations (edges)
# Example genes and relationships (you can add more)
# Nodes are gene or GO terms, edges are relationships like 'enables', 'located in', 'part of'
genes = {
    'Gene1': {'GO:0008152': 'enables'},
    'Gene2': {'GO:0009058': 'part of'},
    'Gene3': {'GO:0044281': 'located in'},
    'Gene4': {'GO:0019319': 'part of'}
}

# Add nodes and edges to the graph
for gene, go_terms in genes.items():
    for go_term, relationship in go_terms.items():
        G.add_node(gene)
        G.add_node(go_term)
        G.add_edge(gene, go_term, label=relationship)

# Step 3: Draw the graph
pos = nx.spring_layout(G)  # Generate a layout for the graph
labels = nx.get_edge_attributes(G, 'label')  # Get edge labels

# Draw nodes and edges
nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=3000, font_size=10, font_weight='bold', arrows=True)
nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)

# Show the graph
plt.show()
