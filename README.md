# **GOATS**
Gene Ontology Associated Term Similarities.

## Usage:  
To find gene terms shared by genes in a geneset.  
This require 3 files :
- A text file containing the Gene set to query. The format is 1 gene per line.
- A .OBO file (contains the Gene Ontology database). You can download this file here [Current Ontology](https://current.geneontology.org/ontology/index.html)
- A .GAF file for your organism (contains the annotations for genes and GO terms for a specific Taxon). You can find such files here [GAF Files Database](https://current.geneontology.org/products/pages/downloads.html) or, if your organism is less common, find informations about where to find it on the web page [Download Annotations](https://geneontology.org/docs/download-go-annotations/).  

**Required arguments**:
- `-g`, `--geneset` Path to geneset file. ONE GENE SYMBOL PER LINE.
- `-obo`, `--obo`   Path to the .obo ontology file
- `-gaf`, `--gaf`   Path to the .gaf annotation file   

**Optional arguments**:  
- `-threshold`, `--threshold`   Threshold value for filtering -> [0,1] (default: 0.5)
- `-sd`, `--show_descriptions`  To show description object in terminal -> 0 or 1 (default: 0)
- `-save`, `--save`             Save the outputs -> 0 or 1 (default: 1)  

**Example**:  
```bash
. run.sh -g geneset.txt -obo go-basic.obo -gaf goa_human.gaf -threshold 0.7
```  

