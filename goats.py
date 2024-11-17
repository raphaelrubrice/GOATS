import mygene as mg

mg = mg.MyGeneInfo()

def get_term_ids(go_list, numeric=False):
    if numeric:
        return [int(term['id'][3:]) for term in go_list]
    return [term['id'] for term in go_list]

def get_term_descriptions(go_list):
    return [term['term'] for term in go_list]

def get_gene_predicate(go_list):
    return [term['qualifier'] for term in go_list]

def get_infos(geneset):
    info_dico = {gene:{} for gene in geneset}
    out = mg.querymany(geneset, scopes='symbol,entrezgene,ensemblgene', species=9606, fields=['entrezgene','go'])
    for i in range(len(out)):
        if 'notfound' in out[i].keys():
            gene = out[i]['query']
            print(f'\n{gene} not found. Removing {gene}.')
            geneset.remove(gene)
    for i in range(len(geneset)):
        goterms = out[i]['go']
        go_mf = goterms['MF'] if 'MF' in goterms.keys() else None
        mf_go_ids = get_term_ids(go_mf, numeric=True)
        mf_go_desc = get_term_descriptions(go_mf)
        mf_go_pred = get_gene_predicate(go_mf)

        go_bp = goterms['BP'] if 'BP' in goterms.keys() else None
        bp_go_ids = get_term_ids(go_bp, numeric=True)
        bp_go_desc = get_term_descriptions(go_bp)
        bp_go_pred = get_gene_predicate(go_bp)

        go_cc = goterms['CC'] if 'CC' in goterms.keys() else None
        cc_go_ids = get_term_ids(go_cc, numeric=True)
        cc_go_desc = get_term_descriptions(go_cc)
        cc_go_pred = get_gene_predicate(go_cc)

        # print("\nMF GO terms:", mf_go_desc)
        # print("\nBP GO terms:", bp_go_desc)
        # print("\nCC GO terms:", cc_go_desc)

        # print("\nMF GO ids:", mf_go_ids)
        # print("\nBP GO ids:", bp_go_ids)
        # print("\nCC GO ids:", cc_go_ids)

        gene_dico = {'ids':{'MF':mf_go_ids, 'BP':bp_go_ids, 'CC':cc_go_ids},
                     'terms':{'MF':mf_go_desc, 'BP':bp_go_desc, 'CC':cc_go_desc},
                     'predicate':{'MF':mf_go_pred, 'BP':bp_go_pred, 'CC':cc_go_pred}}
        info_dico[geneset[i]].update(gene_dico)
    return info_dico

# def common_ids(infos):
#     dico = {'MF':[], 'BP':[], 'CC':[]}
#     for section in ['MF', 'BP', 'CC']:
#         infos[gene]['ids'][section]

def show_relations(infos, gene, format='str', query='predicate', show=True):
    out = {gene:[]} if format=='str' else {}
    if show:
        print(f'\n{gene} :')
    for section in ['MF', 'BP', 'CC']:
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
                    print(f'\n{f"'GO:{pad+str(go_term)}'"} : {go_relation}')
                else:
                    print(f'\n{go_term} : {go_relation}')
    return out

def common_sets(geneset, collection=None):
    assert (collection == None) or ('list' in str(type(collection))), "Collection must be None or a List object."
    infos = get_infos(geneset)
    go_collections = infos[geneset[0]]['ids'].keys()

    if not collection:
        collection = go_collections
    in_common = {}
    for c in collection:
        common_goterms = set([])
        gene_roles = {}
        for gene in infos.keys():
            goterms = set(infos[gene]['ids'][c]) # list of all go term IDS for this gene for this collection
            predicates = show_relations(infos, gene, format='int', query='predicate', show=False)

            if len(common_goterms) == 0:
                # On initialization
                common_goterms = goterms
            else:
                # we update the set to be only IDS which are in common (INTERSECT)
                common_goterms = common_goterms.intersection(goterms)

        for gene in infos.keys():
            gene_roles[gene] = {go_id:predicate for (go_id,predicate) in predicates.items() if go_id in common_goterms}
        descriptions = show_relations(infos, geneset[0], format='int', query='terms', show=False)
        
        term_descriptions = [descriptions[go_id] for go_id in common_goterms]
        in_common[c] = {'ids':common_goterms, 'predicates':gene_roles, 'descriptions':term_descriptions}
    return in_common

def get_role(in_common, collection, gene, go_id):
    return in_common[collection]['predicates'][gene][go_id]

def show_description(in_common, go_id_list=None):
    if go_id_list:
        for go_id in go_id_list:
            for collection in in_common.keys():
                if go_id in in_common[collection]['ids']:
                    i = 0
                    for item in in_common[collection]['ids']:
                        if item == go_id:
                            break
                        i += 1
                    desc = in_common[collection]['descriptions'][i]
                    print(f"\n'{go_id}' : {desc}")
                    for gene in in_common[collection]['predicates'].keys():
                        print(f"{gene} -> {get_role(in_common, collection, gene, go_id)} -> {desc}")
                    break
    else:
        memory = set([])
        for collection in in_common.keys():
            go_id_list = set(list(in_common[collection]['ids']))
            for go_id in go_id_list:
                if go_id in in_common[collection]['ids'] and go_id not in memory:
                    i = 0
                    for item in in_common[collection]['ids']:
                        if item == go_id:
                            break
                        i += 1
                    desc = in_common[collection]['descriptions'][i]
                    print(f"\n'{go_id}' : {desc}")
                    for gene in in_common[collection]['predicates'].keys():
                        print(f"{gene} -> {get_role(in_common, collection, gene, go_id)} -> {desc}")
                    memory.add(go_id)
                    

geneset = ['IL6', 'IL2', 'IFN', 'bonjour']

infos = get_infos(geneset)
#print(infos)
out = show_relations(infos, geneset[0], format='int', query='terms', show=False)
#print(out)
in_common = common_sets(geneset)
print(in_common)
print('\n\n')
show_description(in_common)


