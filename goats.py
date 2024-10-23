import mygene as mg

mg = mg.MyGeneInfo()

geneset = ['FOXA1', 'IL6']
out = mg.querymany(geneset, scopes='symbol,entrezgene,ensemblgene', species=9606, fields=['entrezgene','go'])

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


def show_relations(infos, gene):
    out = {}
    for section in ['MF', 'BP', 'CC']:
        N = len(infos[gene]['ids'][section])
        for i in range(N):
            print(f'GO:{infos[gene]['ids'][section][i]} : {infos[gene]['predicate'][section][i]}')
            out.update({gene:{go_term:}})

infos = get_infos(geneset)
show_relations(infos, 'FOXA1')


