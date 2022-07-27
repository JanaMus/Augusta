# RNA-Seq count table normalization
def normalize_CountTable (countTable, gene_lengths, normalization_type):

    if normalization_type == 'RPKM':
        from bioinfokit.analys import norm
        countTable['gl'] = gene_lengths
        nm = norm()
        nm.rpkm(df=countTable, gl='gl')
        normalized_CountTable = nm.rpkm_norm

    elif normalization_type == 'CPM':
        from bioinfokit.analys import norm
        nm = norm()
        nm.cpm(df=countTable)
        normalized_CountTable = nm.cpm_norm

    elif normalization_type == 'TPM':
        from bioinfokit.analys import norm
        countTable['gl'] = gene_lengths
        nm = norm()
        nm.tpm(df=countTable, gl='gl')
        normalized_CountTable = nm.tpm_norm

    else:
        normalized_CountTable = countTable.copy()
    print('Count table normalization done.')
    return normalized_CountTable
