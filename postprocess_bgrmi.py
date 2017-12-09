import pandas

def load_gene_list(file_name):
    dataframe = pandas.read_csv(
        file_name,
        sep='\t', header=0)
    return list(dataframe.columns.values)[:]

def load_adjacency_matrix(file_name):
    datatframe = pandas.read_csv(
        file_name,
        sep='\t', header=0)
    return datatframe.as_matrix()

def serialize_predictions(adjacency_matrix,gene_list):
    gene_count = len(gene_list)

    source_list = []
    target_list = []
    confidence_list = []

    for source in range(0,gene_count):
        source_gene = gene_list[source]

        for target in range(0,gene_count):
            target_gene = gene_list[target]

            weight = adjacency_matrix[source,target]

            if target_gene != source_gene:
                source_list.append(source_gene)
                target_list.append(target_gene)
                confidence_list.append(weight)

    df = pandas.DataFrame({'source':source_list,'target':target_list,'weight':confidence_list})
    df_sorted = df.sort_values('weight',ascending=False)
    return df_sorted


if __name__ == "__main__":

    # Ecoli 5 genes 2 regulators 4 interactions 25 timepoints
    gene_list = load_gene_list('data/bgrmi/original/ecoli_5_genes_4_interactions_25_time_series/Ecoli_5_genes_4_interactions_wildtype.tsv')
    adjacency_matrix = load_adjacency_matrix('data/bgrmi/predicted/ecoli_5_genes_4_interactions_25_time_series/Ecoli_5_genes_4_interactions_dream4_adjacency.tsv')
    serialized_predictions = serialize_predictions(adjacency_matrix,gene_list)
    serialized_predictions.to_csv('data/bgrmi/predicted/ecoli_5_genes_4_interactions_25_time_series/Ecoli_5_genes_4_interactions_dream4_ouput.txt',
                                  header=False, sep='\t', index=False)

    # Ecoli 5 genes 2 regulators 4 interactions 25 timepoints
    gene_list = load_gene_list('data/bgrmi/original/ecoli_5_genes_4_interactions_50_time_series/Ecoli_5_genes_4_interactions_wildtype.tsv')
    adjacency_matrix = load_adjacency_matrix('data/bgrmi/predicted/ecoli_5_genes_4_interactions_50_time_series/Ecoli_5_genes_4_interactions_dream4_adjacency.tsv')
    serialized_predictions = serialize_predictions(adjacency_matrix,gene_list)
    serialized_predictions.to_csv('data/bgrmi/predicted/ecoli_5_genes_4_interactions_50_time_series/Ecoli_5_genes_4_interactions_dream4_ouput.txt',
                                  header=False, sep='\t', index=False)



