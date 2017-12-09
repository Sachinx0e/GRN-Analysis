% Ecoli Genes = 5, Interactions = 4 TimeSeries = 25
mRNA1 = dlmread('data/bgrmi/original/ecoli_5_genes_4_interactions_25_time_series/Ecoli_5_genes_4_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_5_genes_4_interactions_25_time_series/Ecoli_5_genes_4_interactions_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)

% Ecoli Genes = 5, Interactions = 4 TimeSeries = 50
mRNA1 = dlmread('data/bgrmi/original/ecoli_5_genes_4_interactions_50_time_series/Ecoli_5_genes_4_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_5_genes_4_interactions_50_time_series/Ecoli_5_genes_4_interactions_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)