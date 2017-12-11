% Ecoli Genes = 5, Interactions = 4 TimeSeries = 25
mRNA1 = dlmread('data/bgrmi/original/ecoli_5_genes_4_interactions_25_time_series/Ecoli_5_genes_4_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_5_genes_4_interactions_25_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)


% Ecoli Genes = 5, Interactions = 4 TimeSeries = 50
mRNA1 = dlmread('data/bgrmi/original/ecoli_5_genes_4_interactions_50_time_series/Ecoli_5_genes_4_interactions_dream4_timeseries.tsv');
No_Replicates = 50;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_5_genes_4_interactions_50_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)


% Ecoli Genes = 5, Interactions = 7 TimeSeries = 25
mRNA1 = dlmread('data/bgrmi/original/ecoli_5_genes_7_interactions_25_time_series/Ecoli_5_genes_7_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_5_genes_7_interactions_25_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)


% Ecoli Genes = 5, Interactions = 7 TimeSeries = 50
mRNA1 = dlmread('data/bgrmi/original/ecoli_5_genes_7_interactions_50_time_series/Ecoli_5_genes_7_interactions_dream4_timeseries.tsv');
No_Replicates = 50;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_5_genes_7_interactions_50_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)

% Ecoli Genes = 10, Interactions = 14 TimeSeries = 25
mRNA1 = dlmread('data/bgrmi/original/ecoli_10_genes_14_interactions_25_time_series/Ecoli_10_genes_14_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_10_genes_14_interactions_25_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)


% Ecoli Genes = 10, Interactions = 14 TimeSeries = 50
mRNA1 = dlmread('data/bgrmi/original/ecoli_10_genes_14_interactions_50_time_series/Ecoli_10_genes_14_interactions_dream4_timeseries.tsv');
No_Replicates = 50;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_10_genes_14_interactions_50_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)


% Ecoli Genes = 10, Interactions = 20 TimeSeries = 25
mRNA1 = dlmread('data/bgrmi/original/ecoli_10_genes_20_interactions_25_time_series/Ecoli_10_genes_20_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_10_genes_20_interactions_25_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)


% Ecoli Genes = 10, Interactions = 20 TimeSeries = 50
mRNA1 = dlmread('data/bgrmi/original/ecoli_10_genes_20_interactions_50_time_series/Ecoli_10_genes_20_interactions_dream4_timeseries.tsv');
No_Replicates = 50;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_10_genes_20_interactions_50_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)


% Ecoli Genes = 20, Interactions = 24 TimeSeries = 25
mRNA1 = dlmread('data/bgrmi/original/ecoli_20_genes_24_interactions_25_time_series/Ecoli_20_genes_24_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_20_genes_24_interactions_25_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)

% Ecoli Genes = 20, Interactions = 24 TimeSeries = 50
mRNA1 = dlmread('data/bgrmi/original/ecoli_20_genes_24_interactions_50_time_series/Ecoli_20_genes_24_interactions_dream4_timeseries.tsv');
No_Replicates = 50;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_20_genes_24_interactions_50_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)

% Ecoli Genes = 20, Interactions = 35 TimeSeries = 25
mRNA1 = dlmread('data/bgrmi/original/ecoli_20_genes_35_interactions_25_time_series/Ecoli_20_genes_35_interactions_dream4_timeseries.tsv');
No_Replicates = 25;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_20_genes_35_interactions_25_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)

% Ecoli Genes = 20, Interactions = 35 TimeSeries = 50
mRNA1 = dlmread('data/bgrmi/original/ecoli_20_genes_35_interactions_50_time_series/Ecoli_20_genes_35_interactions_dream4_timeseries.tsv');
No_Replicates = 50;
No_Time_Points = 21;
dt = 50/60;
dt = repmat(dt, No_Time_Points-1,1);
[Adjacency_Matrix]=BGRMI(mRNA1, dt, No_Time_Points, No_Replicates);
file_name = 'data/bgrmi/predicted/ecoli_20_genes_35_interactions_50_time_series_dream4_adjacency.tsv';
dlmwrite(file_name,Adjacency_Matrix,'delimiter','\t','precision',10)