# Libraries ###########################################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn import preprocessing
import random
from sklearn.metrics import roc_curve, roc_auc_score

# Functions ###########################################################################################################
def apply_z_score(counts_table):
    '''
    Scales RNAseq counts table using z-score
    @param counts_table: dataframe with normalized RNAseq counts
    @return counts_zscore: dataframe with standardized counts
    '''
    columns = list(counts_table.columns)
    index = counts_table.index
    counts_zscore = preprocessing.StandardScaler().fit_transform(counts_table)
    counts_zscore = pd.DataFrame(counts_zscore, columns=columns, index=index)
    return counts_zscore

def normalize_data(counts_zscore):
    '''
    Normalizes RNAseq counts table using min-max normalization
    @param counts_zscore: dataframe with standardized RNAseq counts
    @return counts_zscore_minmax: dataframe with min-max normalized counts
    '''
    columns = list(counts_zscore.columns)
    index = counts_zscore.index
    counts_zscore_minmax = preprocessing.MinMaxScaler().fit_transform(counts_zscore)
    counts_zscore_minmax = pd.DataFrame(counts_zscore_minmax, columns=columns, index=index)
    return counts_zscore_minmax

def getGeneList(tested_genes, n, rank_by='FDR', logFC_filter=True):
    '''
    Crates a gene list of interest
    @param tested_genes: dataframe containing gene symbols and their respective FDR and logFC values
    @param n: number of genes to be returned
    @param rank_by: variable to be used for ranking. Possible values: 'FDR', 'logFC_up', 'logFC_down'
    @param logFC_filter: if True, keeps only genes with (logFC <= -1 | logFC >= 1)
    @return gene_list: a list containing the gene symbols of interest
    '''
    if logFC_filter:
        tested_genes = tested_genes.loc[(tested_genes['logFC'] < -1) | (tested_genes['logFC'] > 1)]
    if rank_by == 'logFC_up':
        ranked_genes = tested_genes[tested_genes['FDR'] <= 0.01]  # keep significant genes only
        ranked_genes = ranked_genes.sort_values(by='logFC', ascending=False)
    elif rank_by == 'logFC_down':
        ranked_genes = tested_genes[tested_genes['FDR'] <= 0.01]  # keep significant genes only
        ranked_genes = ranked_genes.sort_values(by='logFC', ascending=True)
    else:
        ranked_genes = tested_genes.sort_values(by='FDR')
    gene_list = ranked_genes.iloc[:n, 1].tolist()
    return(gene_list)

def prepareScoringData(counts_zscore_minmax, sample_info, genes):
    '''
    Filters counts table using genes list and adds mutation status to each sample
    @param counts_zscore_minmax: dataframe with RNAseq counts (standardized and normalized), with patients as rows and
    genes as columns
    @param sample_info: dataframe containing mutation_status for all samples
    @param genes: list of genes to be considered
    @return scoring_counts: dataframe containing the RNAseq counts for the genes of interest and mutation status for all samples
    '''
    # Filter df_rna using gene_list
    # Add mutation status as last columns
    gene_list = genes
    mutation = sample_info['mutation_status'].tolist()
    scoring_counts = counts_zscore_minmax[counts_zscore_minmax.columns.intersection(gene_list)]
    scoring_counts = scoring_counts.assign(mutation = mutation)
    return(scoring_counts)

def kaufmanScore(scoring_counts):
    '''
    Computes Kaufman score (unweighted mean) for each sample and stores them in a "score" column in the counts dataframe
    @param scoring_counts: dataframe with RNAseq counts for the genes of interest (generated with prepareScoringData() function)
    @return df_score: dataframe containing the Kaufman score computed for each sample
    '''
    # Computes Kaufman score for each sample (arithmetic mean of all gene values)
    ncol = len(scoring_counts.columns)
    score = scoring_counts.iloc[:, 1:ncol-1].apply(np.mean, axis=1)
    df_score = scoring_counts.assign(score = score)
    return(df_score)

def plotScores(df_score, title = None, y_lim=None):
    '''
    Plots Kaufman score for mutant vs wt (boxplot)
    @param df_score: dataframe containing Kaufman score and mutation_status for all samples
    @param title: string, the plot title
    @param y_lim: 2 integers, specifying the min and max values along the plot y-axis
    @return fig: figure containing the plot
    '''
    fig_dims = (4, 6)
    fig, ax = plt.subplots(figsize=fig_dims)
    sns.boxplot(x='mutation', y='score', data=df_score, palette=["darkblue", "red"])
    plt.ylim(y_lim)
    plt.xticks(rotation=0, fontsize=15)
    plt.yticks(fontsize=13)
    plt.xlabel('STK11 mutation', fontsize=15, weight='bold')
    plt.ylabel('K-score', fontsize=15, weight='bold')
    plt.title(title, fontsize=17, weight='bold')
    return(fig)

def plotROC(fpr, tpr, title=None):
    '''
    Plots a ROC curve
    @param fpr: array, the false positive rate
    @param tpr: array, the true positive rate
    @param title: string, the plot title
    @param y_lim: 2 integers, specifying the min and max values along the plot y-axis
    @return fig: figure containing the plot
    '''
    fig_dims = (6, 6)
    fig, ax = plt.subplots(figsize=fig_dims)
    plt.plot(fpr, tpr)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel('FPR', fontsize=15, weight='bold')
    plt.ylabel('TPR', fontsize=15, weight='bold')
    plt.title(title, fontsize=17, weight='bold')
    return(fig)


# Load data ###########################################################################################################
# counts: counts_table with RNAseq counts normalized by TMM in edgeR
# protein-coding genes only, tumor samples only, outlier plate samples excluded
counts = pd.read_csv('data/data_for_ML/luad_rna_TMM_counts_pc_genes_tumor_only_no_outliers.csv')
# sample_info: info about samples (mutation_status)
sample_info = pd.read_csv('data/data_for_ML/sample_info_tumor_only_no_outliers.csv')
# genes: list of gene symbols (protein coding genes only, corresponding to df_rna (order preserved)
genes = pd.read_csv('data/data_for_ML/genes.csv')
# exact_test: table with exact test results (edgeR) for all protein-coding genes
tested_genes = pd.read_csv('results/DE_analysis/ExactTest_pc_genes_tumor_only_no_outliers.csv')

# Pre-process data ####################################################################################################
# Transpose df_rna (to have genes as columns, samples as rows)
counts_table = counts.rename(index = genes['Symbol']) # rename rows as gene symbols
counts_table = counts_table.transpose()

# Standardize data with Z-score
counts_zscore = apply_z_score(counts_table)
counts_zscore.to_csv('data/data_for_ML/TMM_counts_zscore.csv', index=False)

# Scale data using zero-one normalization
counts_zscore_minmax = normalize_data(counts_zscore)
counts_zscore_minmax.to_csv('data/data_for_ML/TMM_counts_zscore_minmax.csv', index=False)

# Kaufman scoring #####################################################################################################
# get list of genes of interest #######################################################################################
# 16-signature genes
kaufman_genes = ['DUSP4', 'PDE4D', 'IRS2', 'BAG1', 'HAL', 'TACC2', 'AVPI1', 'CPS1', 'PTP4A1', 'RFK',
                 'SIK1', 'FGA','GLCE', 'TESC', 'MUC5AC', 'TFF1']

# Ranked gene lists
top_5_FDR = getGeneList(tested_genes, 5, rank_by='FDR', logFC_filter=True)
top_16_FDR = getGeneList(tested_genes, 16, rank_by='FDR', logFC_filter=True)
top_50_FDR = getGeneList(tested_genes, 50, rank_by='FDR', logFC_filter=True)
top_100_FDR = getGeneList(tested_genes, 100, rank_by='FDR', logFC_filter=True)
top_1000_FDR = getGeneList(tested_genes, 1000, rank_by='FDR', logFC_filter=True)

top_16_logFC_up = getGeneList(tested_genes, 16, rank_by='logFC_up', logFC_filter=True)
top_100_logFC_up = getGeneList(tested_genes, 100, rank_by='logFC_up', logFC_filter=True)

top_16_logFC_down = getGeneList(tested_genes, 16, rank_by='logFC_down', logFC_filter=True)
top_100_logFC_down = getGeneList(tested_genes, 100, rank_by='logFC_down', logFC_filter=True)

top_5_important = ['NPY', 'EYS', 'MTMR7', 'ODC1', 'SLC16A14']
top_16_important = ['NPY', 'EYS', 'MTMR7', 'ODC1', 'SLC16A14', 'NNAT', 'ZACN', 'CACNB2', 'GREB1', 'HPX',
                    'FURIN', 'TENM1', 'ASPG', 'RPH3AL', 'RET', 'CRB1']

# get random gene list
n = random.sample(range(0, len(tested_genes.index)), 16) # generate 16 random indexes
random_genes = tested_genes.iloc[n, 1]

# Scoring ############################################################################################################
# Prepare data for scoring (select relevant genes and get mutation status)
scoring_counts = prepareScoringData(counts_zscore_minmax, sample_info, top_16_FDR)

# Score samples and plot data
df_score = kaufmanScore(scoring_counts)
fig = plotScores(df_score, title='K-score[top 16 FDR genes]', y_lim=None)
fig.savefig('results/kaufman_score/kscore_top_16_FDR.png', bbox_inches='tight')

# ROC curves to evaluate scoring performance
mutation = df_score['mutation'].tolist()
score = df_score['score'].tolist()
fpr, tpr, thresh = roc_curve(mutation, score, pos_label=1)
fig = plotROC(fpr, tpr, title='K-score[16 random genes]')
fig.savefig('results/kaufman_score/k_score_top_16_FDR_ROC.png')
auc_score = roc_auc_score(mutation, score)

# randomized control - shuffle labels
random.shuffle(mutation)
fpr_rand, tpr_rand, thresh_rand = roc_curve(mutation, score, pos_label=1)
fig = plotROC(fpr_rand, tpr_rand)
fig.savefig('results/kaufman_score/k_score_randomized_labels_ROC.png')
auc_score_random = roc_auc_score(mutation, score)

# Misc ###############################################################################################################
# check intersections between gene lists
np.intersect1d(top_16_FDR, top_16_logFC_up) # ['BLOC1S5-TXNDC5', 'INHA', 'NNAT', 'NPY'] - 4 genes overlap

np.intersect1d(top_100_FDR, top_100_logFC_up)
# 48 genes overlap
# ['ASPG', 'BAALC', 'BLOC1S5-TXNDC5', 'BMP6', 'BPIFA2', 'BPIFB2',
#  'CALCA', 'CALCB', 'CBR1', 'CCDC154', 'CHGB', 'CHRDL2', 'CRB1',
#  'EYS', 'FAM163A', 'FAM9B', 'FGL1', 'GALNTL6', 'GLTPD2', 'GREB1',
#  'HPX', 'INHA', 'KLRC2', 'LCN15', 'LILRA2', 'LRRC26', 'MAFA',
#  'MARCHF4', 'MTMR7', 'NNAT', 'NPAS3', 'NPY', 'NTNG2', 'ODC1',
#  'PAK3', 'PTPRN', 'RELN', 'RERGL', 'RET', 'SLC16A14', 'SLC26A4',
#  'SLC7A2', 'SPP2', 'SRARP', 'TAF1L', 'TENM1', 'TXNDC2', 'ZACN']

np.intersect1d(top_16_FDR, top_16_logFC_down) # no overlap

np.intersect1d(top_100_FDR, top_100_logFC_down)
# 7 genes overlap
# ['ADGRF1', 'ALB', 'FAT2', 'GJA3', 'IVL', 'SHISA3', 'SLC6A20']


