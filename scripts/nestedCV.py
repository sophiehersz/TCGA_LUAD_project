import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_validate

# Functions #########################################################################################################
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
    @return: a list containing the gene symbols of interest
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

def getRandomGeneList(tested_genes, n):
    '''
    Returns a random gene list, as well as a dataframe containing statistics (p-value, FDR...) for the selected genes
    @param tested_genes: dataframe containing gene symbols and their respective FDR and logFC values
    @param n: number of genes to be returned
    @return random genes: a list of random gene symbols
    @return random_genes_stats: a dataframe containing statistics (p-value, FDR...) for the selected genes
    '''
    indexes = random.sample(range(0, len(tested_genes.index)), n) # generate n random indexes
    random_genes_stats = tested_genes.iloc[indexes, :]
    random_genes = tested_genes.iloc[indexes, 1]
    return(random_genes, random_genes_stats)

def prepareData(counts_zscore_minmax, sample_info, gene_list):
    '''
    Filters df_rna using gene_list and creates X (input variables) and y (labels) for ML
    @param counts_zscore_minmax: dataframe with normalized and standardized gene counts for all samples
    @oaram sample_info: dataframe containing sample barcodes and their corresponding mutation status
    @param gene_list: list of the gene symbols of interest
    @return X: array of input variables
    @return y: array of labels
    '''
    y = sample_info['mutation_status'].values
    X = counts_zscore_minmax[counts_zscore_minmax.columns.intersection(gene_list)].values
    return(X, y)

def nestedCV(model, params, counts_zscore_minmax, sample_info, gene_lists):
    '''
    Performs nested cross-validation to optimise model hyperparameters and evaluate performance
    @param model: the model to be trained
    @param params: hyperparameter search space
    @param counts_zscore_minmax: dataframe with normalized and standardized gene counts for all samples
    @param sample_info: dataframe containing sample barcodes and their corresponding mutation status
    @param gene_lists: dictionary containing the gene_list names and corresponding gene symbols
    @return cv_results: dictionary with full cross-validation results for each gene_list (as keys)
    @return df_scores: dataframe with test scores for all outer folds
    @return average_scores: dataframe with mean and std of scores for all gene_lists
    '''
    # configure the cross-validation procedure
    cv_inner = StratifiedKFold(n_splits=5, random_state=0, shuffle=True)
    cv_outer = StratifiedKFold(n_splits=10, random_state=0, shuffle=True)

    # define search
    search = GridSearchCV(model, params,
                          scoring='roc_auc',
                          return_train_score=True,
                          cv=cv_inner,
                          refit=True,
                          n_jobs=-1)

    # create empty dictionary and dataframe to store results
    cv_results = dict.fromkeys(gene_lists.keys(), None)
    scores = pd.DataFrame(columns=gene_lists.keys())

    for i, k in enumerate(gene_lists.keys()):
        # prepare X and y using the list of genes of interest
        X, y = prepareData(counts_zscore_minmax, sample_info, gene_lists[k])

        # execute the nested cross-validation and get results
        cv_out = cross_validate(search, X, y,
                                scoring='roc_auc',
                                cv=cv_outer,
                                n_jobs=-1,
                                return_train_score=True,
                                return_estimator=True)
        cv_results[k] = cv_out
        scores[k] = cv_out['test_score']
        print('gene list ' + str(i + 1) + ' of ' + str(len(gene_lists.keys())))

    # compute average scores and std
    average_scores = pd.DataFrame(index=gene_lists.keys(), columns=['average_roc_auc', 'std'])
    average_scores['average_roc_auc'] = scores.apply(np.mean, axis='index')
    average_scores['std'] = scores.apply(np.std, axis='index')

    return (cv_results, scores, average_scores)

def plot_nested_cv_results(scores, title=None, y_lim=None):
    '''
    Plots nested cv results (boxplot)
    @param scores: dictionary containing nested cv score for all gene lists of interest
    @param title: string, the plot title
    @param y_lim: 2 integers, specifying the min and max values along the plot y-axis
    @return fig: figure containing the plot
    '''
    fig, ax = plt.subplots()
    sns.boxplot(data=scores, palette='Greys')
    sns.despine()
    plt.ylim(y_lim)
    plt.xticks(rotation=90, fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel('Gene sets', fontsize=15, weight='bold')
    plt.ylabel('ROC-AUC score', fontsize=15, weight='bold')
    plt.title(title, fontsize=17, weight='bold')
    return(fig)

def getEstimatorParams(cv_results, gene_list_name, params):
    '''
    Returns dataframe with hyperparameters for the best estimators after cross-validation
    @param cv_results: dictionary containing nested cv results, output from nestedCV() function
    @param gene_list_name: the name of the gene_list of interest (must be a cv_results key)
    @param params: list containing the names of the hyperparameters that were optimised. e.g. ['n_estimators', 'min_samples_split', 'max_depth']
    @return best_params: dataframe containing best values for parameters of interest
    '''
    k = len(cv_results[gene_list_name]['estimator'])  # find number of outer K-folds
    best_params = pd.DataFrame(columns=range(len(params)))
    for i in range(k):
        all_params = cv_results[gene_list_name]['estimator'][i].best_estimator_.get_params()
        values = [all_params[x] for x in params]
        best_params = best_params.append(np.reshape(values, (1, -1)).tolist())
    best_params.columns = params
    return (best_params)

def getFeatureImportance(cv_results, gene_list_name, gene_lists, model_type):
    '''
    Creates dataframe with feature importance from the best estimators for each outer fold
    @param cv_results: dictionary containing nested cv results, output from nestedCV() function
    @param gene_list_name: name of the gene_list of interest (must be a cv_results key)
    @param gene_lists: dictionary containing the gene_list names and corresponding gene symbols
    @param model_type: 'LR' for Logistic regression or 'RF' for Random Forest
    @return feat_importances: a dataframe containing feature importance values for each gene (columns) and outer runs (indexes)
    '''
    k = len(cv_results[gene_list_name]['estimator'])  # find number of outer K-folds
    n = len(gene_lists[gene_list_name])  # find number of genes
    feat_importances = pd.DataFrame(columns=range(n))

    if model_type == 'RF':
        for i in range(k):
            features = cv_results[gene_list_name]['estimator'][i].best_estimator_.feature_importances_
            feat_importances = feat_importances.append(np.reshape(features, (1, -1)).tolist())
        feat_importances.columns = gene_lists[gene_list_name]
    elif model_type == 'LR':
        for i in range(k):
            features = cv_results[gene_list_name]['estimator'][i].best_estimator_.coef_
            feat_importances = feat_importances.append(features.tolist())
        feat_importances.columns = gene_lists[gene_list_name]

    return (feat_importances)

def plotFeatureImportance(feat_importances, title=None, y_lim=None):
    '''
    Plots feature importances (boxplot)
    @param feat_importances: dataframe containing the feature importance values for each gene, at each outer run
    @param title: string, the plot title
    @param y_lim: 2 integers, specifying the min and max values along the plot y-axis
    @return fig: figure containing the plot
    '''
    fig, ax = plt.subplots()
    sns.boxplot(data=feat_importances, palette='Greys')
    sns.despine()
    plt.ylim(y_lim)
    plt.xticks(rotation=90, fontsize=12)
    plt.ylabel('Values', fontsize=12, weight='bold')
    plt.title(title, fontsize=15, weight='bold')
    return(fig)

# Load data ##########################################################################################################
counts = pd.read_csv('data/data_for_ML/luad_rna_TMM_counts_pc_genes_tumor_only_no_outliers.csv') # counts_table with RNAseq counts normalized by TMM in edgeR
genes = pd.read_csv('data/data_for_ML/genes.csv') # list of gene symbols, corresponding to counts table (order preserved)
sample_info = pd.read_csv('data/data_for_ML/sample_info_tumor_only_no_outliers.csv') # dataframe with info about samples (mutation_status)
tested_genes = pd.read_csv('results/DE_analysis/ExactTest_pc_genes_tumor_only_no_outliers.csv') # dataframe with exact test results (edgeR)

# Pre-process data ####################################################################################################
# Transpose counts table  (to have genes as columns, samples as rows)
counts_table = counts.rename(index = genes['Symbol']) # rename rows as gene symbols
counts_table = counts_table.transpose()
# Standardize data with Z-score
counts_zscore = apply_z_score(counts_table)
counts_zscore.to_csv('data/data_for_ML/TMM_counts_zscore.csv', index=False)
# Scale data using zero-one normalization
counts_zscore_minmax = normalize_data(counts_zscore)
counts_zscore_minmax.to_csv('data/data_for_ML/TMM_counts_zscore_minmax.csv', index=False)

# Make gene lists ####################################################################################################
# 16-signature genes
kaufman_genes = ['DUSP4', 'PDE4D', 'IRS2', 'BAG1', 'HAL', 'TACC2', 'AVPI1', 'CPS1', 'PTP4A1', 'RFK',
                 'SIK1', 'FGA','GLCE', 'TESC', 'MUC5AC', 'TFF1']

# Ranked gene lists
kaufman_genes = ['DUSP4', 'PDE4D', 'IRS2', 'BAG1', 'HAL', 'TACC2', 'AVPI1', 'CPS1', 'PTP4A1', 'RFK',
                 'SIK1', 'FGA','GLCE', 'TESC', 'MUC5AC', 'TFF1']

# Ranked gene lists
top_5_FDR = getGeneList(tested_genes, 5, rank_by='FDR', logFC_filter=True)
top_16_FDR = getGeneList(tested_genes, 16, rank_by='FDR', logFC_filter=True)
top_50_FDR = getGeneList(tested_genes, 50, rank_by='FDR', logFC_filter=True)
top_100_FDR = getGeneList(tested_genes, 100, rank_by='FDR', logFC_filter=True)
top_1000_FDR = getGeneList(tested_genes, 1000, rank_by='FDR', logFC_filter=True)

# Non-significant genes
ranked_genes = tested_genes.sort_values(by='FDR', ascending=False)
bottom_16_FDR = ranked_genes.iloc[:16, 1].tolist()
bottom_100_FDR = ranked_genes.iloc[:100, 1].tolist()

# Random gene sets
random_genes_1, random_genes_stats_1 = getRandomGeneList(tested_genes, 16)
random_genes_2, random_genes_stats_2 = getRandomGeneList(tested_genes, 16)
random_genes_3, random_genes_stats_3 = getRandomGeneList(tested_genes, 16)

# Genes selected based on RF feature importance
top_5_important = ['NPY', 'EYS', 'MTMR7', 'ODC1', 'SLC16A14']
top_16_important = ['NPY', 'EYS', 'MTMR7', 'ODC1', 'SLC16A14', 'NNAT', 'ZACN', 'CACNB2', 'GREB1', 'HPX',
                    'FURIN', 'TENM1', 'ASPG', 'RPH3AL', 'RET', 'CRB1']

# prepare gene_lists object
gene_lists = {'Kaufman':kaufman_genes,
              'Top 5 DE':top_5_FDR,
              'Top 16 DE':top_16_FDR,
              'Top 50 DE':top_50_FDR,
              'Top 100 DE':top_100_FDR,
              'Top 1000 DE':top_1000_FDR,
              'Bottom 16 NS':bottom_16_FDR,
              'Bottom 100 NS':bottom_100_FDR,
              '16 random genes 1':random_genes_1,
              '16 random genes 2':random_genes_2,
              '16 random genes 3':random_genes_3}

important_gene_list = {'Kaufman': kaufman_genes,
                       'Top 5 DE':top_5_FDR,
                       'Top 5 important': top_5_important,
                       'Top 16 DE':top_16_FDR,
                       'Top 16 important': top_16_important,
                       'Top 50 DE':top_50_FDR}

# Logistic regression ##############################################################################################
# define the model
LR_model = LogisticRegression(solver='liblinear',
                              class_weight='balanced',
                              random_state=0)

# define hyperparameter search space
LR_params = {'C': [0.01, 0.1, 1, 10, 100],
             'penalty': ['l1', 'l2']}

# run nested CV
LR_cv_results, LR_cv_scores, LR_cv_average_scores = nestedCV(LR_model, LR_params, counts_zscore_minmax, sample_info, important_gene_list)

# plot cv scores
fig = plot_nested_cv_results(LR_cv_scores,
                                 y_lim=[0.2, 1.1],
                                 title='Logistic Regression\nnested cross-validation (k=10)')
fig.savefig('results/DE_analysis/logistic_regression/LR_nested_cv_important_genes.png', bbox_inches='tight')

# get best estimators parameters for each outer run
best_params = getEstimatorParams(LR_cv_results,
                              gene_list_name='Top 50 DE',
                              params=['C', 'penalty'])

# Get coefficients of best estimators for all outer runs
LR_coefficients = getFeatureImportance(LR_cv_results, gene_list_name='Top 50 DE', gene_lists=gene_lists, model_type='LR')

# Plot results
fig = plotFeatureImportance(LR_coefficients, y_lim=[-5,20], title='LR coefficients\ntop 16 FDR genes')
fig.savefig('results/DE_analysis/logistic_regression/LR_coefficients_top_50_FDR.png', bbox_inches='tight')


# Random Forest ##############################################################################################
# define the model
RF_model = RandomForestClassifier(class_weight='balanced',
                                  random_state=0)

# define hyperparameter search space
RF_params = {'n_estimators': [50, 100],
             'max_depth': [2, 3, 4]}

# run nested CV
RF_cv_results, RF_cv_scores, RF_cv_average_scores = nestedCV(RF_model, RF_params, counts_zscore_minmax, sample_info, important_gene_list)

# plot cv scores
fig = plot_nested_cv_results(RF_cv_scores,
                                 y_lim=[0.2, 1.1],
                                 title='Random Forest\nnested cross-validation (k=10)')
fig.savefig('results/DE_analysis/random_forest/RF_nested_cv_important_genes.png', bbox_inches='tight')

# get best estimators hyperparameters for each outer run
RF_param_df = getEstimatorParams(RF_cv_results,
                                 gene_list_name='Top 50 DE',
                                 params=['n_estimators', 'max_depth'])

# Get feature importances of best estimators for all outer runs
RF_feat_importance = getFeatureImportance(RF_cv_results, gene_list_name='Top 50 DE', gene_lists=gene_lists, model_type='RF')

# Plot results
fig = plotFeatureImportance(RF_feat_importance, y_lim=None, title='Random Forest feature importance\nTop 50 DE genes')
fig.savefig('results/DE_analysis/random_forest/RF_feature_importance_top_50_DE_genes.png', bbox_inches='tight')