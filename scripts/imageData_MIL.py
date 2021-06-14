import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import train_test_split
from mil.trainer.trainer import Trainer
from mil.metrics import AUC
from mil.models import AttentionDeepPoolingMil
from mil.utils.padding import Padding
from mil.validators import StratifiedKFold
import matplotlib.pyplot as plt
import os
from sklearn.metrics import roc_curve, roc_auc_score


# Functions #########################################################################################################
def pickleSave(object, file_name):
    '''
    Saves object in pickle format
    @param object: the object to be saved
    @param file_name: the output path
    '''
    open_file = open(file_name, "wb")
    pickle.dump(object, open_file)
    open_file.close()

def pickleOpen(file_name):
    '''
    Opens pickle object
    @param file_name: the file's path
    '''
    open_file = open(file_name, "rb")
    loaded_list = pickle.load(open_file)
    open_file.close()
    return(loaded_list)

def filterWSI(wsi_list, sample_info):
    '''
    Creates dataframe with WSI file names containing only patients in patient_list. Also adds mutation status to each sample.
    @param wsi_list: dataframe containing the paths to each sample bag
    @param sample_info: dataframe with patient IDs to be considered and their corresponding mutation status
    @return wsi_filtered: dataframe containing file names and mutation status for all embedding bags belonging to patients in patient_list
    '''
    wsi_list['cases.submitter_id'] = wsi_list['0'].str[20:32]  # adds column with patient ID
    wsi_list['file_name'] = wsi_list['0'].str[20:]  # adds column with embeddings file names
    patient_list = sample_info['cases.submitter_id'].unique()  # gets list of patient ID
    wsi_filtered = wsi_list.loc[wsi_list['cases.submitter_id'].isin(patient_list)]  # keep only patients in patients_list
    wsi_filtered = pd.merge(wsi_filtered, sample_info[['cases.submitter_id', 'mutation_status']].drop_duplicates(),
                            how='left', on='cases.submitter_id')  # adds mutation status to wsi samples
    wsi_filtered = wsi_filtered.rename(columns={'mutation_status': 'label'})
    return (wsi_filtered)

def balanceData(X_train, y_train):
    '''
    Selects subset of training data to get balanced classes
    @param X_train: the training bags
    @param y_train: the training bag labels
    @return: X_train_bal: balanced training bags
    @return y_train_bal: labels corresponding to X_train_bal
    '''
    # get indexes of mutant and wt samples
    all_mut = [idx for idx, element in enumerate(y_train) if element == 1]
    all_wt = [idx for idx, element in enumerate(y_train) if element == 0]
    n_samples = min(len(all_mut), len(all_wt)) # number of samples to be selected

    idx_mut = np.random.choice(all_mut, size=n_samples, replace=False)  # select random unique indexes
    idx_wt = np.random.choice(all_wt, size=n_samples, replace=False)
    idx = np.concatenate((idx_mut, idx_wt))
    np.random.shuffle(idx) # shuffle indexes

    X_train_bal = X_train[idx]
    y_train_bal = y_train[idx]

    return(X_train_bal, y_train_bal)

def subsetData(X_train, y_train, n_samples):
    '''
    Selects subset of training data
    @param X_train: the training bags
    @param y_train: the training bag labels
    @param n_samples: number of samples to be returned (for each class)
    @return: X_train_subset: subset of training bags
    @return y_train_subset: corresponding labels
    '''
    all_class_0 = [idx for idx, element in enumerate(y_train) if element == 0]
    all_class_1 = [idx for idx, element in enumerate(y_train) if element == 1]

    idx_class_0 = np.random.choice(all_class_0, size=n_samples, replace=False)  # select random unique indexes
    idx_class_1 = np.random.choice(all_class_1, size=n_samples, replace=False)

    idx = np.concatenate((idx_class_0, idx_class_1))
    np.random.shuffle(idx) # shuffle indexes

    X_train_subset = X_train[idx]
    y_train_subset = y_train[idx]

    return(X_train_subset, y_train_subset)

def prepareMILdata(path, wsi_filtered, use_smaller_bags=True, n_instances_max=100, balanced=False, subset=False, n_samples=None):
    '''
    Prepares training and testing data
    @param path: the path where the image embedding bags are stored
    @param wsi_filtered: dataframe containing file_names and sample mutation status
    @param use_smaller_bags: boolean, if True trims the bags to contain a maximum of n_patches_max instances (to speed up training)
    @param n_instances_max: integer, the maximum number of instances to be kept in each bag
    @param balanced: boolean, if True returns an equal number of samples per class in the training dataset
    @param subset: boolean, whether or not to select only a subset of the training data. If True, n_samples must be specified
    @param n_samples: integer, number of samples per class to be selected from training data (if subset is True)
    @return: X_train, X_test, y_train, y_test
    '''
    file_names = np.array(wsi_filtered['file_name'])
    labels = np.array(wsi_filtered['label'])
    X_train, X_test, y_train, y_test = train_test_split(file_names, labels, stratify=labels) # use file_names to make train/test split

    if balanced:
        X_train, y_train = balanceData(X_train, y_train)

    if subset:
        X_train, y_train = subsetData(X_train, y_train, n_samples=n_samples)

    train_bags=[]
    test_bags=[]
    for i, f in enumerate(X_train):
        bag = pd.read_csv(path + f + '.csv')
        bag = np.array(bag)
        # to speed up training, one can limit max number of instances per bag:
        if use_smaller_bags and len(bag) > n_instances_max:
            idx = np.random.choice(len(bag), size=n_instances_max, replace=False) # select random unique instances
            new_bag = bag[idx]
        else:
            new_bag = bag
        train_bags.append(new_bag)
        print('Preparing X_train... ' + str(i + 1) + '/' + str(len(X_train)))

    for i, f in enumerate(X_test):
        bag = pd.read_csv(path + f + '.csv')
        bag = np.array(bag)
        test_bags.append(bag)
        print('Preparing X_test... ' + str(i + 1) + '/' + str(len(X_test)))

    X_train = train_bags
    X_test = test_bags

    return(X_train, X_test, y_train, y_test)

def train_MIL(X_train, X_test, y_train, n_epochs=5, batch_size=5):
    '''
    Trains an Attention-based Deep Multiple Instance Learning model on image embeddings using stratified kfold cross-validation
    @param X_train: the training bags
    @param X_test: the test bags
    @param y_train: the train labels
    @param n_epochs: number of epochs for training
    @param batch_size: the batch size for training
    @return: model: the trained model
    @return history: the training history
    '''
    # get maximum number of instances per bag (for padding)
    max_len_train = np.max([len(bag) for bag in X_train])
    max_len_test = np.max([len(bag) for bag in X_test])
    max_ = np.max([max_len_train, max_len_test])

    trainer = Trainer()
    metrics = [AUC]
    model = AttentionDeepPoolingMil(gated=True, threshold=0.4)
    pipeline = [('padding', Padding(max_len=max_))]
    trainer.prepare(model, preprocess_pipeline=pipeline, metrics=metrics)
    valid = StratifiedKFold(n_splits=5, shuffle=True)

    history = trainer.fit(X_train, y_train,
                          model__epochs=n_epochs,
                          model__batch_size=batch_size,
                          sample_weights='balanced',
                          validation_strategy=valid,
                          verbose=1)

    return (trainer, history)

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


# STK11 mut vs wt classification ###################################################################################
# Load data ########################################################################################################
sample_info = pd.read_csv('data/data_for_ML/sample_info_tumor_only_no_outliers.csv') # dataframe containing list of patients to be useds
wsi_list = pd.read_csv('data/tcga-dataset/LUAD.csv') # list of wsi file_names

# Pre-processing ###################################################################################################
# Filter wsi_list to keep only samples from patients in patient_list and add mutation status to samples
wsi_filtered = filterWSI(wsi_list, sample_info)

# Prepare data for MIL
X_train, X_test, y_train, y_test = prepareMILdata(path='data/tcga-dataset/tcga_lung_data_feats/',
                                                  wsi_filtered=wsi_filtered,
                                                  use_smaller_bags=True,
                                                  n_instances_max=100,
                                                  balanced=True)
# Save data
pickleSave(object = X_train, file_name = 'data/data_for_ML/train_image_bags_100_balanced.pkl')
pickleSave(object = y_train, file_name = 'data/data_for_ML/train_image_labels_100_balanced.pkl')
pickleSave(object = X_test, file_name = 'data/data_for_ML/test_image_bags_100_balanced.pkl')
pickleSave(object = y_test, file_name = 'data/data_for_ML/test_image_labels_100_balanced.pkl')

# MIL training #####################################################################################################
X_train = pickleOpen('data/data_for_ML/train_image_bags_100.pkl')
y_train = pickleOpen('data/data_for_ML/train_image_labels_100.pkl')
X_test = pickleOpen('data/data_for_ML/test_image_bags_100.pkl')
y_test = pickleOpen('data/data_for_ML/test_image_labels_100.pkl')

trainer, history = train_MIL(X_train, X_test, y_train, n_epochs=5, batch_size=5)

# Evaluate model
y_pred = trainer.predict(X_test)
y_pred_train = trainer.predict(X_train)
fpr, tpr, thresh = roc_curve(y_train, y_pred_train, pos_label=1)
fig = plotROC(fpr, tpr, title='STK11 MUT vs WT\nTrain Set')
fig.savefig('results/deep_attention_MIL/STK11_100_ROC_train.png')
auc_score = roc_auc_score(y_train, y_pred_train)

# Check prediction probabilities
fig, ax = plt.subplots()
plt.hist(y_pred_train)
plt.title('Probability distribution for predictions\nTrain set', weight='bold')
fig.savefig('results/deep_attention_MIL/STK11_100_pred_prob_train.png')


# LUAD vs LUSC classification ######################################################################################
# Load data ########################################################################################################
LUAD_list = pd.read_csv('data/tcga-dataset/LUAD.csv')
LUSC_list = pd.read_csv('data/tcga-dataset/LUSC.csv')
sample_list = pd.concat([LUAD_list, LUSC_list])

# Pre-processing ##################################################################################################
# generate labels (LUAD = 0, LUSC = 1)
LUAD_labels = np.repeat(0, len(LUAD_list))
LUSC_labels = np.repeat(1, len(LUSC_list))
labels = np.append(LUAD_labels, LUSC_labels)

# generate dataframe containing file names and their corresponding labels
wsi_df = pd.DataFrame()
wsi_df['file_name'] = sample_list['0'].str[20:]
wsi_df['label'] = labels

# check if all files exist and keep only the names of the ones that do
idx = [os.path.isfile('data/tcga-dataset/tcga_lung_data_feats/' + i + '.csv') for i in wsi_df['file_name']]
wsi_df = wsi_df[idx]

# Prepare data for MIL
X_train, X_test, y_train, y_test = prepareMILdata(path='data/tcga-dataset/tcga_lung_data_feats/',
                                                  wsi_filtered=wsi_df,
                                                  use_smaller_bags=True,
                                                  n_instances_max=100,
                                                  subset=True,
                                                  n_samples=52)

# Save data
pickleSave(object = X_train, file_name = 'data/data_for_ML/train_image_bags_LUSC_LUAD_100_subset.pkl')
pickleSave(object = y_train, file_name = 'data/data_for_ML/train_image_labels_LUSC_LUAD_100_subset.pkl')
pickleSave(object = X_test, file_name = 'data/data_for_ML/test_image_bags_LUSC_LUAD_100_subset.pkl')
pickleSave(object = y_test, file_name = 'data/data_for_ML/test_image_labels_LUSC_LUAD_100_subset.pkl')

# MIL training #####################################################################################################
X_train = pickleOpen('data/data_for_ML/train_image_bags_LUSC_LUAD_100_subset.pkl')
y_train = pickleOpen('data/data_for_ML/train_image_labels_LUSC_LUAD_100_subset.pkl')
X_test = pickleOpen('data/data_for_ML/test_image_bags_LUSC_LUAD_100_subset.pkl')
y_test = pickleOpen('data/data_for_ML/test_image_labels_LUSC_LUAD_100_subset.pkl')

# Train model
trainer_2, history_2 = train_MIL(X_train, X_test, y_train, n_epochs=5, batch_size=5)

# Evaluate model
y_pred = trainer_2.predict(X_test)
y_pred_train = trainer_2.predict(X_train)
fpr_2, tpr_2, thresh_2 = roc_curve(y_train, y_pred_train, pos_label=1)
fig = plotROC(fpr_2, tpr_2, title='LUAD vs LUSC\nTrain set')
fig.savefig('results/deep_attention_MIL/LUAD_LUSC_100_subset_ROC_train.png')
auc_score = roc_auc_score(y_train, y_pred_train)

# Check prediction probabilities
fig, ax = plt.subplots()
plt.hist(y_pred_train)
plt.title('Probability distribution for predictions\nTrain set', weight='bold')
fig.savefig('results/deep_attention_MIL/LUAD_LUSC_100_subset_pred_prob_train.png')

# Control - label shuffling ########################################################################################
np.random.shuffle(y_train)
trainer_3, history_3 = train_MIL(X_train, X_test, y_train, n_epochs=5, batch_size=5)

y_pred = trainer_3.predict(X_test)
y_pred_train = trainer_3.predict(X_train)
fpr_3, tpr_3, thresh_3 = roc_curve(y_train, y_pred_train, pos_label=1)
fig = plotROC(fpr_3, tpr_3, title='LUAD vs LUSC\nTrain set')
fig.savefig('results/deep_attention_MIL/LUAD_LUSC_100_subset_shuffled_ROC_train.png')
auc_score = roc_auc_score(y_train, y_pred_train)

# Check prediction probabilities
fig, ax = plt.subplots()
plt.hist(y_pred)
plt.title('Probability distribution for predictions\nTest set', weight='bold')
fig.savefig('results/deep_attention_MIL/LUAD_LUSC_100_subset_shuffled_pred_prob_test.png')