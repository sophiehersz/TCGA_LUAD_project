import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_validate
import matplotlib.pyplot as plt
import seaborn as sns

# Functions #########################################################################################################
def prepareData(image_embeddings, mutation_status):
    '''
    @param image_embeddings: dataframe with single image embeddings per patient
    @oaram mutation_status: dataframe containing STK11 mutation status
                            sample order must be the same as for image_embeddings
    '''
    y = mutation_status.values.ravel()
    X = image_embeddings.values
    return(X, y)

def nestedCV_images(model, params, embeddings_list, mutation_status):
    '''
    Performs nested cross-validation to optimise model hyperparameters and evaluate performance
    @param model: the model to be trained
    @param params: hyperparameter search space
    @param: image_embeddings: dataframe with single image embeddings per patient
    @oaram mutation_status: dataframe containing STK11 mutation status
                            sample order must be the same as for image_embeddings
    Returns cv_results: full cross-validation results
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
    cv_results = dict.fromkeys(embeddings_list.keys(), None)
    df_scores = pd.DataFrame(columns=embeddings_list.keys())

    for i, k in enumerate(embeddings_list.keys()):
        # prepare X and y
        X, y = prepareData(embeddings_list[k], mutation_status)

        # execute the nested cross-validation and get results
        cv_out = cross_validate(search, X, y,
                                scoring='roc_auc',
                                cv=cv_outer,
                                n_jobs=-1,
                                return_train_score=True,
                                return_estimator=True)

        cv_results[k] = cv_out
        df_scores[k] = cv_out['test_score']
        print('embeddings ' + str(i + 1) + ' of ' + str(len(embeddings_list.keys())))

    # compute average scores and std
    average_scores = pd.DataFrame(index=embeddings_list.keys(), columns=['average_roc_auc', 'std'])
    average_scores['average_roc_auc'] = df_scores.apply(np.mean, axis='index')
    average_scores['std'] = df_scores.apply(np.std, axis='index')

    return (cv_results, df_scores, average_scores)

def plot_nested_cv_results(scores, y_lim=None, title=None, rotation=0):
    fig, ax = plt.subplots()
    sns.boxplot(data=scores, palette='Greys')
    sns.despine()
    plt.ylim(y_lim)
    plt.xticks(rotation=rotation, fontsize=13)
    plt.yticks(fontsize=13)
    #plt.xlabel('Gene sets', fontsize=15, weight='bold')
    plt.ylabel('ROC-AUC score', fontsize=15, weight='bold')
    plt.title(title, fontsize=17, weight='bold')
    return(fig, ax)

# Load data ##########################################################################################################
mean_embeddings = pd.read_csv('data/data_for_ML/mean_image_embeddings.csv')
median_embeddings = pd.read_csv('data/data_for_ML/median_image_embeddings.csv')
mutation_status = pd.read_csv('data/data_for_ML/mutation_status_embeddings.csv')

# Random Forest ##############################################################################################
# define the model
RF_model = RandomForestClassifier(class_weight='balanced',
                                  random_state=0)

# define hyperparameter search space
RF_params = {'n_estimators': [50, 100],
             'max_depth': [2, 3, 4]}
             #'min_samples_split': [4, 6, 8]}

embeddings_list = {'Mean aggregation': mean_embeddings,
                   'Median aggregation': median_embeddings}

# run nested CV
RF_cv_results, RF_scores, RF_average_scores = nestedCV_images(RF_model, RF_params, embeddings_list, mutation_status)
fig, ax = plot_nested_cv_results(RF_scores, y_lim=[0,1], title='Random Forest\nnested cross-validation (k=10)')
fig.savefig('results/random_forest/RF_nested_cv_image_embeddings.png', bbox_inches='tight')