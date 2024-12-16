#!/usr/bin/env python3

import pandas as pd
import argparse
import pickle
import os.path
import sys
import numpy
from sklearn.model_selection import cross_validate, cross_val_score
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_fscore_support
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

import warnings
warnings.filterwarnings('ignore') # to get rid of UndefinedMetricWarning

def fit_predict_on_top_features(df_rebate, ref_feature, df_all_features, num_top_features, scikit_model):
    
    top_features = df_rebate.iloc[:,0].values[:num_top_features].tolist()
    print ("Selected top feature according to reference feature importance {} are:\n {}".format(ref_feature, top_features))
    
    print('shape all features:', df_all_features.shape)

    
    df_selected_features = df_all_features[['microRNA.candidate_x',ref_feature]+top_features].set_index('microRNA.candidate_x').copy()
    print('orig shape:', df_selected_features.shape)

    # reference feature must not be NULL
    df_selected_features = df_selected_features[df_selected_features[ref_feature].notnull()]
    print('shape after ref-features null removal:',df_selected_features.shape)

    # Other feature are expected to not be NULL!?

    df_selected_features.dropna(inplace=True)
    print('shape after all-features null removal:',df_selected_features.shape)

    # Duplicated entries (due to different experient sources?) are removed
    df_selected_features.drop_duplicates(inplace=True)
    print('shape after row deduplication:', df_selected_features.shape)
    df_selected_features
    
    
    Y = df_selected_features[ref_feature].values
    X = df_selected_features[df_selected_features.columns.difference([ref_feature])]
    scikit_model.fit(X, Y)
    print("Training score:{%.3f}".format(scikit_model.score(X, Y)))
    
    cv_results = cross_val_score(scikit_model, X, Y, cv=3,
                               scoring='f1')
    print('3-fold CV F1', cv_results)

    #print (''.join([ '{}, {}{}'.format(s[0],s[1],'\n') for s in zip(Y,scikit_model.predict(X))] ))

def fit_dfs_trains(df_features, df_labels, scikit_model, top_features=None, validate_indexes=False, unique_indexes=False):

    if validate_indexes and (not df_labels.index.equals(df_features.index)):
        raise RuntimeError("Error: label and feature keys don't match\n\t{}\n\t{}".format(df_labels.index,df_features.index))
    if unique_indexes is True:
        if not (df_labels.index.is_unique and df_features.index.is_unique):
            raise RuntimeError ("Error: expected unique_indexes but found: \n\t{}\n\t{}".format(df_labels.index,df_features.index))
    if len(df_labels.columns) != 1:
        raise RuntimeError("Expected one column for labels csv file, found: {}".format(len(df_labels.columns)))
    if top_features is not None: 
        df_selected_features = df_features[top_features].copy()
    else:
        df_selected_features = df_features.copy()

    
    print('shape train features:', df_selected_features.shape)

    # reference feature must not be NULL
    if df_selected_features.isnull().values.any():
        df_selected_features.dropna(inplace=True)
        print('Warning: features with NA values detected and discarded\n\tshape after null removal:',
            df_selected_features.shape)

    #if len(df_selected_features[df_selected_features.duplicated(keep=False)]) != 0:
    
    #    # Duplicated entries (due to different experient sources?) are removed
    #    df_selected_features.drop_duplicates(inplace=True)
    #    print('Warning: duplicated rows are detected and discarded\n\tshape after row deduplication:', df_selected_features.shape)
    
    

 

    y = df_labels.iloc[:,0].values
    X = df_selected_features.values
    #print("X", X)
    #print("Y", Y)
    n_splits = 2
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=10, random_state=2222)
    #skf = StratifiedKFold(n_splits=3, random_state=None, shuffle=True)
    arr_PRFS = list()#numpy.array()
    for train_index, test_index in rskf.split(X, y):
        # print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]    
        y_test_pred = scikit_model.fit(X_train, y_train).predict(X_test)
        PRFS = precision_recall_fscore_support(y_test, y_test_pred, average='binary')[:3]
        arr_PRFS.append(PRFS)

    arr_PRFS = numpy.array(arr_PRFS)

    df_metrics = pd.DataFrame(numpy.vstack((numpy.mean(arr_PRFS, axis=0), numpy.std(arr_PRFS, axis=0))),
    columns=["Precision","Recall","F1"], index=['metric','std']).round(3)

    print("10-rep {}-fold CV:".format(n_splits))
    print(df_metrics)
    
    scikit_model.fit(X, y)
    print("Training score: {0:.3f}".format(scikit_model.score(X, y)))


def is_valid_file(file_name):
    if os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        raise FileNotFoundError(os.path.abspath(file_name))
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Train a machine learning model for the tabular CSV inputs, and predict on optional CSV input'\
        '\nSample call: \"python fit_predict.py --features-train test_features.csv  --labels-train test_labels.csv --features-train-header 0 --model SVM-rbf\"'\
        '"python fit_predict.py --predict --features-predict test_features.csv --features-predict-header 0 --load-model model.pkl"'\
        '"python fit_predict.py --features-train test_features.csv  --labels-train test_labels.csv --features-train-header 0  --save-model --out-model jj.pkl --predict --features-predict test_features.csv --features-predict-header 0 --load-model model.pkl"')
    parser.add_argument('--features-train', required='--load-model' not in sys.argv, type=is_valid_file, help='Input features to train model in Comma separated format (CSV)')
    parser.add_argument('--features-train-index-col', default=0, help='The column number with index keys are provided')
    parser.add_argument('--features-train-header', default=None, help='Whether features CSV has a header line ')
    parser.add_argument('--labels-train', required='--load-model' not in sys.argv, type=is_valid_file, help='Input reference labels to train model as one-column Comma separated format (CSV)')
    parser.add_argument('--labels-train-index-col', default=0, help='The column number with index keys are provided')
    parser.add_argument('--labels-train-header', default=None, type=int,help='Whether labels CSV has a header line')
    parser.add_argument('--model-choice', choices=['SVM-rbf','SVM-linear','Logistic-liblinear','Logistic-lbfgs'],
        default='Logistic-liblinear', help='Included classifiers from scikit, see https://scikit-learn.org/stable/modules/svm.html#classification')
    parser.add_argument('--save-model', action='store_true', help='Save the trained scikit model into file'),
    parser.add_argument('--load-model',required=False, type=is_valid_file, help='Load model from file, instead of training on CSV input'),    
    parser.add_argument('--out-model', type=str, default='model.pkl', help='Save the trained scikit model into file'),
    parser.add_argument('--predict', action='store_true', help='Optionally predict using the trained model, requires --features-predict option'),
    parser.add_argument('--features-predict', required=False, type=is_valid_file, help='Input features to predict in comma separated format (CSV)')
    parser.add_argument('--features-predict-index-col', default=0, help='The column number with index keys are provided')
    parser.add_argument('--features-predict-header', default=None, help='Whether features CSV has a header line ')
    parser.add_argument('--out-predict-labels', type=str, default='predicted.csv', help='File to store the predicted labels'),
    parser.add_argument('--not-validate-indexes', action="store_true", help='Check keys are consistently the same in lables and features')
    parser.add_argument('--not-unique-indexes', action="store_true", help='Check keys are not duplicated')
    parser.add_argument('--store-weights', action="store_true", help='Save SVM coefficient weights')
    parser.add_argument('--standardize-scaling', action="store_true", help='Do not scale and standardize features')


# Save to file in the current working directory

    args = parser.parse_args()
    #print(args)
    
    if args.predict is True and not args.features_predict:
        raise RuntimeError ("--predict option requires assigning the csv file --features-predict")

    if args.load_model and (args.features_train or args.labels_train):
        print("Warning: using the loaded model, train CSV inputs are ignored")
    
    
    if args.load_model:
        # Load from file
        with open(args.load_model, 'rb') as inmodel:
            model = pickle.load(inmodel)
    else:
        if args.model_choice == 'SVM-rbf':
            model = svm.SVC(kernel='rbf', gamma='scale', probability=True, class_weight="balanced")
        if args.model_choice == 'SVM-linear':
            model = svm.SVC(kernel='linear', gamma='scale', probability=True, class_weight="balanced")
        elif args.model_choice == 'Logistic-liblinear':
            model = LogisticRegression(solver='liblinear', class_weight="balanced")    
        elif args.model_choice == 'Logistic-lbfgs':
            model = LogisticRegression(solver='lbfgs', class_weight="balanced")    

        if args.standardize_scaling:
            # build pipe: first standardize by substracting mean and dividing std
            # next do classificaiton
            model = make_pipeline(StandardScaler(), model) # model is now a actually a pipeline that first scales the data

        print(args.features_train, "index_col=",args.features_train_index_col, "header=", args.features_train_header)
        fit_dfs_trains(
            pd.read_csv(args.features_train, index_col=args.features_train_index_col, header=int(args.features_train_header) if args.features_train_header else None), 
            pd.read_csv(args.labels_train, index_col=args.labels_train_index_col, header=int(args.labels_train_header) if args.labels_train_header else None), 
            scikit_model =  model, validate_indexes= not args.not_validate_indexes, unique_indexes= not args.not_unique_indexes)
    if args.store_weights is True:
            
            print('\n'.join([str(model.coef_[0][i:i+4]) for i in range(0,len(model.coef_[0]),4)]))

    if args.predict is True:
        df_features_predict = pd.read_csv(args.features_predict, index_col=args.features_predict_index_col, header=int(args.features_predict_header) if args.features_predict_header else None)
        X_predict = df_features_predict.values
        y_predict = model.predict(X_predict)
        class_probability_predict = model.predict_proba(X_predict)
        stacked_arr = numpy.column_stack((y_predict, class_probability_predict))

        print('Predicted for {}'.format(args.features_predict))
        # with open(args.out_predict_labels, 'w') as outlabels:
        #     outlabels.write('predicted_class,prob_{},prob_{}\n'.format(model.classes_[0],model.classes_[1]))
        #     y_predict.tofile(outlabels,sep="\n")
        df_out_predict = pd.DataFrame(stacked_arr, index=df_features_predict.index)
        df_out_predict.sort_values(by=[df_out_predict.columns[1],df_out_predict.columns[2]], ascending=False, inplace=True)
        df_out_predict['rank'] = df_out_predict.iloc[:,2].rank(method="first",ascending=False)
        df_out_predict.to_csv(args.out_predict_labels, index=True,
        header=['predicted_class','prob_{}'.format(model.classes_[0]),'prob_{}'.format(model.classes_[1]),'rank'],float_format='%.3F')

    if args.save_model is True:
        with open(args.out_model, 'wb') as outf:
            pickle.dump(model, outf)
            print ("Model saved to: {}".format(args.out_model))


