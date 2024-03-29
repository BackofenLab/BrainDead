# BrainDead - RNA classification via accessible k-mers

## Setup

The recommended way for running BrainDead locally is via Conda.

### Installing miniconda

If you don't have conda on your system, please follow the installtion insturction for the minimal conda setup here: https://docs.conda.io/en/latest/miniconda.html

### Importing the dependencies

This command creates the conda environment for all the neccessary dependencies: `conda env create --file conda-environment.yml`

### Starting the environment

The imported environment can be activate with the command `conda activate BrainDead `.

Afterwards you should be able to call BrainDead's scripts (see below), e.g. running the [example calls](#sample-calls-on-test-data).


## Command helps

### Generate kmer features with `generate_kmer_features.py`

```
usage: generate_kmer_features.py [-h] --kmers KMERS --fasta FASTA [--threads THREADS]
                                 [--batchsize BATCHSIZE] [--report-counts]
                                 [--out-csv OUT_CSV] [--minE-subopt MINE_SUBOPT]
                                 [--minE-intarna MINE_INTARNA]
                                 [--feature-context FEATURE_CONTEXT]

Generate kmer-based count features based on sequence plus RNAsubopt and IntaRNA position-
wise energies.. Sample calls: "python generate_kmer_features.py --kmers "AGA,GC,GGG"
--fasta test.fa --out-csv counts.csv""python generate_kmer_features.py --kmers
"AGA,GC,GGG" --fasta test.fa --out-csv "stdout" --report-counts --minE-subopt -5 --minE-
intarna -2"

optional arguments:
  -h, --help            show this help message and exit
  --kmers KMERS         List of kmers as a comma separated string e.g. "AGG,GA,GG"
  --fasta FASTA         Sequences to extract features from as a FASTA file
  --threads THREADS     Number of threads used for processing (default: 1) (WARNING:
                        threads > 1 will impair stdout prints
  --batchsize BATCHSIZE
                        If the number of processed fasta sequences is greater than batch
                        size batch processing will be applied. This will lower memory
                        consumption (default: 10000)
  --report-counts       Whether to report counts as integer, default is binary
                        nohit(0)-hit(1)
  --out-csv OUT_CSV     CSV File name to write counts, pass "stdout" for stdout
  --minE-subopt MINE_SUBOPT
                        Minimum free energy of the position on RNAsubopt result
  --minE-intarna MINE_INTARNA
                        Minimum free energy of the position on IntaRNA result
  --feature-context FEATURE_CONTEXT
                        feature groups (contexts) are to be generated by case-insensitive
                        single letter a - any context (just k-mer occurrence) s -
                        unpaired in stable intra-molecular structures h - unpaired in
                        stable inter-molecular homo-duplex RRIs u - unpaired in both in
                        (s) and (h)
```

### Fit and predict ML model with `fit_predict.py`

```
usage: fit_predict.py [-h] --features-train FEATURES_TRAIN
                      [--features-train-index-col FEATURES_TRAIN_INDEX_COL]
                      [--features-train-header FEATURES_TRAIN_HEADER]
                      --labels-train LABELS_TRAIN
                      [--labels-train-index-col LABELS_TRAIN_INDEX_COL]
                      [--labels-train-header LABELS_TRAIN_HEADER]
                      [--model-choice {SVM-rbf,SVM-linear,Logistic-liblinear,Logistic-lbfgs}]
                      [--save-model] [--load-model LOAD_MODEL]
                      [--out-model OUT_MODEL] [--predict]
                      [--features-predict FEATURES_PREDICT]
                      [--features-predict-index-col FEATURES_PREDICT_INDEX_COL]
                      [--features-predict-header FEATURES_PREDICT_HEADER]
                      [--out-predict-labels OUT_PREDICT_LABELS]
                      [--not-validate-indexes] [--not-unique-indexes]
                      [--store-weights] [--standardize-scaling]

Train a machine learning model for the tabular CSV inputs, and predict on
optional CSV input Sample call: "python fit_predict.py --features-train
test_features.csv --labels-train test_labels.csv --features-train-header 0
--model SVM-rbf""python fit_predict.py --predict --features-predict
test_features.csv --features-predict-header 0 --load-model model.pkl""python
fit_predict.py --features-train test_features.csv --labels-train
test_labels.csv --features-train-header 0 --save-model --out-model jj.pkl
--predict --features-predict test_features.csv --features-predict-header 0
--load-model model.pkl"

optional arguments:
  -h, --help            show this help message and exit
  --features-train FEATURES_TRAIN
                        Input features to train model in Comma separated
                        format (CSV)
  --features-train-index-col FEATURES_TRAIN_INDEX_COL
                        The column number with index keys are provided
  --features-train-header FEATURES_TRAIN_HEADER
                        Whether features CSV has a header line
  --labels-train LABELS_TRAIN
                        Input reference labels to train model as one-column
                        Comma separated format (CSV)
  --labels-train-index-col LABELS_TRAIN_INDEX_COL
                        The column number with index keys are provided
  --labels-train-header LABELS_TRAIN_HEADER
                        Whether labels CSV has a header line
  --model-choice {SVM-rbf,SVM-linear,Logistic-liblinear,Logistic-lbfgs}
                        Included classifiers from scikit, see https://scikit-
                        learn.org/stable/modules/svm.html#classification
  --save-model          Save the trained scikit model into file
  --load-model LOAD_MODEL
                        Load model from file, instead of training on CSV input
  --out-model OUT_MODEL
                        Save the trained scikit model into file
  --predict             Optionally predict using the trained model, requires
                        --features-predict option
  --features-predict FEATURES_PREDICT
                        Input features to predict in comma separated format
                        (CSV)
  --features-predict-index-col FEATURES_PREDICT_INDEX_COL
                        The column number with index keys are provided
  --features-predict-header FEATURES_PREDICT_HEADER
                        Whether features CSV has a header line
  --out-predict-labels OUT_PREDICT_LABELS
                        File to store the predicted labels
  --not-validate-indexes
                        Check keys are consistently the same in lables and
                        features
  --not-unique-indexes  Check keys are not duplicated
  --store-weights       Save SVM coefficient weights
  --standardize-scaling
                        Do not scale and standardize features

```

## Sample calls on test data:

1. Generating k-mer sequence and structure features
```
python src/generate_kmer_features.py --kmers "AGA,GC,GGG" --fasta data/test.fa --out-csv data/counts.csv
```
2. Training and saving the ML model:
```
$ python src/fit_predict.py --features-train data/190411-mTLR7-features.csv  --labels-train data/190411-mTLR7-labels.csv --features-train-header 0  --save-model --out-model data/190411-mTLR7.pkl
```
3. Training and predicting:
```
$ python src/fit_predict.py --features-train data/190411-mTLR7-features.csv  --labels-train data/190411-mTLR7-labels.csv --features-train-header 0  --save-model --out-model data/190411-mTLR7.pkl --predict --features-predict data/190411-mTLR7-features.csv --features-predict-header 0 --out-predict-labels data/predict.cvs
```
4. Predicting from the saved model:
```
$ python src/fit_predict.py --predict --features-predict data/190411-mTLR7-features.csv --features-predict-header 0 --out-predict-labels predict.cvs --load-model data/190411-mTLR7.pkl 
```


