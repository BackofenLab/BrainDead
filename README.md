# BrainDead - RNA classification via accessible k-mers

## Setup

To install and run BrainDead, you should use a commandline interface like `bash`.

### (1) Downloading BrainDead script files from Github

First, go to a folder where you want to "install" BrainDead, i.e. downloading the script files listed in this github repository. 

```sh
# download all files from 
wget https://github.com/BackofenLab/BrainDead/archive/refs/heads/main.zip
# uncompress the archive (might need an additional installation of "unzip")
unzip main.zip
# enter the unpacked and created folder
cd BrainDead-main
```

Afterwards, you find in

- `src` the script files (see below for usage)
- `data` the data files used for training the models etc.

The recommended way for running BrainDead locally is via Conda.

### (2) Installing miniconda (if not already done)

If you don't have conda on your system (check if `conda --version` gives any output or an error message), please follow the installtion instruction for the minimal conda setup available via the following link.

- [miniconda installation](https://docs.conda.io/en/latest/miniconda.html)

`conda` is a commandline tool to simplify the installation and local management of tools and their dependencies in linux-based systems.

### (3) Importing the dependencies via conda

Next, we want to install all files and tools needed to run BrainDead via conda.
To this end, we will create a "conda environment", i.e. kind of a folder that contains all tools not part of this github repository.

The needed files and tools are listed in the `conda-environment.yml` file within the `BrainDead-main` folder extracted in step (1).

The following command creates the environment and imports BrainDead's dependencies.

```sh
conda env create --file conda-environment.yml
```

### (4) Starting the conda environment

The created environment can be activate with the following command 

```sh
conda activate BrainDead
```

Afterwards you should be able to call BrainDead's scripts (see below).

*Note: the script files are located in the `src` subdirectory!`

At the end of the page, you find the [example calls](#sample-calls-on-test-data).


## Command helps

### Generate kmer features with `generate_kmer_features.py`

```
usage: generate_kmer_features.py [-h] --kmers KMERS --fasta FASTA [--report-counts]
                                 [--out-csv OUT_CSV] [--minE-subopt MINE_SUBOPT]
                                 [--minE-intarna MINE_INTARNA]

Generate kmer-based count features based on sequence plus RNAsubopt and IntaRNA position-
wise energies.. 
Sample calls: 
"python generate_kmer_features.py --kmers "AGA,GC,GGG" --fasta test.fa --out-csv counts.csv" 
"python generate_kmer_features.py --kmers "AGA,GC,GGG" --fasta test.fa --out-csv "stdout" --report-counts --minE-subopt -5 --minE-intarna -2"

optional arguments:
  -h, --help            show this help message and exit
  --kmers KMERS         List of kmers as a comma separated string e.g. "AGG,GA,GG"
  --fasta FASTA         Sequences to extract features from as a FASTA file
  --report-counts       Whether to report counts as integer, default is binary
                        nohit(0)-hit(1)
  --out-csv OUT_CSV     CSV File name to write counts, pass "stdout" for stdout
  --minE-subopt MINE_SUBOPT
                        Minimum free energy of the position on RNAsubopt result
  --minE-intarna MINE_INTARNA
                        Minimum free energy of the position on IntaRNA result
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
2. Training and saving the ML model (here an `SVM-rbf` model):
```
$ python src/fit_predict.py --model-choice 'SVM-rbf' --features-train data/190411-mTLR7-features.csv  --labels-train data/190411-mTLR7-labels.csv --features-train-header 0  --save-model --out-model data/190411-mTLR7.pkl
```
3. Predicting from the saved model:
```
$ python src/fit_predict.py --predict --features-predict data/190411-mTLR7-features.csv --features-predict-header 0 --out-predict-labels predict.cvs --load-model data/190411-mTLR7.pkl 
```
Alternatively:

2+3. Training and predicting in one step:
```
$ python src/fit_predict.py --features-train data/190411-mTLR7-features.csv  --labels-train data/190411-mTLR7-labels.csv --features-train-header 0  --save-model --out-model data/190411-mTLR7.pkl --predict --features-predict data/190411-mTLR7-features.csv --features-predict-header 0 --out-predict-labels data/predict.cvs
```


