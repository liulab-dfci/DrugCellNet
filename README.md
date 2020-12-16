# DrugCellNet: Predicting Anticancer Drug Responses Using a Dual-Layer Integrated Cell Line-Drug Network Model

Author: Naiqian Zhang

Contact: naiqian@tongji.edu.cn

Apr 1, 2019

Dual Layer Network bridges cell-cell and drug-drug similarities to predict cancer cell drug responses, and derives its performance from the observation that drugs with similar molecular features elicit similar responses across cells, and cells with similar signatures respond in kind. Given datasets of cells, drugs, and drug-cell interactions, DLN can infer responses of unknown drugs and unknown cells based on similarity to examples in the training data. The algorithm, described in <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4587957/">Naiqian Zhang, et al.</a> optimizes three parameters while outperforming parameter-heavy machine learning approaches, and provides a computationally-efficient framework for computing pairwise drug-cell interactions.

There are three R files in the attachment. If you want to test our model, you need only run demo.R. The result of our model is the predicted value of the drug response. Because we normalized drug response value in the data precessing step, so the predicted value is a relative one. If you want to get the absolute predicted value of a given drug, the relative value should be multiply by standard deviation plus the mean of the raw drug response values.

# Tutorial

## download training and test data
```R
trainData <- "https://zenodo.org/record/2621169/files/Train.Rdata.zip"
temp <- tempfile()
download.file(trainData,temp)
load(unz(temp, "Train.Rdata"))


testData <- "https://zenodo.org/record/2621169/files/Test.Rdata"
temp <- tempfile()
download.file(testData,temp)
load(temp)

source("./duallayer_network_training.R")
source("./duallayer_network_predict.R")
```

## Train model from training set

input: 
 - Training.DrugResponse is the drug sensitivity matrix(column is drug,row is cellline).
 - Training.Expression is cell lines expression matrix (col is cell line,row is gene).
 - raining.DrugChemical is drug chemical feature matrix(col is feature, row is drug).
 - All the three data constitutes the training data set

output: the parameters of our dual-layer network model

```R
model.train = get.train.main(Training.DrugResponse,Training.Expression,Training.DrugChemical) 
```

## Predict

We use the a dual-layer integrated cell line-drug network model to predict the drug response of a given cell line and drug

input: 
 - Object.Expression is one column matrix represent the object cell line's gene expression (col is cel line,row is gene)
 - Object.DrugChemical is a matrix which only has one row represents the object drug's chemical structure

output: the predicted drug sensitivity of our dual-layer network model

Senerio I: Both cell line and drug are available in the training set, so we use dual-layer model to predict drug response
```R
# normalize test data
data.processed1 = get.processed.data(Training.DrugResponse,Training.Expression,Training.DrugChemical,Object1.Expression,Object1.DrugChemical)
# predict test sample using trained model 
result1 = get_Predict(data.processed1,model.train)
```

Senerio II: Drug is available in the training set but the cell line is new. We can use the upper layer model to predict
```R
# normalize test data
data.processed2 = get.processed.data(Training.DrugResponse,Training.Expression,Training.DrugChemical,Object2.Expression,Object2.DrugChemical)
# predict test sample using trained model 
result2 = get_Predict(data.processed2,model.train)
```

Senerio III: Cell line is in the training set but the drug is new. We use the lower layer model to predict
```R
# normalize test data
data.processed3 = get.processed.data(Training.DrugResponse,Training.Expression,Training.DrugChemical,Object3.Expression,Object3.DrugChemical)
# predict test sample using trained model 
result3 = get_Predict(data.processed3,model.train)
```
