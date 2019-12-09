
rm(list = ls())

## download training and test data
trainData <- "https://zenodo.org/record/2621169/files/Train.Rdata.zip"
temp <- tempfile()
download.file(trainData,temp)
load(unz(temp, "Train.Rdata"))

testData <- "https://zenodo.org/record/2621169/files/Test.Rdata"
temp <- tempfile()
download.file(testData,temp)
load(temp)

# load("Test.Rdata")
# load("Train")

source("./duallayer_network_training.R")
source("./duallayer_network_predict.R")

## Train model from training set
model.train = get.train.main(Training.DrugResponse,Training.Expression,Training.DrugChemical) 

## Senerio I: Both cell line and drug are available in the training set, so we use dual-layer model to predict drug response
## normalize test data
data.processed1 = get.processed.data(Training.DrugResponse,Training.Expression,Training.DrugChemical,Object1.Expression,Object1.DrugChemical)
## predict test sample using trained model 
result1 = get_Predict(data.processed1,model.train)


## Senerio II: Drug is available in the training set but the cell line is new. We can use the upper layer model to predict
## normalize test data
data.processed2 = get.processed.data(Training.DrugResponse,Training.Expression,Training.DrugChemical,Object2.Expression,Object2.DrugChemical)
## predict test sample using trained model 
result2 = get_Predict(data.processed2,model.train)


## Senerio III: Cell line is in the training set but the drug is new. We use the lower layer model to predict
## normalize test data
data.processed3 = get.processed.data(Training.DrugResponse,Training.Expression,Training.DrugChemical,Object3.Expression,Object3.DrugChemical)
## predict test sample using trained model 
result3 = get_Predict(data.processed3,model.train)



