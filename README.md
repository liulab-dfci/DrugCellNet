# DrugCellNet
Predicting Anticancer Drug Responses Using a Dual-Layer Integrated Cell Line-Drug Network Model

- There are three R files in the attachment. If you want to test our model, you need only run demo.R.
- The result of our model is the predicted value of the drug response. Because we normalized drug response value in the data precessing step, so the predicted value is a relative one. If you want to get the absolute predicted value of a given drug, the relative value should be multiply by standard deviation plus the mean of the raw drug response values.
