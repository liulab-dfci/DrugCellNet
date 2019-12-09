
####################################################
## Titile: Predicting Anticancer Drug Responses Using a Dual-Layer Integrated
## Cell Line-Drug Network Model
## Author: Naiqian Zhang
## Contact: naiqian@tongji.edu.cn
## Apr 1, 2019
####################################################

# We propose a dual-layer integrated cell line-drug network model, which uses both cell line similarity network (CSN) data and drug similarity network (DSN) data to predict the drug response of a given cell line using a weighted model. 

# input: Training.DrugResponse is the drug sensitivity matrix(column is drug,row is cellline).
#        Training.Expression is cell lines expression matrix (col is cell line,row is gene).
#        Training.DrugChemical is drug chemical feature matrix(col is feature, row is drug).
#        All the three data constitutes the training data set
# output: the parameters of our dual-layer network model

####################################################

# process and normalize the training data 
#input:training data set
#output: normalized drug response, corrlation of cells based on gene expression,correlation of drugs based on drug chemical characters
get.Training.processed.data=function(Training.DrugResponse,Training.Expression,Training.DrugChemical){
  
  DrugResponse.nor=scale(Training.DrugResponse)
  CellExpression.nor=t(scale(t(as.matrix(log(Training.Expression+1)))))
  DrugChemical.nor=t(scale(Training.DrugChemical))
  DrugChemical.nor[1:4,1:4]
  dim(DrugChemical.nor)
  Drugs=intersect(colnames(DrugResponse.nor),colnames(DrugChemical.nor))
  Cells=intersect(rownames(DrugResponse.nor),colnames(CellExpression.nor))
  
  TR.DrugResponse=DrugResponse.nor[Cells,Drugs]
  CellExpression=CellExpression.nor[,Cells]
  DrugChemical=DrugChemical.nor[,Drugs]

  Cell_correlation=cor(CellExpression,use = "pairwise.complete.obs")
  Drug_correlation=cor(DrugChemical,use = "pairwise.complete.obs")
  
  return(list(TR.DrugResponse,Cell_correlation,Drug_correlation))
}


# train sigma which is the parameter of upper layer model named Cellline similarity network model

get.OneDrug.CNS.Predict=function(drugname,TR.DrugResponse,sigma,Cell_correlation){
  true=TR.DrugResponse[,drugname,drop=F]
  predict=c()
  cells=rownames(true)
  for(j in 1:nrow(true)){
    cell=cells[j]
    user_act=as.vector(true[-which(rownames(true)==cell),])
    user_cor=Cell_correlation[-which(rownames(Cell_correlation)==cell),cell]
    user_wight=exp(-((1-user_cor)^2/(2*(sigma^2))))
    pre=sum(user_act*user_wight,na.rm = T)/sum(user_wight,na.rm = T)
    predict=append(predict,pre)
  }
  names(predict)=rownames(true)
  return(predict)
}


get.AllDrug.CNS.Predict=function(TR.DrugResponse,sigma,Cell_correlation){
  Drugs=colnames(TR.DrugResponse)
  AllDrug.CNS.Predict=c()
  for(i in 1:length(Drugs)){
    drugname=Drugs[i]
    
    OneDrug.CNS.Predict=get.OneDrug.CNS.Predict(drugname,TR.DrugResponse,sigma,Cell_correlation)
    AllDrug.CNS.Predict=cbind(AllDrug.CNS.Predict,OneDrug.CNS.Predict)
  }
  colnames(AllDrug.CNS.Predict)=Drugs
  rownames(AllDrug.CNS.Predict)=rownames(TR.DrugResponse)
  return(AllDrug.CNS.Predict)
}



get_sigma=function(TR.DrugResponse,Cell_correlation){
  ERROR=c()
  for(sigma in seq(0.1,1,0.01)){
    AllDrug.CNS.Predict=get.AllDrug.CNS.Predict(TR.DrugResponse,sigma,Cell_correlation)
    Error=sum((AllDrug.CNS.Predict-TR.DrugResponse)^2,na.rm = T)
    ERROR=append(ERROR,Error)
  }
  names(ERROR)=seq(0.1,1,0.01)
  sigma=as.numeric(names(which.min(ERROR)))
  return(sigma)
}

# train delt which is the parameter of lower layer model named drug similarity network model

get.OneCell.DNS.Predict=function(cellname,TR.DrugResponse,delt,Drug_correlation){
  true=TR.DrugResponse[cellname,,drop=F]
  predict=c()
  drugs=colnames(true)
  for(i in 1:ncol(true)){
    drug=drugs[i]
    user_act=as.vector(true[,-which(colnames(true)==drug)])
    user_cor=Drug_correlation[-which(rownames(Drug_correlation)==drug),drug]
    user_wight=exp(-((1-user_cor)^2/(2*(delt^2))))
    pre=sum(user_act*user_wight,na.rm = T)/sum(user_wight,na.rm = T)
    predict=append(predict,pre)
  }
  names(predict)=rownames(true)
  return(predict)
}


get.AllCell.DNS.Predict=function(TR.DrugResponse,delt,Drug_correlation){
  Cells=rownames(TR.DrugResponse)
  Drugs=colnames(TR.DrugResponse)
  AllCell.DNS.Predict=c()
  for(i in 1:length(Cells)){
    cellname=Cells[i]
    
    OneCell.DNS.Predict=get.OneCell.DNS.Predict(cellname,TR.DrugResponse,delt,Drug_correlation)
    AllCell.DNS.Predict=rbind(AllCell.DNS.Predict,OneCell.DNS.Predict)
  }
  colnames(AllCell.DNS.Predict)=Drugs
  rownames(AllCell.DNS.Predict)=rownames(TR.DrugResponse)
  return(AllCell.DNS.Predict)
}


get_delt=function(TR.DrugResponse,Drug_correlation){
  ERROR=c()
  for(delt in seq(0.1,1,0.01)){
    AllCell.DNS.Predict=get.AllCell.DNS.Predict(TR.DrugResponse,delt,Drug_correlation)
    Error=sum((AllCell.DNS.Predict-TR.DrugResponse)^2,na.rm = T)
    ERROR=append(ERROR,Error)
  }
  names(ERROR)=seq(0.1,1,0.01)
  delt=as.numeric(names(which.min(ERROR)))
  return(delt)
}


# train lammda which is the parameter of dual network model, we use lammda joint the upper and lower layer network models

get_lammda=function(TR.DrugResponse,CNS.Predict,DNS.Predict){
  drugs=colnames(TR.DrugResponse)
  Lambda=c()
  for(i in 1:length(drugs)){
    drug=drugs[i]
    CNS.Predict.drugi=CNS.Predict[,drug]
    DNS.Predict.drugi=DNS.Predict[,drug]
    Ture.drugi=TR.DrugResponse[,drug]
    Error=c()
    for(lambda in seq(0,1,0.01)){
      DualN.predict.drugi=(1-lambda)*CNS.Predict.drugi+lambda*DNS.Predict.drugi
      error=sum((Ture.drugi-DualN.predict.drugi)^2,na.rm=T)
      Error=append(Error,error)
      }
    names(Error)=seq(0,1,0.01)
    lambda_drugi=as.numeric(names(which.min(Error)))
    Lambda=append(Lambda,lambda_drugi)
  }
  names(Lambda)=drugs
  return(Lambda)
}




#train all the parameters in our dual-layer network model

Training.model=function(TR.DrugResponse,Cell_correlation,Drug_correlation)
{
  parameter.sigma=get_sigma(TR.DrugResponse,Cell_correlation)
  parameter.delt=get_delt(TR.DrugResponse,Drug_correlation)
  
  CNS.Predict=get.AllDrug.CNS.Predict(TR.DrugResponse,parameter.sigma,Cell_correlation)
  DNS.Predict=get.AllCell.DNS.Predict(TR.DrugResponse,parameter.delt,Drug_correlation)
  
  parameter.lammda=get_lammda(TR.DrugResponse,CNS.Predict,DNS.Predict)
  
  #Dual.Predict=get_Dual.Predict(TR.DrugResponse,CNS.Predict,DNS.Predict,parameter.lammda)
  parameter=list(parameter.sigma,parameter.delt,parameter.lammda)
  names(parameter)=c("sigma","delt","lammda")
  return(parameter)
  
}

# train the model in training data set, get the parameters in dual-layer network model

get.train.main=function(Training.DrugResponse,Training.Expression,Training.DrugChemical){
  Training.Processed.data=get.Training.processed.data(Training.DrugResponse,Training.Expression,Training.DrugChemical)
  TR.DrugResponse=Training.Processed.data[[1]]
  Cell_correlation=Training.Processed.data[[2]]
  Drug_correlation=Training.Processed.data[[3]]
  Training.parameter=Training.model(TR.DrugResponse,Cell_correlation,Drug_correlation)
  return(Training.parameter)
}

