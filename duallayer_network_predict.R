####################################################
## Titile: Predicting Anticancer Drug Responses Using a Dual-Layer Integrated
## Cell Line-Drug Network Model
## Author: Naiqian Zhang
## Contact: naiqian@tongji.edu.cn
## Apr 1, 2019
####################################################

# We use the a dual-layer integrated cell line-drug network model to predict the drug response of a given cell line and drug
# input: Training.DrugResponse is the drug sensitivity matrix(column is drug,row is cellline).
#        Training.Expression is cell lines expression matrix (col is cell line,row is gene).
#        Training.DrugChemical is drug chemical feature matrix(col is feature, row is drug).
#        The three data sets constitute the training data
#        Object.Expression is one column matrix represent the object cell line's gene expression (col is cel line,row is gene)
#        Object.DrugChemical is a matrix which only has one row represents the object drug's chemical structure
#
#
# output: the predicted drug sensitivity of our dual-layer network model


# process and normalize the training data and test data 

get.processed.data=function(Training.DrugResponse,Training.Expression,Training.DrugChemical,Object.Expression,Object.DrugChemical){

  genes=intersect(rownames(Training.Expression),rownames(Object.Expression))
  chemicalFeatures=intersect(colnames(Training.DrugChemical),colnames(Object.DrugChemical))
  
  celllines=intersect(colnames(Training.Expression),rownames(Training.DrugResponse))
  drugs=intersect(colnames(Training.DrugResponse),rownames(Training.DrugChemical))
  
  TR.DrugResponse=scale(Training.DrugResponse)[celllines,drugs]
  

  CellExpression=cbind(Training.Expression[genes,],Object.Expression[genes,,drop=F])
  normlized_CellExpression=as.data.frame(t(scale(t(as.matrix(log(CellExpression+1))))))
  
  DrugChemical=rbind(Training.DrugChemical[,chemicalFeatures],Object.DrugChemical[,chemicalFeatures,drop=F])
  normlized_DrugChemical=as.data.frame(t(scale(DrugChemical)))
    

  TR.CellExpression=as.data.frame(normlized_CellExpression[,celllines])
  TR.DrugChemical=as.data.frame(normlized_DrugChemical[,drugs])
  
  OB.CellExpression=normlized_CellExpression[,colnames(Object.Expression),drop=F]
  OB.DrugChemical=normlized_DrugChemical[,rownames(Object.DrugChemical),drop=F]

  
  return(list(TR.DrugResponse,TR.CellExpression,TR.DrugChemical,OB.CellExpression,OB.DrugChemical))
}


# predict drug sensitivity using upper layer model named Cellline similarity network model
get_upper_predict=function(TR.CellExpression,OB.CellExpression,sigma,DrugName,TR.DrugResponse){
  correlation.cells=cor(OB.CellExpression,TR.CellExpression,use = "pairwise.complete.obs")
  cells_com=setdiff(intersect(colnames(TR.CellExpression),rownames(TR.DrugResponse)),colnames(OB.CellExpression))
  correlation.cells_com=correlation.cells[,cells_com]
  

  DrugResponse_com=TR.DrugResponse[cells_com,DrugName]
  wight=exp(-((1-correlation.cells_com)^2/(2*(sigma^2))))
  predict=sum(DrugResponse_com*wight,na.rm = T)/sum(wight,na.rm = T)
  return(predict)
  }


# predict drug sensitivity using lower layer model named Drug similarity network model
get_lower_predict=function(TR.DrugChemical,OB.DrugChemical,delt,CellName,TR.DrugResponse){
  correlation.drugs=cor(OB.DrugChemical,TR.DrugChemical,use = "pairwise.complete.obs")
  drugs_com=setdiff(intersect(colnames(TR.DrugChemical),colnames(TR.DrugResponse)),colnames(OB.DrugChemical))
  correlation.drugs_com=correlation.drugs[,drugs_com]
  
  DrugResponse_com=TR.DrugResponse[CellName,drugs_com]
  wight=exp(-((1-correlation.drugs_com)^2/(2*(delt^2))))
  predict=sum(DrugResponse_com*wight,na.rm = T)/sum(wight,na.rm = T)
  return(predict)
}


# predict drug sensitivity using dual-layer model 
get_Predict=function(Processed.data,Training.parameter){
  
  TR.DrugResponse=Processed.data[[1]]
  TR.CellExpression=Processed.data[[2]]
  TR.DrugChemical=Processed.data[[3]]
  TR.drugs=colnames(TR.DrugResponse)
  TR.cells=rownames(TR.DrugResponse)
  
  OB.CellExpression=Processed.data[[4]]
  OB.DrugChemical=Processed.data[[5]]
  
  sigma=Training.parameter$sigma
  delt =Training.parameter$delt

  OB.drug=colnames(OB.DrugChemical)
  OB.cell=colnames(OB.CellExpression)
  
  if(OB.drug%in%TR.drugs){
    DrugName=OB.drug
    if(OB.cell%in%TR.cells){
      Upper_predict=get_upper_predict(TR.CellExpression,OB.CellExpression,sigma,DrugName,TR.DrugResponse)
      CellName=OB.cell
      Lower_predict=get_lower_predict(TR.DrugChemical,OB.DrugChemical,delt,CellName,TR.DrugResponse)
      lammda=Training.parameter$lammda[OB.drug]
      Dual.predict=(1-lammda)*Upper_predict+lammda*Lower_predict
      Predict=Dual.predict } else {
        Upper_predict=get_upper_predict(TR.CellExpression,OB.CellExpression,sigma,DrugName,TR.DrugResponse)
          Predict=Upper_predict}  
    } else {
        if(OB.cell%in%TR.cells){
          CellName=OB.cell
          Lower_predict=get_lower_predict(TR.DrugChemical,OB.DrugChemical,delt,CellName,TR.DrugResponse)
          Predict=Lower_predict
          } else {
          Predict=NA
        }}
  return(paste(OB.drug,OB.cell,Predict,sep = " "))
}


