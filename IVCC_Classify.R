# functions for IVCC classification

library(caret)
library(randomForest)
library(data.table)
library(bit64)
library(R2HTML)

TrainTest_SetBuilder=function(data,PositiveLabels,NegativeLabels,p.test.split){ # Building Training and Prediction Data for the Official Classifier
  #data - dataframe containing the full feature set for all instances
  #PositiveLabels - list of positive cases from data 
  #NegativeLabels - list of negative cases from data 
  #p.test.split - train/test split percentage
  
  
  
  temp=data
  #temp$creation_date=as.integer(temp$creation_date)
  row.names(temp)=make.names(temp$ScreenName,unique=TRUE)
  
  pos=temp[temp$ScreenName %in% PositiveLabels,]
  class=array(1,dim=dim(pos)[1])
  pos=cbind(class,pos)
  pos.test.ind=sample(1:dim(pos)[1],floor(dim(pos)[1]*p.test.split))
  pos.test=pos[pos.test.ind,]
  pos.train=pos[as.logical(1- (1:dim(pos)[1]) %in% pos.test.ind    ),]
  
  neg=temp[temp$ScreenName %in% NegativeLabels,]
  class=array(0,dim=dim(neg)[1])
  neg=cbind(class,neg)
  neg.test.ind=sample(1:dim(neg)[1],floor(dim(neg)[1]*p.test.split))
  neg.test=neg[neg.test.ind,]
  neg.train=neg[as.logical(1- (1:dim(neg)[1]) %in% neg.test.ind    ),]
  
  d.train=rbind(pos.train,neg.train)
  d.test=rbind(pos.test,neg.test)
  d.predict=temp[as.logical(1-(temp$ScreenName %in% c(d.train$ScreenName,d.test$ScreenName))),]
  return(list(d.train,d.test,d.predict))}

classifier=function(data_list,
                    algorithm,# 'svmRadial'  # 'lssvmPoly' 'lssvmLinear 'JRip' 'J48' 
                    feature_set, # as.factor(class)~followingCount + tweetCount + creation_date + ...
                    evaluation=TRUE,
                    metric,# 'Kappa'
                    ratio=c(.2,.8),
                    t){ 
  
  data_list[[1]]$class=as.factor(data_list[[1]]$class)
  data_list[[2]]$class=as.factor(data_list[[2]]$class)
  temp.features=c('class',feature_set)
  if(evaluation==TRUE)  train.data=data_list[[1]][,temp.features]
  if(evaluation==FALSE) train.data=rbind(data_list[[1]][,temp.features],data_list[[2]][,temp.features])
  
  if(evaluation==TRUE)  result.data=d[[2]][,temp.features]
  if(evaluation==FALSE) result.data=d[[3]][,temp.features[-1]]
  
  if(algorithm != 'randomForest') classifier=train(class ~ .,data = train.data,na.action=na.omit,method=algorithm,metric='Kappa')
  if(algorithm == 'randomForest') classifier=randomForest(class ~ .,data = train.data,na.action=na.omit,cutoff=ratio)
  
  if(evaluation==TRUE & algorithm != 'randomForest') test=predict.train(classifier,newdata=result.data,type='raw')
  if(evaluation==TRUE & algorithm == 'randomForest') test=predict(classifier,newdata=result.data,type='response')
  if(evaluation==TRUE) model.obj=confusionMatrix(test,d[[2]]$class)
  if(evaluation==TRUE){capture.output(cat(paste('ALGORITHM: ',algorithm,sep=''),'\n'),cat('\n'),
                                      cat(paste('FEATURE SET: ',paste(feature_set,collapse=', '),sep=''),'\n'),
                                      cat('\n'),cat(paste('TIME: ',as.character(Sys.time()),sep=''),'\n'),
                                      cat('\n'),model.obj,file=file.path('output',paste(t,'_model_eval.txt',sep='')))}
  if(evaluation==TRUE) return(model.obj)
  
  if(evaluation==FALSE & algorithm == 'randomForest' )  predicted=predict(classifier,newdata=result.data,type='response')
  if(evaluation==FALSE & algorithm != 'randomForest' )  predicted=predict.train(classifier,newdata=result.data,type='raw')
  if(evaluation==FALSE){                                
    predicted.names=d[[3]]$ScreenName[predicted==1]
    positive.sample=sample(d[[3]]$ScreenName[predicted==1],15)
    negative.sample=sample(d[[3]]$ScreenName[predicted==0],15)
    df=data.frame(c(
      paste(sum(as.numeric(predicted)-1),' predicted positive case, ',dim(d[[3]])[1]-sum(as.numeric(predicted)-1),' predicted negative case.',sep=''),
      '','Sample of Positive Cases',
      paste('<a href= "https://twitter.com/',positive.sample,'" target="_blank"> ',positive.sample,'</a>',sep=''),
      '',"Sample of Negative Cases",
      paste('<a href= "https://twitter.com/',negative.sample,'" target="_blank"> ',negative.sample,'</a>',sep='')),
      stringsAsFactors=FALSE)
    target <- HTMLInitFile(file.path('output'),filename=paste(t,'_Model_Validation',sep=''),
                           Title=paste(t,'_Model_Validation',sep=''))
    HTML(df,file=target,append=TRUE, Border=0, classtable='tablesorter',sortableDF=TRUE)
    HTMLEndFile()
  }
  if(evaluation==FALSE) return(predicted.names)
  
}