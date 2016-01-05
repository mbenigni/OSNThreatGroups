# author: Matthew Benigni
# date initiated: 1 March, 2014
# last update: 16 September, 2014
# This package facilitates iterative vertex clustering and classification (IVCC) as per 'add reference to paper'
# specifically this module develops the feature space for a clustering or classification task

library(igraph)
library(Matrix)
library(rARPACK)
library(data.table) # test comment
library(bit64)
library(dplyr)
library(wordcloud)
library(Cairo)
library(R2HTML)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(jsonlite)
library(sp)
library(rworldmap)

# buildFeatureSpace takes the edgelist and attribute output from the twitter_dm python library and builds a feature space
# for clustering or classification
buildFeatureSpace=function(output_path,
                           #seed_list,
                           attribute_path,
                           label_path,
                           lang_path,
                           friend_edge_path,
                           mention_edge_path, 
                           user_hashtag_path,
                           ht_alpha=.05,
                           user_alpha=.1,
                           hist_date='1970-01-01',
                           k=2){
  # tag simply provides a label for RData files outputted to disk
  # seed_list is the list of the userIDs for all seed agents used for the snowball search
  # attribute_path refers to the file path for the node attribute file
  # friend_edge_path refers to the file path for the friend edge list
  # mention_edge_path refers to teh file path for the mention edge list
  # hist_date a historical cut-off from which to trim mention edges in format '1970-01-01'
  # k defines the number of eigen vectors to extract from the graph laplatian, the default is 2 for 
  # classification problems, but should be larger for clustering applications
  l=ceiling(log2(k))
  hist_epoch=as.integer(as.POSIXct(as.Date(hist_date,origin='1970-01-01'),origin='1970-01-01'))
  attributes=fread(attribute_path,stringsAsFactors=FALSE,data.table=FALSE,integer64='numeric')
  if(anyDuplicated(attributes$userID)>0) attributes=attributes[-1*anyDuplicated(attributes$userID),]
  attributes=attributes[order(attributes$userID),]
  nameKey=data.frame(attributes[,c('userID','ScreenName')],stringsAsFactors=FALSE)
  temp.lang=fread(lang_path,stringsAsFactors=FALSE,integer64='numeric')
  temp.lang=cast(temp.lang, userID~lang,value='count')
  temp.lang[is.na(temp.lang)]=0
  temp.denom=array(apply(temp.lang[,-1],1,sum),dim=t(dim(temp.lang[,-1])))
  temp.lang[,-1]=temp.lang[,-1]/temp.denom
  col.lang=colnames(temp.lang)[-1]
  t=temp.lang[,-1]
  lang=sapply(1:dim(temp.lang)[1],function(x) col.lang[which.max(t[x,])])
  temp.lang=data.frame(temp.lang,lang)
  temp.lang=merGe(nameKey,temp.lang)
  urlRatio=attributes$urlCount/attributes$tweetCount
  urlRatio[is.infinite(urlRatio)]=0
  mentionRatio=attributes$mentionCount/attributes$tweetCount
  mentionRatio[is.infinite(mentionRatio)]=0
  attributes=data.frame(attributes[,1:7],urlRatio,mentionRatio,temp.lang,stringsAsFactors=FALSE)

  
  OutputSet=fread(label_path,stringsAsFactors=FALSE,colClasses=c('numeric','character',rep('logical',6)))
  isisIDs=as.character(OutputSet$userID[OutputSet$l.isis | OutputSet$p.isis])
  #load all node based features from networks
  # load friend edges and develop friend tie based features
  temp.edges=fread(friend_edge_path,data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Target %in% attributes$userID,]
  temp.edges=cbind(as.character(temp.edges$Source),as.character(temp.edges$Target))
  F=graph.edgelist(temp.edges)
  rm(temp.edges)
  FnodeAttributes=data.frame(userID=names(V(F)),
                             FinDegree=degree(F,v=V(F),mode='in'),
                             FoutDegree=degree(F,v=V(F),mode='out'),
                             stringsAsFactors=FALSE)
  FnodeAttributes=merGe(nameKey,FnodeAttributes)
  FiP=FnodeAttributes$FoutDegree/attributes$followingCount
  FiP[is.na(FiP)]=0
  FiP[is.infinite(FiP)]=0
  evRF=getEigenVectors(input_network=F,input_network_name='RF',k=l,w='LM',max.iter=9000000,nameKey)
  FnodeAttributes=data.frame(FnodeAttributes,FiP,evRF[,-1],stringsAsFactors=FALSE)
  rm(F,evRF)
  
  #load mention edges and develop mention based features
  temp.edges=fread(mention_edge_path,data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Target %in% attributes$userID,]
  if(hist_date != '1970-01-01') temp.edges=trim_by_epoch(temp.edges,hist_epoch)
  M=graph.edgelist(cbind(as.character(temp.edges$Source),as.character(temp.edges$Target)))
  rm(temp.edges)
  MnodeAttributes=data.frame(names(V(M)),degree(M,v=V(M),mode='in'),degree(M,v=V(M),mode='out'),stringsAsFactors=FALSE)
  colnames(MnodeAttributes)=c('userID','MinDegree','MoutDegree')
  MnodeAttributes=merGe(nameKey,MnodeAttributes)
  
  evRM=getEigenVectors(input_network=M,input_network_name='RM',k=l,w='LM',max.iter=9000000,nameKey)
  MnodeAttributes=data.frame(MnodeAttributes,evRM[,-1],stringsAsFactors=FALSE)
  rm(M,evRM)
  
  # load user x hashtag edges and develop features

  temp.edges=fread(user_hashtag_path,stringsAsFactors=FALSE,encoding='UTF-8',
                   data.table=FALSE,integer64='numeric')
  temp.edges$Target=iconv(temp.edges$Target,to='UTF-8')
  ht_freq=as.data.frame(ftable(temp.edges$Target),stringsAsFactors=FALSE)
  if(ht_alpha>0){ht_bounds=quantile(ht_freq$Freq,c(ht_alpha,1-ht_alpha))
  ht.list=as.character(ht_freq$Var1)[ht_freq$Freq>ht_bounds[1]&ht_freq$Freq<ht_bounds[2]]}
  
  if(user_alpha>0){user_freq=as.data.frame(ftable(temp.edges$Source),stringsAsFactors=FALSE)
  user_bounds=quantile(user_freq$Freq,1-user_alpha)
  user.list=as.character(user_freq$Var1)[user_freq$Freq<user_bounds]}
  
  if(ht_alpha>0){temp.edges=temp.edges[temp.edges$Target %in% ht.list,]}
  
  if(user_alpha>0){temp.edges=temp.edges[as.character(temp.edges$Source) %in% user.list,]}
  
  if(hist_date != '1970-01-01') temp.edges=trim_by_epoch(temp.edges,hist_epoch)
  temp.edges$Target=paste('#',temp.edges$Target,sep='')
  temp.edges$Source=as.character(temp.edges$Source)
  
  UHT=graph.edgelist(cbind(as.character(temp.edges$Source),temp.edges$Target))
  
  A.uht=get.adjacency(UHT)
  mapping=bipartite.mapping(UHT)$type
  A=A.uht[,mapping]
  A=A[as.logical(1-mapping),]
  HTcount=rowSums(A)
  D1=Diagonal(n=dim(A)[1],x=1/sqrt(rowSums(A)))
  D2=Diagonal(n=dim(A)[2],x=1/sqrt(colSums(A)))
  An=D1%*%A%*%D2
  obj=svds(An,l,nu=l,nv=0,ncv=l+5)
  vecMat=as.matrix(D1%*%obj$u)
  user.Vectors=data.frame(rownames(An),HTcount,vecMat,stringsAsFactors=FALSE)
  colnames(user.Vectors)=c('userID','HTcount',paste('UHTev',1:l,sep=''))
  UHTnodeAttributes=merGe(nameKey,user.Vectors)
  
  
  FeatureSet=cbind(attributes,FnodeAttributes[,-1],MnodeAttributes[,-1],UHTnodeAttributes[,-1])
  FeatureSet[is.na(FeatureSet)]=0
  FeatureSet=FeatureSet[complete.cases(FeatureSet),]
  feature_set=names(FeatureSet)[!(names(FeatureSet) %in% c('userID','ScreenName','lang'))]
  write.csv(feature_set,file=file.path('input','featureSetNames.csv'),row.names=FALSE)
  return(FeatureSet)
}



# builds feature space for final isis report
bispectralClustering=function(attribute_path=file.path('input','attribute.tsv'),
                              label_path=file.path('input','OutputSet_2015-11-19.csv'),
                              lang_path=file.path('input','langfile.tsv'),
                              friend_edge_path=file.path('input','friend_edgefile.tsv'),
                              mention_edge_path=file.path('input','mention_edgefile.tsv'),
                              user_hashtag_path=file.path('input','user_ht_edgefile.tsv'),
                              output_path=file.path('input','output'),
                              ht_alpha=.05,
                              user_alpha=.1,
                              hist_date='1970-01-01',
                              k=100){
  # tag simply provides a label for RData files outputted to disk
  # seed_list is the list of the userIDs for all seed agents used for the snowball search
  # attribute_path refers to the file path for the node attribute file
  # friend_edge_path refers to the file path for the friend edge list
  # mention_edge_path refers to teh file path for the mention edge list
  # hist_date a historical cut-off from which to trim mention edges in format '1970-01-01'
  # k defines the number of eigen vectors to extract from the graph laplatian, the default is 2 for 
  # classification problems, but should be larger for clustering applications
  l=ceiling(log2(k))
  hist_epoch=as.integer(as.POSIXct(as.Date(hist_date,origin='1970-01-01'),origin='1970-01-01'))
  attributes=fread(attribute_path,stringsAsFactors=FALSE,data.table=FALSE,integer64='numeric')
  posIDs=fread(label_path,stringsAsFactors=FALSE,data.table=FALSE,integer64='numeric')
  if(anyDuplicated(attributes$userID)>0){attributes=attributes[-1*anyDuplicated(attributes$userID),]}
  attributes=attributes[order(attributes$userID),]
  date=matrix(unlist(sapply(attributes$creation_date,strsplit,split=" ")),nrow=dim(attributes)[1],byrow=T)[,1]
  attributes$creation_date=as.Date(date,origin='1070-01-01')
  date=matrix(unlist(strsplit(attributes$lastTweet,split=" ")),ncol=2,byrow=TRUE)[,1]
  attributes$lastTweet=as.Date(date,origin='1070-01-01')
  nameKey=data.frame(attributes[,c('userID','ScreenName')],stringsAsFactors=FALSE)
  temp.lang=fread(lang_path,stringsAsFactors=FALSE,integer64='numeric')
  temp.lang=cast(temp.lang, userID~lang,value='count')
  temp.lang[is.na(temp.lang)]=0
  
  
  png(file=file.path(shiny_path,'www','Jihadist_Langage_Pie.png'),res=150,width=800)
  data=apply(temp.lang,2,sum)
  df=data.frame(lang=names(data),value=data)
  df=df[order(df$value,decreasing=TRUE),]
  df$lang=paste(1:length(df$lang),'. ',df$lang,sep='')
  df_low=data.frame(lang='other',value=sum(df$value[9:dim(df)[1]]))
  df=rBind(df[1:8,],df_low)
  bp<- ggplot(df, aes(x="", y=value, fill=lang))+
    
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0)+
    theme(axis.text.x=element_blank()) +
    ggtitle('Jihadist Tweet Volume by Language')
  pie
  dev.off()
  
  temp.denom=array(apply(temp.lang[,-1],1,sum),dim=t(dim(temp.lang[,-1])))
  temp.lang[,-1]=temp.lang[,-1]/temp.denom
  col.lang=colnames(temp.lang)[-1]
  t=temp.lang[,-1]
  lang=sapply(1:dim(temp.lang)[1],function(x) col.lang[which.max(t[x,])])
  temp.lang=data.frame(temp.lang,lang,stringsAsFactors=FALSE)
  temp.lang=merGe(nameKey,temp.lang)
  urlRatio=attributes$urlCount/attributes$tweetCount
  urlRatio[is.infinite(urlRatio)]=0
  mentionRatio=attributes$mentionCount/attributes$tweetCount
  mentionRatio[is.infinite(mentionRatio)]=0
  attributes=data.frame(attributes[,1:7],urlRatio,mentionRatio,temp.lang,stringsAsFactors=FALSE)
  
  
  
  #load all node based features from networks
  
  # load friend edges and develop friend tie based features
  temp.edges=fread(friend_edge_path,data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Source %in% attributes$userID,]
  temp.edges=temp.edges[temp.edges$Target %in% attributes$userID,]
  temp.edges=cbind(as.character(temp.edges$Source),as.character(temp.edges$Target))
  F=graph_from_edgelist(temp.edges,directed=TRUE)
  A.f=get.adjacency(F)
  PosFollows=rowSums(A.f[,colnames(A.f) %in% as.character(posIDs$V1)])
  RF=as.undirected(F,mode='mutual')
  RFvec=data.frame(userID=names(V(RF)),RFDegree=degree(RF,V(RF),mode='in'),stringsAsFactors=FALSE)
  RFvec=merGe(nameKey,RFvec)
  rm(temp.edges)
  FnodeAttributes=data.frame(userID=names(V(F)),
                             FinDegree=degree(F,v=V(F),mode='in'),
                             FoutDegree=degree(F,v=V(F),mode='out'),
                             JihadFollowing=PosFollows,
                             stringsAsFactors=FALSE)
  FnodeAttributes=merGe(nameKey,FnodeAttributes)
  FnodeAttributes=data.frame(FnodeAttributes,RFDegree=RFvec[,2],stringsAsFactors=FALSE)
  FnodeAttributes$JihadFollowing=FnodeAttributes$JihadFollowing/attributes$followingCount
  FnodeAttributes$JihadFollowing[is.na(FnodeAttributes$JihadFollowing)]=0
  # wrong numbers should be between zero and one 
  rm(F)
  
  
  #load mention edges and develop mention based features
  
  temp.edges=fread(mention_edge_path,data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Source %in% attributes$userID,]
  if(hist_date != '1970-01-01') temp.edges=trim_by_epoch(temp.edges,hist_epoch)
  temp.edges=temp.edges[temp.edges$Target %in% attributes$userID,]
  M=graph.edgelist(cbind(as.character(temp.edges$Source),as.character(temp.edges$Target)))
  RM=as.undirected(M,mode='mutual')
  RMvec=data.frame(userID=names(V(RM)),RMDegree=degree(RM,V(RM),mode='total'),stringsAsFactors=FALSE)
  RMvec=merGe(nameKey,RMvec)
  rm(temp.edges)
  MnodeAttributes=data.frame(userID=names(V(M)),
                             MinDegree=degree(M,v=V(M),mode='in'),
                             MoutDegree=degree(M,v=V(M),mode='out'),
                             stringsAsFactors=FALSE)
  
  MnodeAttributes=merGe(nameKey,MnodeAttributes)
  MnodeAttributes=data.frame(MnodeAttributes,RMDegree=RMvec[,2],stringsAsFactors=FALSE)
  rm(M)
  
  
  # load user x hashtag edges and develop features
  temp.edges=fread(user_hashtag_path,stringsAsFactors=FALSE,encoding='UTF-8',
                   data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Source %in% attributes$userID,]
  temp.edges$Target=iconv(temp.edges$Target,to='UTF-8')
  ht_freq=as.data.frame(ftable(temp.edges$Target),stringsAsFactors=FALSE)
  if(ht_alpha>0){ht_bounds=quantile(ht_freq$Freq,c(ht_alpha,1-ht_alpha))
                 ht.list=as.character(ht_freq$Var1)[ht_freq$Freq>ht_bounds[1]&ht_freq$Freq<ht_bounds[2]]}
  
  if(user_alpha>0){user_freq=as.data.frame(ftable(temp.edges$Source),stringsAsFactors=FALSE)
                   user_bounds=quantile(user_freq$Freq,1-user_alpha)
                   user.list=as.character(user_freq$Var1)[user_freq$Freq<user_bounds]}
  
  if(ht_alpha>0){temp.edges=temp.edges[temp.edges$Target %in% ht.list,]}
  
  if(user_alpha>0){temp.edges=temp.edges[as.character(temp.edges$Source) %in% user.list,]}
  
  if(hist_date != '1970-01-01') temp.edges=trim_by_epoch(temp.edges,hist_epoch)
  temp.edges$Target=paste('#',temp.edges$Target,sep='')
  temp.edges$Source=as.character(temp.edges$Source)
  
  UHT=graph.edgelist(cbind(as.character(temp.edges$Source),temp.edges$Target))
  
  A.uht=get.adjacency(UHT)
  mapping=bipartite.mapping(UHT)$type
  A=A.uht[,mapping]
  A=A[as.logical(1-mapping),]
  Ucount=rowSums(A)
  HTcount=colSums(A)
  D1=Diagonal(n=dim(A)[1],x=1/sqrt(rowSums(A)))
  D2=Diagonal(n=dim(A)[2],x=1/sqrt(colSums(A)))
  An=D1%*%A%*%D2
  obj=svds(An,l,nu=l,nv=l,ncv=l+5)
  uMat=data.frame(rownames(An),Ucount,as.matrix(D1%*%obj$u),stringsAsFactors=FALSE)
  htMat=data.frame(colnames(An),HTcount,as.matrix(D2%*%obj$v),stringsAsFactors=FALSE)
  colnames(uMat)=colnames(htMat)=c('userID','degree',paste('UHTev',1:l,sep=''))
  uhtMat=as.matrix(rBind(uMat,htMat))
  
  kobj=kmeans(uhtMat[,c(-1,-2)],k,iter.max=10000,algorithm='Lloyd')
  uhtDF=data.frame(uhtMat,kobj$cluster,stringsAsFactors=FALSE)
  colnames(uhtDF)[dim(uhtDF)[2]]='cluster'
  
  
  target.indecies=which(uhtDF$userID %in% as.character(posIDs$V1) )
  userOutput=which(!(grepl('#',uhtDF$userID)))
  htOutput=which((grepl('#',uhtDF$userID)))
  
  userDF=uhtDF[userOutput,]
  userDF$userID=as.integer64(userDF$userID)
  htDF=uhtDF[htOutput,]
  
  
  cluster=1:k
  target.clusters=kobj$cluster[target.indecies]
  tcount=sapply(1:k,function(x) sum(target.clusters==x))
  user.clusters=kobj$cluster[userOutput]
  size=sapply(1:k,function(x) sum(user.clusters==x))
  ntcount=size-tcount
  clusterMat=data.frame(cluster,ntcount,tcount,size)
  clusterMat=clusterMat[clusterMat$size>0,]
  clusterMat=clusterMat[clusterMat$tcount>0,]
  clusterMat=clusterMat[order(clusterMat$tcount,decreasing=TRUE),]
  n=20
  barMat=as.matrix(clusterMat[1:n,])
  barMat=barMat[order(barMat[,3]),]
  # barplot for taget id population by cluster
  png(file.path(shiny_path,'www','cluster_barplot.png'),height=600,width=1000,res=150,pointsize=10)
  barplot(t(barMat[,2:3]),names.arg=barMat[,1],horiz=TRUE,las=2,ylab='cluster',
          xlab='number of users',
          main=paste('Top ',n,' of ',k,
                     ' pro-ISIS Topical Clusters wrt validated ISIS supporters: \n light grey indicates manually validated ISIS supporters',
                     sep=''))
  dev.off()
  # id clusters of interest and remove members of those clusters
  htMat=htDF[,c('userID','degree','cluster')]
  htMat$userID=gsub('#','',htMat$userID)
  htMat$userID=gsub('_',' ',htMat$userID)
  htList=lapply(1:k,htDataCleaner,htMat)
  hashtags='hashtag'
  degree='degree'
  cluster='cluster'
  for(i in 1:k){hashtags=c(hashtags,htList[[i]]$hashtag)
                degree=c(degree,htList[[i]]$degree)
                cluster=c(cluster,array(i,dim=dim(htList[[i]])[1]))
  }
  htData=data.frame(hashtags,degree,cluster)
  write.csv(htData,
            file=file.path(output_path,'hashtags_for_translation.csv'),quote=FALSE,row.names=FALSE)
  
  
  userDFtemp=userDF[,c('userID','cluster')]
  HTclusterAttribute=merGe(nameKey,userDFtemp)
  features=data.frame(attributes,FnodeAttributes[,-1],MnodeAttributes[,-1],topic_cluster=HTclusterAttribute[,-1])
  save(features,htData,file=file.path(output_path,paste('isisFeatures_',Sys.Date(),'.RData',sep='')))
  userData=features[,c('ScreenName','topic_cluster','JihadFollowing','RMDegree','RFDegree','lang','ar','en','ps','ur','followingCount','followerCount','tweetCount',
                       'urlRatio','mentionRatio')]
  userData$ScreenName=hyperlink(userData$ScreenName)
  save(userData,file=file.path(shiny_path,'data','userData.RData'))
}



# builds feature space for final isis report
reportFeatureSpace=function(attribute_path=file.path('input','attribute.tsv'),
                            current_label_path_path=file.path('input','OutputSet_2015-11-19.csv'),
                            hist_label_path=file.path('input','OutputSet_2014-11-24.csv'),
                            lang_path=file.path('input','langfile.tsv'),
                            friend_edge_path=file.path('input','friend_edgefile.tsv'),
                            mention_edge_path=file.path('input','mention_edgefile.tsv'),
                            user_hashtag_path=file.path('input','user_ht_edgefile.tsv'),
                            output_path=file.path('input','output'),
                            ht_alpha=.05,
                            user_alpha=.1,
                            hist_date='1970-01-01',
                            k_ht=150,
                            k_m=300){
  # tag simply provides a label for RData files outputted to disk
  # seed_list is the list of the userIDs for all seed agents used for the snowball search
  # attribute_path refers to the file path for the node attribute file
  # friend_edge_path refers to the file path for the friend edge list
  # mention_edge_path refers to teh file path for the mention edge list
  # hist_date a historical cut-off from which to trim mention edges in format '1970-01-01'
  # k defines the number of eigen vectors to extract from the graph laplatian, the default is 2 for 
  # classification problems, but should be larger for clustering applications
  l_ht=ceiling(log2(k_ht))
  l_m=ceiling(log2(k_m))
  hist_epoch=as.integer(as.POSIXct(as.Date(hist_date,origin='1970-01-01'),origin='1970-01-01'))
  attributes=fread(attribute_path,stringsAsFactors=FALSE,data.table=FALSE,integer64='numeric')
  attributes=attributes[order(attributes$userID),]
  OutputSet_c=fread(current_label_path,stringsAsFactors=FALSE,data.table=FALSE,integer64='numeric')
  OutputSet_h=fread(hist_label_path,stringsAsFactors=FALSE,data.table=FALSE,integer64='numeric')
  isis1=OutputSet_c$userID[OutputSet_c$p.isis|OutputSet_c$l.isis]
  isis2=OutputSet_h$userID[OutputSet_h$p.isis|OutputSet_h$l.isis]
  isisIDs=unique(c(isis1,isis2))
  
  if(anyDuplicated(attributes$userID)>0){attributes=attributes[-1*anyDuplicated(attributes$userID),]}
  date=matrix(unlist(sapply(attributes$creation_date,strsplit,split=" ")),nrow=dim(attributes)[1],byrow=T)[,1]
  attributes$creation_date=as.Date(date,origin='1070-01-01')
  date=matrix(unlist(strsplit(attributes$lastTweet,split=" ")),ncol=2,byrow=TRUE)[,1]
  attributes$lastTweet=as.Date(date,origin='1070-01-01')
  nameKey=data.frame(attributes[,c('userID','ScreenName')],stringsAsFactors=FALSE)
  temp.lang=fread(lang_path,stringsAsFactors=FALSE,integer64='numeric',data.table=FALSE)
  temp.lang=cast(temp.lang, userID~lang,value='count')
  temp.lang[is.na(temp.lang)]=0
  
  
  png(file=file.path(shiny_path,'www','Jihadist_Langage_Pie.png'),res=150,width=800)
  data=apply(temp.lang,2,sum)
  df=data.frame(lang=names(data),value=data)
  df=df[order(df$value,decreasing=TRUE),]
  df$lang=paste(1:length(df$lang),'. ',df$lang,sep='')
  df_low=data.frame(lang='other',value=sum(df$value[9:dim(df)[1]]))
  df=rBind(df[1:8,],df_low)
  bp<- ggplot(df, aes(x="", y=value, fill=lang))+
    
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0)+
    theme(axis.text.x=element_blank()) +
    ggtitle('Jihadist Tweet Volume by Language')
  pie
  dev.off()
  
  temp.denom=array(apply(temp.lang[,-1],1,sum),dim=t(dim(temp.lang[,-1])))
  temp.lang[,-1]=temp.lang[,-1]/temp.denom
  col.lang=colnames(temp.lang)[-1]
  t=temp.lang[,-1]
  lang=sapply(1:dim(temp.lang)[1],function(x) col.lang[which.max(t[x,])])
  temp.lang=data.frame(temp.lang,lang,stringsAsFactors=FALSE)
  temp.lang=merGe(nameKey,temp.lang)
  urlRatio=attributes$urlCount/attributes$tweetCount
  urlRatio[is.infinite(urlRatio)]=0
  mentionRatio=attributes$mentionCount/attributes$tweetCount
  mentionRatio[is.infinite(mentionRatio)]=0
  attributes=data.frame(attributes[,1:7],urlRatio,mentionRatio,temp.lang,stringsAsFactors=FALSE)
  
  
  
  #load all node based features from networks
  
  # load friend edges and develop friend tie based features
  temp.edges=fread(friend_edge_path,data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Source %in% attributes$userID,]
  mostFollowed=as.data.table(ftable(temp.edges$Target))
  colnames(mostFollowed)=c('userID','isisFollowed')
  mostFollowed=mostFollowed[!(mostFollowed$userID %in% OutputSet_c$userID),]
  mostFollowed=mostFollowed[order(mostFollowed$isisFollowed,decreasing=TRUE),]
  
  temp.edges=temp.edges[temp.edges$Target %in% attributes$userID,]
  write.csv(temp.edges,file=file.path(output_path,paste('isis_friend_edgefile_',Sys.Date(),'.csv',sep='')),
            quote=FALSE,row.names=FALSE)
  temp.edges=cbind(as.character(temp.edges$Source),as.character(temp.edges$Target))
  F=graph_from_edgelist(temp.edges,directed=TRUE)
  RF=as.undirected(F,mode='mutual')
  RFvec=data.frame(userID=names(V(RF)),RFDegree=degree(RF,V(RF),mode='in'),stringsAsFactors=FALSE)
  RFvec=merGe(nameKey,RFvec)
  rm(temp.edges)
  FnodeAttributes=data.frame(userID=names(V(F)),
                             FinDegree=degree(F,v=V(F),mode='in'),
                             FoutDegree=degree(F,v=V(F),mode='out'),
                             stringsAsFactors=FALSE)
  FnodeAttributes=merGe(nameKey,FnodeAttributes)
  FiP=FnodeAttributes$FoutDegree/attributes$followingCount
  FiP[is.na(FiP)]=0
  FiP[is.infinite(FiP)]=0
  FnodeAttributes=data.frame(FnodeAttributes,FiP,RFDegree=RFvec[,2],stringsAsFactors=FALSE)
  # wrong numbers should be between zero and one 
  rm(F)
  
  
  #load mention edges and develop mention based features
  
  temp.edges=fread(mention_edge_path,data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Source %in% attributes$userID,]
  if(hist_date != '1970-01-01') temp.edges=trim_by_epoch(temp.edges,hist_epoch)
  
  mostMentioned=as.data.table(ftable(temp.edges$Target))
  colnames(mostMentioned)=c('userID','isisMentioned')
  mostMentioned=mostMentioned[!(mostMentioned$userID %in% OutputSet_c$userID),]
  mostMentioned=mostMentioned[order(mostMentioned$isisMentioned,decreasing=TRUE),]
  
  searchList=merge(mostMentioned,mostFollowed,by='userID')
  searchList=searchList[order(searchList$isisFollowed,decreasing=TRUE),]
  write.csv(searchList,file=file.path(output_path,paste('searchList_',Sys.Date(),'.csv',sep='')),
            row.names=FALSE,quote=FALSE)
  mention.edges=temp.edges
  temp.edges=temp.edges[temp.edges$Target %in% attributes$userID,]
  write.csv(temp.edges,file=file.path(output_path,paste('isis_mention_edgefile_',Sys.Date(),'.csv',sep='')),
            quote=FALSE,row.names=FALSE)
  
  M=graph.edgelist(cbind(as.character(temp.edges$Source),as.character(temp.edges$Target)))
  evRM=getEigenVectors(input_network=M,input_network_name='RM',k=l_m,w='SR',max.iter=9000000,nameKey)
  
  
 
  MnodeAttributes=data.frame(userID=names(V(M)),
                             MinDegree=degree(M,v=V(M),mode='in'),
                             MoutDegree=degree(M,v=V(M),mode='out'),
                             stringsAsFactors=FALSE)
  
  MnodeAttributes=merGe(nameKey,MnodeAttributes)
  MnodeAttributes=data.frame(MnodeAttributes,evRM[,-1],stringsAsFactors=FALSE)
  mention_kobj=kmeans(evRM[,c(-1,-1*dim(evRM)[2])],k_m,iter.max=10000,algorithm='Lloyd')
  MnodeAttributes=data.frame(MnodeAttributes[,c(1:3,(dim(MnodeAttributes)[2]-1):dim(MnodeAttributes)[2])],mention_cluster=mention_kobj$cluster)
  
  rm(M)
  
  
  # load user x hashtag edges and develop features
  temp.edges=fread(user_hashtag_path,stringsAsFactors=FALSE,encoding='UTF-8',
                   data.table=FALSE,integer64='numeric')
  temp.edges=temp.edges[temp.edges$Source %in% attributes$userID,]
  temp.edges$Target=iconv(temp.edges$Target,to='UTF-8')
  ht_freq=as.data.frame(ftable(temp.edges$Target),stringsAsFactors=FALSE)
  if(ht_alpha>0){ht_bounds=quantile(ht_freq$Freq,c(ht_alpha,1-ht_alpha))
                 ht.list=as.character(ht_freq$Var1)[ht_freq$Freq>ht_bounds[1]&ht_freq$Freq<ht_bounds[2]]}
  
  if(user_alpha>0){user_freq=as.data.frame(ftable(temp.edges$Source),stringsAsFactors=FALSE)
                   user_bounds=quantile(user_freq$Freq,1-user_alpha)
                   user.list=as.character(user_freq$Var1)[user_freq$Freq<user_bounds]}
  
  if(ht_alpha>0){temp.edges=temp.edges[temp.edges$Target %in% ht.list,]}
  
  if(user_alpha>0){temp.edges=temp.edges[as.character(temp.edges$Source) %in% user.list,]}
  
  if(hist_date != '1970-01-01') temp.edges=trim_by_epoch(temp.edges,hist_epoch)
  ht.edges=temp.edges
  temp.edges$Target=paste('#',temp.edges$Target,sep='')
  temp.edges$Source=as.character(temp.edges$Source)
  
  UHT=graph.edgelist(cbind(as.character(temp.edges$Source),temp.edges$Target))
  
  A.uht=get.adjacency(UHT)
  mapping=bipartite.mapping(UHT)$type
  A=A.uht[,mapping]
  A=A[as.logical(1-mapping),]
  Ucount=rowSums(A)
  HTcount=colSums(A)
  D1=Diagonal(n=dim(A)[1],x=1/sqrt(rowSums(A)))
  D2=Diagonal(n=dim(A)[2],x=1/sqrt(colSums(A)))
  An=D1%*%A%*%D2
  obj=svds(An,l_ht,nu=l_ht,nv=l_ht,ncv=l_ht+5)
  uMat=data.frame(rownames(An),Ucount,as.matrix(D1%*%obj$u),stringsAsFactors=FALSE)
  htMat=data.frame(colnames(An),HTcount,as.matrix(D2%*%obj$v),stringsAsFactors=FALSE)
  colnames(uMat)=colnames(htMat)=c('userID','degree',paste('UHTev',1:l_ht,sep=''))
  uhtMat=as.matrix(rBind(uMat,htMat))
  
  ht_kobj=kmeans(uhtMat[,c(-1,-2)],k_ht,iter.max=10000,algorithm='Lloyd')
  uhtDF=data.frame(uhtMat,ht_kobj$cluster,stringsAsFactors=FALSE)
  colnames(uhtDF)[dim(uhtDF)[2]]='topic_cluster'
  
  target.ids=unique(c(OutputSet_c$userID[OutputSet_c$l.isis],OutputSet_h$userID[OutputSet_h$l.isis]))
  target.indecies=which(uhtDF$userID %in% as.character(target.ids) )
  userOutput=which(!(grepl('#',uhtDF$userID)))
  
  
  userDF=uhtDF[userOutput,]
  userDF$userID=as.integer64(userDF$userID)
  userDF=userDF[,c('userID','topic_cluster')]
  HTclusterAttribute=merGe(nameKey,userDF)
  HTclusterAttribute=data.frame(HTclusterAttribute,targetIDs=HTclusterAttribute$userID %in% target.ids)
  ht_tgtByCluster=sapply(1:k_ht,function(x) sum(HTclusterAttribute$targetIDs[HTclusterAttribute$topic_cluster==x]))
  ht_ClusterSize=sapply(1:k_ht,function(x) sum(HTclusterAttribute$topic_cluster==x))
  ht_tgtPct=ht_tgtByCluster/ht_ClusterSize
  htBarMat=data.frame(cluster=1:k_ht,unlabelled=ht_ClusterSize-ht_tgtByCluster,labelled=ht_tgtByCluster,ht_tgtPct)
  htBarMat=htBarMat[order(htBarMat$ht_tgtPct),]
  d_ht=dim(htBarMat)[1]
  n=20
  
  Mtemp=MnodeAttributes[MnodeAttributes$userID %in% target.ids,]
  mention_tgtByCluster=sapply(1:k_m,function(x) sum(Mtemp$mention_cluster==x))
  mention_tgtPct=mention_tgtByCluster/mention_kobj$size
  mentionBarMat=data.frame(cluster=1:k_m,unlabelled=mention_kobj$size-mention_tgtByCluster,labelled=mention_tgtByCluster,mention_tgtPct)
  mentionBarMat=mentionBarMat[order(mentionBarMat$mention_tgtPct),]
  d_mention=dim(mentionBarMat)[1]
  

  # barplot for taget id population by cluster
  png(file.path(shiny_path,'www','cluster_barplot.png'),height=1500,width=1000,res=150,pointsize=10)
  par(mfcol=c(2,1))
  barplot(t(htBarMat[,2:3][(d_ht-n):d_ht,]),horiz=TRUE,las=2,main='Hashtag Clusters Most Likely to be Jihadist',xlab='number of users',ylab='cluster')
  barplot(t(mentionBarMat[,2:3][(d_mention-n):d_mention,]),horiz=TRUE,las=2,main='Mention Clusters Most Likely to be Jihadist',xlab='number of users',ylab='cluster')
  dev.off()
  # id clusters of interest and remove members of those clusters
  htOutput=which((grepl('#',uhtDF$userID)))
  htDF=uhtDF[htOutput,]
  
  htMat=htDF[,c('userID','degree','topic_cluster')]
  htMat$userID=gsub('#','',htMat$userID)
  htMat$userID=gsub('_',' ',htMat$userID)
  htList=lapply(1:k_ht,htDataCleaner,htMat)
  hashtags='hashtag'
  degree='degree'
  cluster='cluster'
  for(i in 1:k_ht){hashtags=c(hashtags,htList[[i]]$hashtag)
                degree=c(degree,htList[[i]]$degree)
                cluster=c(cluster,array(i,dim=dim(htList[[i]])[1]))
  }
  htData=data.frame(hashtags,degree,cluster)
  write.csv(htData,
            file=file.path(output_path,'hashtags_for_translation.csv'),quote=FALSE,row.names=FALSE)
  
  
  rmecdf=ecdf(MnodeAttributes$RMDegree)
  features=data.frame(attributes,FnodeAttributes[,-1],MnodeAttributes[,-1],HTclusterAttribute[,-1],jihadInfluence=FnodeAttributes$FiP+rmecdf(MnodeAttributes$RMDegree))
  save(features,htData,ht.edges,mention.edges,file=file.path(output_path,paste('isisFeatures_',Sys.Date(),'.RData',sep='')))
  
  userData=features[,c('ScreenName','topic_cluster','mention_cluster','jihadInfluence','FiP','RMDegree','RFDegree','lang','ar','en','ps','ur','fr','followingCount','followerCount','tweetCount',
                   'urlRatio','mentionRatio')]
  userData$ScreenName=hyperlink(userData$ScreenName)
  write.csv(features[,c('userID','ScreenName','topic_cluster','mention_cluster','jihadInfluence','FiP','RMDegree','RFDegree','lang','ar','en','ps','ur','fr','followingCount','followerCount','tweetCount',
                        'urlRatio','mentionRatio')],file=file.path(output_path,'user_attribute.csv'),quote=FALSE,row.names=FALSE)
  save(userData,features,k_ht,k_m,htMat,ht.edges,mention.edges,file=file.path(shiny_path,'data','userData.RData'))
  }



# trim_by_Date allows the user to remove edges that are older than a given date
trim_by_Date=function(edge.list,date){
  edge.list$date=as.Date(edge.list$date,origin='1970-01-01')
  edge.list=edge.list[edge.list$date>=as.Date(date,origin='1970-01-01'),]
  return(edge.list)}

trim_by_epoch=function(edge.list,epoch){
  edge.list=edge.list[edge.list$date>=epoch,]
  return(edge.list)}

trim_by_2_dates=function(edge.list,s_epoch,e_epoch){
  edge.list=edge.list[edge.list$date>=s_epoch & edge.list$date<=e_epoch ,]
  return(edge.list)}


hyperlink=function(screenNames){
  linkList=paste('\'<a href= "https://twitter.com/',screenNames,'" target="_blank"> ',screenNames,'</a>\'',sep='')
  return(linkList)}

# extracts eigen vectors for spectral clustering
getEigenVectors=function(input_network,input_network_name,k=2,w='SR',max.iter=9000000,nameKey){        
  input_network=as.undirected(input_network,mode='mutual')
  L=graph.laplacian(input_network,normalized=TRUE,sparse=getIgraphOpt("sparsematrices"))
  ev <- eigs(L,k,which=w,lower=TRUE,ncv=k+5,opts=list(maxitr=max.iter))
  total.degree=degree(input_network,v=V(input_network),mode='total')
  temp=data.frame(rownames(L),ev$vectors,total.degree,stringsAsFactors=FALSE)
  names(temp)[1]='userID'
  temp$userID=as.numeric(temp$userID)
  temp=merGe(nameKey,temp)
  temp[is.na(temp)]=0
  colnames(temp)=c('userID',paste(input_network_name,'ev',1:k,sep=''),paste(input_network_name,'Degree',sep=''))
  return(temp)
}

clusterExplorer=function(cluster,features,htDF,n){
  pG=brewer.pal(9,"Greys")[6:8]
  features=features[features$cluster==cluster,]
  features=features[,c(2:8,dim(features)[2])]
  features$ScreenName=paste('<a href= "https://twitter.com/',features$ScreenName,'" target="_blank"> ',features$ScreenName,'</a>',sep='')
  hashtags=htDF[htDF$cluster==cluster,]
  png(file.path('output_html','webpage','images',paste('clusterCloud',cluster,'.png',sep='')),pointsize=6,res=300,bg='transparent',width=1000,height=800)
  cex_limit=4
  wordcloud(substring(hashtags$userID,2)[1:n],as.numeric(hashtags$degree[1:n]),bg='transparent',
            use.r.layout=FALSE,scale=c(1.5,.15),rot.per=.15,colors=pG)
  dev.off()
  
  output=features[order(features$followerCount,decreasing=TRUE),]
  outputJSON=toJSON(output)
  write(outputJSON,file=file.path('output_html','webpage','clusterExplorer.json'))
  target <- HTMLInitFile(file.path('output_html','webpage'),
                         filename=paste('ClusterExplorer_',cluster,sep=''),
                         Title='User x Hash Tag Bispectral Co-Cluster Explorer')
  HTML(output[sample(1:dim(output)[1],30),],file=target,append=TRUE, Border=0, classtable='tablesorter',sortableDF=TRUE)
  HTMLEndFile()
}

attributeExplorer=function(attribute,value,features,htDF){
  col.index=which(colnames(features)==attribute)
  features=features[features[,col.index]==value,]
  names=features$ScreenName
  nlinks=paste('<a href= "https://twitter.com/',names,'" target="_blank"> ',names,'</a>',sep='')
  output=data.frame(nlinks,features[,3:8])
  output[,2:4]=sapply(output[,2:4],as.numeric)
  output=output[order(output$followerCount,decreasing=TRUE),]
  outputJSON=toJSON(output)
  write(outputJSON,file=file.path('output_html','webpage','attributeExplorer.json'))
  target <- HTMLInitFile(file.path('output_html','webpage'),
                         filename='attributeExplorer',
                         Title='User x Hash Tag Bispectral Co-Cluster Explorer')
  HTML(output,file=target,append=TRUE, Border=0, classtable='tablesorter',sortableDF=TRUE)
  HTMLEndFile()
}

merGe=function(nameKey,data){
  
  if(!(class(data$userID)=='numeric')) data$userID=as.numeric(data$userID)
  data=data[order(data$userID),]
  if(sum(nameKey$userID %in% data$userID)<dim(nameKey)[1]){
  newRows=nameKey$userID[!(nameKey$userID %in% data$userID)]
  
  m2=data.frame(newRows,matrix(c(rep(0,length(newRows)*(dim(data)[2]-1))),nrow=length(newRows),byrow=FALSE))
  colnames(m2)=colnames(data)
  m2=rBind(data,m2)
  m2=m2[order(m2$userID),]
  return(m2)}
  if(sum(nameKey$userID %in% data$userID)==dim(nameKey)[1]){return(data)}
}

htDataCleaner=function(k,htMat){
  temp=htMat[htMat$topic_cluster==k,]
  temp=temp[order(temp$degree,decreasing=TRUE),]
  d=dim(temp)[1]
  if(d>100) temp=temp[1:100,]
  output=data.frame(hashtag=temp$hashtag,degree=temp$degree,cluster=temp$topic_cluster,stringsAsFactors=FALSE)
}

coords2country = function(points)
{  
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # convert our list of points to a SpatialPoints object
  
  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  # return the ADMIN names of each country
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}

geoList=function(geo_path){
  
  geotags=fread(geo_path,data.table=FALSE,integer64='numeric')
  userIDs=unique(geotags$userID)
  countries=as.character(coords2country(cbind(geotags$lon,geotags$lat)))
  geotags=data.frame(userID=geotags$userID,country=countries,date=geotags$date,stringsAsFactors=FALSE)
  
  mostFrequent=function(userID,geotags){
    vec=geotags$country[geotags$userID==userID]
    vec=vec[!is.na(vec)]
    l=length(unique(vec))
    c=ifelse(length(vec)>1, names(which.max(table(vec))),vec)
    return(c(c,l))
  }
  
  
  pcountry=matrix(as.vector(sapply(userIDs,mostFrequent,geotags)),ncol=2,byrow=TRUE)
  pcountry=data.frame(userID=userIDs,country=pcountry[,1],ncountries=as.numeric(pcountry[,2]),stringsAsFactors=FALSE)
  
  return(list(pcountry,geotags))  
}

barPlotBuilder=function(userData,
                        targetIDs,
                        clusterType,
                        sort_type,
                        top_n){
  
  if(clusterType=='mention') clusterType='mention_cluster'
  if(clusterType=='hashtag') clusterType='topic_cluster'
  if(clusterType=='follow') clusterType='follow_cluster'
  if(clusterType=='ivcc') clusterType='ivcc_cluster'
  
  cSizeFunction=function(x,tempVec){
    s=sum(tempVec==x,na.rm=TRUE) 
    return(s)
  }
  
  tSizeFunction=function(x,temp,clusterType){
    
    s=sum(temp$targetIDs[temp[,clusterType]==x],na.rm=TRUE)
    return(s)
  }
  
  
  temp=userData
  
  cSeries=sort(unique(temp[,clusterType]))
  cSize=sapply(cSeries,cSizeFunction,temp[,clusterType])
  tSize=sapply(cSeries,tSizeFunction,temp,clusterType)
  ntSize=cSize-tSize
  barDF=data.frame(cluster=cSeries,seeds=tSize,non_seeds=ntSize,seed_proportion=tSize/cSize)
  barDF=barDF[order(barDF[,sort_type],decreasing=TRUE),]
  barDF=barDF[1:top_n,]
  barDF=barDF[order(barDF[,sort_type]),]
  barplot(height=as.matrix(t(barDF[,3:2])),horiz=TRUE,las=2,names.arg=barDF[,1],main=paste('Clusters of Interest\n cluster type: ',clusterType,'\n sort type: ',ifelse(sort_type=='seeds','seed count','seed proportion'),sep=''),ylab='cluster',xlab='size')
  
  
}