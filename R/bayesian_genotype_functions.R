
CI<-function(X,alpha){
  PI<-X/sum(X)
  ConfidentI<-(-qnorm(alpha/8)*sqrt(PI*(1-PI)/sum(X))+1/2/sum(X))
  ConfidentI[X==0]<-(-log(alpha/4)/sum(X))
  return(ConfidentI)
}

get_probabilites <- function(X,alpha_dirichlet=c(0.5,0.5,0.5,0.5)*2,epsilon=0.01){
  ## Hypotheses 
  X<-sort(X,decreasing=TRUE)
  
  H1<-c(1,0,0,0)
  H2_1<-c(0.75,0.25,0,0)
  H2_2<-c(0.5,0.5,0.0,0)
  H2_3<-c(0.6666,0.33333,0.0,0)
  H3_1<-c(0.5,0.25,0.25,0)
  H3_2<-c(0.33333,0.333333,0.333333,0)
  H4_1<-c(0.25,0.25,0.25,0.25)
  
  E1<-ddirichlet((H1+epsilon)/sum(H1+epsilon),alpha_dirichlet+X)
  E2<-max(ddirichlet((H2_1+epsilon)/sum(H2_1+epsilon),alpha_dirichlet+X),
          ddirichlet((H2_2+epsilon)/sum(H2_2+epsilon),alpha_dirichlet+X),
          ddirichlet((H2_3+epsilon)/sum(H2_3+epsilon),alpha_dirichlet+X))
  E3<-max(ddirichlet((H3_1+epsilon)/sum(H3_1+epsilon),alpha_dirichlet+X),
          ddirichlet((H3_2+epsilon)/sum(H3_2+epsilon),alpha_dirichlet+X))
  E4<-ddirichlet((H4_1+epsilon)/sum(H4_1+epsilon),alpha_dirichlet+X)
  
  
  
  while(sort(c(E1,E2,E3,E4),decreasing=TRUE)[2] == 0 ){
    
    X <- X/10
    E1<-ddirichlet((H1+epsilon)/sum(H1+epsilon),alpha_dirichlet+X)
    E2<-max(ddirichlet((H2_1+epsilon)/sum(H2_1+epsilon),alpha_dirichlet+X),
            ddirichlet((H2_2+epsilon)/sum(H2_2+epsilon),alpha_dirichlet+X),
            ddirichlet((H2_3+epsilon)/sum(H2_3+epsilon),alpha_dirichlet+X))
    E3<-max(ddirichlet((H3_1+epsilon)/sum(H3_1+epsilon),alpha_dirichlet+X),
            ddirichlet((H3_2+epsilon)/sum(H3_2+epsilon),alpha_dirichlet+X))
    E4<-ddirichlet((H4_1+epsilon)/sum(H4_1+epsilon),alpha_dirichlet+X)
    
  }
  return(log10(c(E1,E2,E3,E4)))
}

calcK <- function(X,islog=TRUE,Kthresh=2){

  ## Alternative
  Res <- rep(0,4)
  Res[which.max(X)] <- 1  
  ########################
  #return(c(H,D,Tr,Q))   
  return(Res)  
  
}

getMutation <- function(db,germline_db,sequnece_gap='SEQUENCE_IMGT',vcall='V_CALL'){
  
  mut <- sapply(1:nrow(db),function(x){mutCount(read = db[x,paste(sequnece_gap)],
                                                ref=germline_db[which(names(germline_db)==strsplit(x = db[x,paste(vcall)],',')[[1]][1])])})
  return(mut)
}


getGenesWithMutataions <- function(db, mut, vcall='V_CALL', mutNum=0 ,minReads = 1){
  
  ### Only for genes that have one V assignmnet
  IND.ONE <- grep(",",invert = T,value = F,x = db[,paste(vcall)])
  db.NEW <- db[IND.ONE,]
  mut.list <- table(db.NEW[mut <= mutNum,paste(vcall)])
  
  ### ADD partial counts for sequences with multiple assignmnets
  IND.MORE <- grep(",",invert = F,value = F,x = db[,paste(vcall)])
  VCALLS.MORE <- db[IND.MORE,paste(vcall)]
  VCALLS.MORE <- strsplit(VCALLS.MORE,',')
  
  sapply(VCALLS.MORE, function(x){m <- 1/length(x); 
  sapply(x,function(y){if(length(mut.list[names(mut.list)==y])==0){mut.list <<- c(mut.list,m);names(mut.list)[length(mut.list)] <<- paste(y)} 
    else { mut.list[names(mut.list)==y] <<- mut.list[names(mut.list)==y]+m}})})
  
  genes <- unique(sapply(strsplit(names(mut.list),'*',fixed=T),'[',1))
  
  alleles.db  <- list()
  
  sapply(genes,function(x){alleles.db[[x]] <<- mut.list[grep(paste0(x,'*'),names(mut.list),fixed=T)]}) 
  alleles.db <- sapply(alleles.db,sort,decreasing=T)
  
  sapply(1:length(alleles.db),function(x){names(alleles.db[[x]]) <<- sapply(strsplit(names(alleles.db[[x]]),"*",fixed=T),'[',2)})
  
  alleles.db[sapply(alleles.db,sum) < minReads] <- NULL
  alleles.db <- lapply(alleles.db,sort,decreasing=T)
  
  #   if(length(alleles.db) == 1) {
  #     alleles.table <- data.frame(GENE = colnames(alleles.db), ALELLE = paste(rownames(alleles.db),collapse = ','),
  #                                 COUNT= paste((alleles.db),collapse = ',') )
  #   } else {
  alleles.table <- data.frame(GENE = names(alleles.db), ALELLE = sapply(alleles.db,function(x){paste(names(x),collapse = ',')}),
                              COUNT=sapply(alleles.db,function(x){paste((x),collapse = ',')}) )
  #}
  
  
  
  return(list(genes = genes,allels_db=alleles.db,alleles_table=alleles.table))
  
}

mutCount <- function(read,ref){
  mutnum <- 0
  
  N <- min(nchar(ref),nchar(read))
  ref <- toupper(substr(ref,1,N))
  read <- toupper(substr(read,1,N))
  
  ref.split <-  unlist(strsplit(read,''))
  last.dot.idx <- ifelse(ref.split[1] == '.',which(ref.split != '.')[1]-1,0)
  
  mutnum <- sum(unlist(sapply(1:nchar(read),function(x){if((substr(read,x,x) != substr(ref,x,x)) 
                                                           && x >  last.dot.idx
                                                           && substr(read,x,x) != "N" )
  {mutnum=mutnum+1}   })))
  return(mutnum)
  
}


extractForPlot <- function(tab){
  IND <- grep('KDiff',names(tab))
  Kdiff <- as.numeric(as.character(tab[[IND]]))
  Ks <- sapply(Kdiff,function(x){ if(0 <= x && x <= log10(3)) return('none');
    if(log10(3) < x && x <= log10(20)) return('positive');
    if(log10(20) < x &&  x <= log10(150)) return('strong');
    if(log10(150) < x) return('very_strong')})
  IND_A <- grep('ALELLE',names(tab))
  IND_D <- grep('Decisions',names(tab))
  tab[,IND_D] <- sapply(strsplit(as.character(tab[,IND_D]),'_'),'[',1)
  alleles <- apply(tab,1,function(x){allele <- unlist(strsplit(as.character(x[IND_A]),','));
  if(x[IND_D] == 'H') return(paste(allele[1]));
  if(x[IND_D] == 'D') return(paste(allele[1:2],collapse = ','));
  if(x[IND_D] == 'T') return(paste(allele[1:3],collapse = ','));  
  if(x[IND_D] == 'Q') return(paste(allele[1:4],collapse = ',')); 
  })
  IND_C <-  grep('COUNT',names(tab))
  al1 <- round(as.numeric(sapply(strsplit(as.character(tab[,IND_C]),","),"[",1)),digits = 2)
  al2 <- round(as.numeric(sapply(strsplit(as.character(tab[,IND_C]),","),"[",2)),digits = 2)
  al3 <- round(as.numeric(sapply(strsplit(as.character(tab[,IND_C]),","),"[",3)),digits = 2)
  al4 <- round(as.numeric(sapply(strsplit(as.character(tab[,IND_C]),","),"[",4)),digits = 2)
  
  raw.count <- sapply(paste(al1,al2,al3,al4, sep=','), function(x){sapply(gregexpr(',NA',x,fixed = T),function(y){if(y[1] != -1){substr(x,1,y[1]-1)} else (x)})})
  
  return(data.frame(GENE=tab[,grep('GENE',names(tab))],RAW_COUNT=raw.count,ALLELES=alleles,EVIDENCE=Ks,K=Kdiff))
}

# PARAMETERS
# db - clip_db
# germline_db 
# thresh - max number of mutations allowed to count the sequnce for genotyping
# vcall - name of the column of the V assignmnet
# alpha_dirichlet and epsilon - Dirichlet distribution parametes - not to be changed for now
getGenotype <- function(db,germline_db,thresh=3,vcall='V_CALL',
                        alpha_dirichlet=c(0.5,0.5,0.5,0.5)*2, epsilon=0.01){
  
  all_genes = germline_db %>% names %>% getGene(strip_d=FALSE) %>% unique %>%
    data.frame(GENE = .)
  
  ## number of mutations for each sequence relative to germline assigned
  mut = getMutCount(db$SEQUENCE_IMGT, db$V_CALL, germline_db) %>% 
    sapply(function(x) min(unlist(x), na.rm=T))
  
  genotype <- c()

  # minReads - the minimal number of reads for a gene to consider for genotyping
  # mutNum -  max number of mutations allowed to count the sequnce for genotyping
  genes <- getGenesWithMutataions(db = db,mut = mut,vcall = vcall, mutNum = thresh,minReads = 5)
    
  alleles.db  <- genes$allels_db
  alleles.table <- genes$alleles_table
    
  if(nrow(alleles.table)!=0){  
    probs <- t(sapply(alleles.db,function(x){len=min(length(x),4);get_probabilites(sort(c(x,rep(0,4-len)),decreasing = T)[1:4])}))
    probs[probs==-Inf] <- -1000
    colnames(probs) <- c('H','D','T','Q')
      
    find_H <- t(apply(probs,1,function(x)calcK(x,islog=TRUE,Kthresh=0)))
    colnames(find_H) <- c('H','D','T','Q')
      
    tmp.table <- cbind(apply(find_H,1,function(x){paste(colnames(find_H)[(x==max(x)) & x > 0],collapse ="_")}),probs)
    colnames(tmp.table)[2:5] <- c(paste0(colnames(tmp.table)[2:5]," (log)"))
    colnames(tmp.table)[1] <- paste0('Decisions')
    if(nrow(tmp.table)==1) { 
      tmp.table <- cbind(tmp.table,apply(matrix(tmp.table[,2:5],nrow=nrow(tmp.table)),1,function(x){k <- sort(as.numeric(x),decreasing = T);k[1]-k[2]}))
    } else {
      tmp.table <- cbind(tmp.table,apply(tmp.table[,2:5],1,function(x){k <- sort(as.numeric(x),decreasing = T);k[1]-k[2]})) 
    }
      
      colnames(tmp.table)[6] <- paste0('KDiff')
      
      names(alleles.table)[2:3] <- paste0(names(alleles.table)[2:3])
      alleles.table <- cbind(alleles.table,tmp.table)
      genotype <- cbind(extractForPlot(alleles.table))
      row.names(genotype) <- NULL
    }
    
return(genotype)
}