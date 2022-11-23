##### runfactoization runs all the considered multi-omics factorization
### the required inputs are:
### "folder" corresponding to the path to the folder where the input files are contained, in the idea that all omics matrices are organized inside a unique folder
### "file.names" corresponding to a vector containing the names of all the omics files
### "num.factors" containing the number of factors in which we desire to decompose the matrices
### "sep=" "" corresponding to the separator used in the omics files required to properly read them
### "single.cell" indicating if the data are single cell data. In this case the filtering of the data will be more intense in respect to other data types
### "filtering"

### the input files need to be log2 transformed before running the analysis

library("RGCCA")
library("r.jive")
library("IntNMF")
library("omicade4")
library("GPArotation")
source("scripts/tICA.R")
library("iCluster")

# skipping these methods for now, skipped by authors
# deleted the code that runs these methods
#library("MSFA")

# skipping these methods for now, difficult to install
# deleted the code that runs these methods
#library("MOFAtools")
#library("MOFA")
#library("MOFA2")
#library("tensorBSS")

runfactorization <- function(folder, file.names, num.factors, sep=" ", filtering="none"){
  factorizations<-list()
  t<-1
  method<-numeric(0)
  num.factors<-as.numeric(num.factors)
  
  ##creating list of omics
  omics <- list()
  for(i in 1:length(file.names)){
    fn = paste(folder, file.names[i], sep="/")
    print(fn)
    omics[[i]] <- read.table(fn, sep=sep, row.names=1, header=T)
    omics[[i]] <- as.matrix(omics[[i]])
  }
  
  ####
  #omics<-lapply(omics, function(x) t(x))
  ######
  
  ##restricting to common samples and filtering
  samples <- colnames(omics[[1]])
  print(paste0('samples initially: ', samples))
  print(paste0('num samples initially: ', length(samples)))
  for(j in 1:length(omics)){
    samples <- intersect(samples,colnames(omics[[j]]))
  }
  print(paste0('samples post-intersect: ', samples))
  print(paste0('num samples post-intersect: ', length(samples)))



  for(j in 1:length(omics)){
      omics[[j]]<-omics[[j]][,samples]

      print(head(omics[[j]]))

      if(filtering!="none"){
          x<-apply( omics[[j]],1,sd)
          x<-as.matrix(sort(x, decreasing = T))
          w<-which(x>0)
          if(filtering=="stringent"){
              selected <- rownames(x)[1:min(w[length(w)],5000)]
          }
          else{
              selected <- rownames(x)[1:min(w[length(w)],6000)]
          }
          m <- match(rownames(omics[[j]]),selected)
          w <- which(!is.na(m))
          omics[[j]] <- omics[[j]][w,]
      }
      else{
          #omics[[j]] <- omics[[j]][which(apply(omics[[j]],2,sd)>0),]
          omics[[j]] <- omics[[j]]
      }

      fn = paste0(folder, '/', file.names[j], '.harmonized.tsv')
      write.table(omics[[j]], fn, sep='\t')

      print(paste('Filtering on rows where the std is low', dim(omics[[j]])))
    }  


  ##### JIVE
  print("##### JIVE #####")
  factorizations_jive <- jive(omics,
                              rankJ=num.factors,
                              rankA = rep(num.factors, length(omics)),
                              method = "given",
                              conv = "default",
                              maxiter = 100,
                              showProgress=FALSE)

  rankJV <- factorizations_jive$rankJ;
  rankIV.v <- factorizations_jive$rankA;

  J<-numeric(0)
  ng<-0
  metagenes_jive <- list();
  for(j in 1:length(omics)){
    J <- rbind(J,factorizations_jive$joint[[j]]);
    ng<-c(ng,dim(factorizations_jive$joint[[j]])[1])
  }

  svd.o <- svd(J);
  jV <- svd.o$v %*% diag(svd.o$d);

  for(j in 1:length(omics)){
    metagenes_jive[[j]] <- svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]; ###error in dimension
    rownames(metagenes_jive[[j]])<-rownames(omics[[j]])
    colnames(metagenes_jive[[j]])<-1:num.factors
  }

  factors_jive=jV[,1:rankJV]
  rownames(factors_jive)<-colnames(omics[[1]])
  colnames(factors_jive)<-1:num.factors
  factorizations[[t]]<-list(factors_jive,metagenes_jive)
  t<-t+1
  method<-c(method,"JIVE")
  
  out<-list(factorizations=factorizations, 
            method=method)
  
  return(out)
}
