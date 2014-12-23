#Ka, omega, fitness ~ RLS and svm 
setwd("~/projects/0.network.aging.prj/4.svm")

### read life span data
tb = read.csv("RLS.of.564.gene.deletion.BY.Managbanag08Plos.csv")
tb$ORF = as.character(tb$ORF)
row.names(tb) = as.character( tb$ORF )
tb$ratio = tb$RLS_Del_alpha / tb$BY4742 #R2 = 0.83

# partition lifespan into categories
cutoffs = quantile(tb$ratio, prob=c(0.20, 0.35, 0.65, 0.80), na.rm=T); 
tb$score = NA;
tb$score[tb$ratio< cutoffs[1] ] = "short"
tb$score[tb$ratio> cutoffs[4] ] = "long"
tb$score[tb$ratio<cutoffs[3]&tb$ratio>cutoffs[2] ] = "neutral"
table(tb$score)

### read Ka
Ktb = read.csv("Sce.Spa.KaKs.csv")
Ktb$orfname = as.character(Ktb$orfname)
row.names(Ktb) = Ktb$orfname

### read fitness
Ftb = read.csv( "growth.fitness.hom.csv")
Ftb$orf = as.character( Ftb$orf )

#### read morphology data
scmd = read.table( "mutant_analysis.tab", sep="\t", header=T)
row.names(scmd) = as.character( scmd$name )

#add Ka and Omega
scmd$Ka    = Ktb$Ka[ match( scmd$name, Ktb$orfname  ) ]
scmd$Omega = Ktb$Omega[ match( scmd$name, Ktb$orfname  ) ]

#add fitness
scmd$YPE = Ftb$YPE[ match( scmd$name, Ftb$orf) ]


#normalization all vectors
#for( j in 2:505 ){
for( j in 2:length(scmd[1,]) ){    
  my.mean = mean(scmd[,j],na.rm=T) 
  scmd[,j] = ( scmd[,j] - my.mean)/ sqrt( var( scmd[,j], na.rm=T)) #normalize by columns
  scmd[ is.na(scmd[,j]),j] = my.mean #Fill the missing value with column-means
}

################
library(e1071) #for svm
#library(pcurve); #for PCA, pcure disapeared in R 3.02
#scmdPCA = pca(scmd[,2:505])  #pca 
#row.names(scmdPCA$pcs) = as.character( scmd[,1] )
#scmdPCS = scmdPCA$pcs

# switch to princomp, 20141222
scmdPCA = princomp(scmd[,2:505], cor=TRUE)  #pca 
row.names(scmdPCA$scores) = as.character( scmd[,1] )
scmdPCS = scmdPCA$scores

p1 = scmdPCS[,1] #first pc
p2 = scmdPCS[,2] #first pc

######svm on morphology data alone, begin
shared.orfs = intersect(tb$ORF, row.names(scmd)) #shared orfs bw RLS and Gasch datasets
y = factor( tb[shared.orfs, "score"])
#x = scmdPCS[shared.orfs, ]
x = scmdPCS[shared.orfs, ]

model = svm(x,y, scale=F)
pred.train.scmd = predict( model, x)
table(pred.train.scmd, y)

pred.all.scmd = predict(model, scmdPCS); #run on all data
table( pred.all.scmd)

#find out novel predictions
known.long = tb$ORF[tb$score=="long"]
known.long = known.long[! is.na(known.long)]
pred.long = names( pred.all.scmd[ pred.all.scmd=="long"] )
shared = intersect(known.long, pred.long) #there are known long-genes
novel.long.orfs = setdiff(  pred.long, known.long ) #these are predict long-genes
######svm on morphology data alone, end

tmp = data.frame(novel.long.orfs)
write.table( tmp , "novel.long.orfs.verions2.tab", row.names=F, col.names=F, quote=F)


