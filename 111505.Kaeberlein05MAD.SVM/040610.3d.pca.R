#Ka, omega, fitness ~ RLS and svm 

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


#normalization
for( j in 2:505 ){
  my.mean = mean(scmd[,j],na.rm=T) 
  scmd[,j] = ( scmd[,j] - my.mean)/ sqrt( var( scmd[,j], na.rm=T)) #normalize by columns
  scmd[ is.na(scmd[,j]),j] = my.mean #Fill the missing dave with column-means
}

################
library(e1071)
library(pcurve);

scmdPCA = pca(scmd[,2:505])  #pca 
row.names(scmdPCA$pcs) = as.character( scmd[,1] )
scmdPCS = scmdPCA$pcs

### 3D scatter plot of the first 3 PC components
library(scatterplot3d);

PC1 = scmdPCS[,1]
PC2 = scmdPCS[,2]
PC3 = scmdPCS[,3]

 postscript("3Dscater.PC123.040610.ps", width=8,height=8); 
 s3d <- scatterplot3d( PC1, PC2, PC3, type='p',color="blue", pch=16,
  main="First Three Princial Components", 
  angle=45
  );  
 dev.off();
