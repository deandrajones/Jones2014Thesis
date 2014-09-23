library(e1071)
library(pcurves);

gasch = read.csv("gasch00.tab", sep="\t", header=T );
row.names(gasch) = gasch[,1]

sub.g = gasch[, seq(4,176)]
sub.g[ is.na(sub.g) ] = 0;
row.names(sub.g) = gasch[,1]
gpca = pca(sub.g );
row.names(gpca$pcs) = gasch[,1]

write.table(gpca$pcs, "gasch.pca.tab", col.name=F, sep="\t", quote=F);
gpcs = data.frame( gpca$pcs) ;

tb = read.table("kaeberlein.tab" );
names(tb) = c("orf", "name", "arls", "ln2Ratio", "score");
k.orf = as.character( tb$orf );

gpcs[ k.orf, 1:2] # test the selection, good


###svm prediction
y = factor( tb$score )
x = gpcs[k.orf, 1:43]

model = svm( x, y);

     # test with train data
     pred <- predict(model, x)
     # (same as:)
#     pred <- fitted(model)

     # Check accuracy:
     table( pred, y)  #Ha, this is 100% correct.

 pred.gasch = predict(model, gpcs[,1:43]);

d.gasch[ pred.gasch=="increase"] #damn, no predicted increase at all !!!


####################loose the different cutoff values

#decide to use 10% as cutoff 
tb$score = "moderate";
tb$score[ (tb$arls/26.5)>1.1 ] = "increase";
tb$score[ (tb$arls/26.5)<0.9 ] = "decrease";
factor(tb$score);
y = factor( tb$score )
model = svm( x, y);
pred <- predict(model, x)
table( pred, y)  #Ha, this is 100% correct.

 pred.gasch = predict(model, gpcs[,1:43]);

summary(pred.gasch); #shit 4573 decrease genes, only 5 increase genes from training set.


#### try again
tb$score = "moderate";
tb$score[ (tb$arls/26.5)>1.075 ] = "increase";  ##this give 26 increase from whole genome
tb$score[ (tb$arls/26.5)<0.85 ] = "decrease";
factor(tb$score);
y = factor( tb$score )
model = svm( x, y);
pred <- predict(model, x)
table( pred, y)  #Ha, this is 100% correct.

 pred.gasch = predict(model, gpcs[,1:43]);

summary(pred.gasch); 
n = row.names(gpcs);
plus.orf = n[ pred.gasch=="increase"]

gasch[plus.orf, 2]

write.table( tb, "kaeberlein.2.tab", quote=F, col.names=F, row.names=F,sep="\t");

write.table( gasch[plus.orf, 2], "111505.predicted.rls.extension.orfs.tab" , quote=F, row.name=F,col.names=F);



