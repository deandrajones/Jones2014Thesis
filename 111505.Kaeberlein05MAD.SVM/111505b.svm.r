library(e1071)
library(pcurves);

kg = read.table("k.and.gasch.tab", sep="\t" );
#names( kg[,c(1,5)] ) = c( "orf", "score" ) #does not work

#model = svm( score ~ ., data=tb);

x = kg[, seq(9,181)]

x[ is.na(x) ] = 0;
xpca = pca(x);

y = factor( kg[,5]);
#model = svm( xpca$pcs[,seq(1,10)], y); #this is poor
model = svm( xpca$pcs, y);

     # test with train data
#     pred <- predict(model, x)
     # (same as:)
     pred <- fitted(model)

     # Check accuracy:
     table(factor(pred), y)  #Ha, this is 100% correct.

#yy = as.character(y)
#predpred = as.character(pred)
#cbind( as.character(pred), as.character(y));
#table( as.character(pred), as.character(y));

