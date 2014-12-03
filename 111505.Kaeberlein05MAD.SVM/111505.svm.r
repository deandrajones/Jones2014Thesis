library(e1071)

tb = read.table("kaeberlein.tab" );
names(tb) = c("orf", "name", "arls", "ln2Ratio", "score");

hist(tb$arls/26.5 )

#decide to use 20% as cutoff 
tb$score = "moderate";
tb$score[ (tb$arls/26.5)>1.2 ] = "increase";
tb$score[ (tb$arls/26.5)<0.8 ] = "decrease";
factor(tb$score);

write.table( tb, "kaeberlein.tab", row.name=F, col.name=F, quote=F, sep="\t");

attach(tb);

#model = svm( score ~ ., data=tb);
x = subset( tb, select = -c(orf,name,score) )
y = factor( tb$score );
model = svm( x, y);

     # test with train data
     pred <- predict(model, x)
     # (same as:)
     pred <- fitted(model)

     # Check accuracy:
     table(pred, y)  #this is 100% correct as expected. 
