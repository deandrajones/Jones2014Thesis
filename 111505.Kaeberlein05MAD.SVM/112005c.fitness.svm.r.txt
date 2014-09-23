library(e1071)
library(pcurve);

stein = read.csv("steinmetz.pca.tab", sep="\t", header=F );
row.names(stein) = stein[,1]

sub.s = stein[, seq(2,26)]
row.names(sub.s) = stein[,1]

tb = read.table("kaeberlein.short.tab" );
names(tb) = c("orf", "name", "arls", "ln2Ratio", "score");
k.orf = as.character( tb$orf );
row.names(tb) = tb$orf;

train = sub.s[ k.orf, ] # test the selection, good

###svm prediction
y = factor( tb$score )
x = train[,seq(1,10)]

model = svm( x, y);

     # test with train data
     pred <- predict(model, x)
     # (same as:)
#     pred <- fitted(model)

     # Check accuracy:
     table( pred, y)  #Ha, this is very poor. 


# d.gasch[ pred.gasch=="increase"] 


####################loose the different cutoff values

tb$score = "moderate";
tb$score[ (tb$arls/26.5)>1.07 ] = "increase";
tb$score[ (tb$arls/26.5)<0.9 ] = "decrease";
factor(tb$score);
y = factor( tb$score )
model = svm( x, y);
pred <- predict(model, x)
table( pred, y)  

 pred.stein = predict(model, sub.s[,seq(1,10)] );

summary(pred.stein);  #man this give 1917 increase genes


