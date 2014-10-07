%load the training data 
load('data/quickanddirty_traindata'); 
%train the classifier 
model = train(Xtr,ytr, 'nsigma', 10, 'nlambda', 10); 
%load the test data 
load('data/quickanddirty_testdata'); 
%predict the labels for the test set and asses prediction accuracy 
[yhat,acc] = test(model, Xte,yte); 

disp(sprintf('\nPrediction accuracy is:'))
disp(sprintf('\tClass %d\t',1:max(ytr)))
disp(sprintf('\t%2.2f%%\t',acc*100))
