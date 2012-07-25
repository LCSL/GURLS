%load the training data 
load('data/quickanddirty_traindata'); 
%train the classifier 
[opt] = gurls_train(Xtr,ytr); 
%load the test data 
load('data/quickanddirty_testdata'); 
%predict the labels for the test set and asses prediction accuracy 
[yhat,acc] = gurls_test(Xte,yte,opt); 

disp(sprintf('\nPrediction accurcay is:'))
disp(sprintf('\tClass %d\t',1:max(ytr)))
disp(sprintf('\t%2.2f%%\t',acc*100))
