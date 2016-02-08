# GURLS - Grand Unified Regularized Least Squares



## Table of Contents

- Introduction
- Installation
- Quickstart
- Documentation
- Further Informations



## Introduction

GURLS - (Grand Unified Regularized Least Squares) is a MATLAB software package for regression and 
(multiclass) classifiers based on the Regularized Least Squares (RLS) loss function. 

## Installation

Open MATLAB and execute:

```matlab
>> run('PACKAGEROOT/gurls/utils/gurls_install.m');
```

where `PACKAGEROOT` is where you unzipped the GURLS package. This will add all the important directories 
to your path. Run 'savepath' if you want the installation to be permanent.




## Quick start

A straightforward example of GURLS usage can be found in gurls_helloworld.m, 
which employs gurls_train and gurls_test functions. 
The most important functionalities of the library can be controlled just
by using gurls_train and gurls_test with the appropriate options.

Below we describe more in details how to use them on your data.

Put your training data in a n-by-d matrix X with each row being a data sample.
Put the training labels in a y n-by-one vector. Let the labels go from 1 to T with T 
being the number of classes.

Run:

```matlab
>> model = gurls_train(X, y);
```

This will train a non-linear model based on a gaussian kernel and validated by hold out cross validation.

Put your test data in a n-by-d matrix Xte like above and your test labels in 
a n-by-one vector yte. Let the labels go from 1 to T with T being the number of classes.

Run: 

```matlab
>> [ypred, acc] = gurls_test(model, Xte, yte)
```

This will compute the outputs predicted by the model with the associated accuracy. 


## Documentation

The functions `gurls_train` and `gurls_test` solve different kinds of supervised learning problems: regression, classification, multiclass classification. The functions automatically recognize the kind of problem, find the most suitable model and learning algorithm for the given dataset and select the appropriate metaparameters, if not specified in the options.

### `gurls_train`
The function train has different signatures:
```matlab
model = gurls_train(Xtrain, ytrain);
model = gurls_train(Xtrain, ytrain, options);
model = gurls_train(Xtrain, ytrain, 'optname1', optval1, 'optname2', optvalue2, ...);
```
`Xtrain` can be a n x n matrix of the form `[x1; x2;... xn]` where `xi` are n training inputs of d dimensions. Otherwise it can be a precomputed Gram matrix (thus a n x n positive semidefinite matrix).

`ytrain` is a n x T matrix where the i-th row is the vector of labels associated to the i-th sample point `xi` in `Xtrain`. `ytrain` is a matrix of real numbers in the regression tasks, while a matrix of -1,+1 in the classification or multiclassification ones, where the `ytrain(i,j)` is 1 if `xi` belongs to the class j otherwise it is -1. This notation is general and covers the case of overlapping classes. Note that for the multiclass case another allowed format for `ytrain` is a n x 1 vectors of numbers in {1,2,..,T}, where `ytrain(i)=j` if `xi` belongs to the class j.

#### Training `options`

`options` is a list of options and is defined as follows
```matlab
options = struct('optname1', optvalue1, 'optname2', optvale2, ...)
```

if `options` is defined in this way the second and the third signatures are equivalent.

`'optname1', optval1,...` these are options with the associated values that are passed to the train function. Note that all the options specified here will be passed to the underlying modules of gurls. The following is a list of options with the possible values.


*  `datatype`: (automatically deduced, if not specified) can be `'vector'` when the input type is a matrix with the format specified for `Xtrain`, or `'kernel'` when ` Xtrain` is a Gram matrix. 
* `'problem'`: (automatically deduced, if not specified) can be `'regression'` or `'classification'`, see the definition of `'ytrain'`
* `'algorithm'`: it specifies the algorithm to use for solving the problem, it can be 

    * `'lrls'`: (Only when the datatype is vector) it is linear regularized least squares,
    * `'krls'`: (default) it is kernel regularized least squares,
    * `'krlsrf'`: it is kernel regularized least squares with random features,
    * `'gpr'`: it is gaussian processes regression.


* `'filter'`: (automatically deduced, if not specified) it specifies the filter used by the algorithm, it can be 

    * `'tikh'`: (default) tikhonov filter,
    * `'land'`: landweber iterative filter. To use this, specify the parameter ` regrange` with the range of metaparameters to test (e.g. ` 1:100`),
    * `'nu'`: nu-method iterative filter. To use this, specify the parameter ` regrange` with the range of metaparameters to test (e.g. ` 1:100`),
    * `'conjgrad'`: conjugate gradient filter,
    * `'randtikh'`: randomized tikhonov filter.


* `'kernelfun'`: it specifies the kernel function used by the algorithm. It applies for all the algorithms except for lrls that is linear by default. Possible values are

    * `'linear'`: linear kernel,
    * `'rbf'`: (default) gaussian kernel,
    * `'quasiperiodic'`: quasiperiodic kernel, (the options ` paramsel.alpha` and ` period` must be specified).
    * `'chisquared'`: chi-square kernel,
    * `'datatype'`: used when `'Xtrain'` is a precomputed Gram matrix.


* `'partuning'`:  specifies the cross validation approach to be used for automatically selecting the metaparameters. It can be `'loo'` for leave one out, or `'ho'` for hold-out. The default is ` ho`. In case the user wants to provide a separate test set, this can instead be provided as below
* `'Xva'` and `'yval'`: (both need to be provided) user-provided input and output validation sets. They need to be matrices with arbitrary number of rows but respectively same number of columns to `X` and `y`.
* `'perfm'`: measure of the performances of the model. It can be

    * `'rmse'`: root mean square error (default when the problem is of regression type),
    * `'macroavg'`: macro average (default when the problem is of classification type),
    * `'precrec'`: precision recall,
    * `'gpregr'`: specific performance measure for the gaussian process regression (mandatory when the algorithm is gpr).
    
    
* `'pars'`: (automatically deduced, if not specified) it tells the system which metaparameters must be selected.

    * `'none'`: both the regularization and the kernel metaparameters are specified in the options,
    * `'reg'`: the regularization metaparameter has to be found by cross validation,
    * `'ker'`: the kernel metaparameter has to be found by cross validation,
    * `'all'`: both the regularization and the kernel metaparameter has to be found by cross validation.

* `'regrange'`: the range where the automatic cross validation will search for the best regularization parameter. It must be specified for the nu and land filters, otherwise it is automatically computed by the system.
* `'kerrange'`: the range where the automatic cross validation will search for the best kernel parameter. It is automatically computed by the system, if not specified.
* `'regpar'`: specific value for the regularization metaparameter.
* `'kerpar'`: specific value for the kernel metaparameter.


### `gurls_test`

Available signatures 

```matlab
ypredicted = gurls_test(model, Xtest)
[ypredicted, accuracy] = gurls_test(model, Xtest, ytest)
[ypredicted, accuracy] = gurls_test(model, Xtest, ytest, perfmeas)
```

The first signature computes the predicted labels `ypredicted` associated to the test dataset. The second computes the predicted labels and their accuracy with respect to the test label `ytest`. In this case the accuracy is computed by the performance measure used in the training phase. In order to compute it with a different performance measure, use the third signature. Note that the format of `Xtest` and of `ytest` must be the same of `Xtrain` and `ytrain`.

### `GurlsOption`
The `GurlsOptions` is an internal class intended for collecting all the variables needed by the function `gurls`. It is the class of `model` that is the model learned by `gurls_train` and is the class of `opt` that is the object automatically produced by `gurls_defopt`. Let `namevar` refers to a variable of interest, use 

* `isprop(opt, 'namevar')` to check if it is in `opt`  
* `obj.namevar` to read its value 
* `obj.namevar = value` to write its value (when the variable is in `obj`)
*  `opt.newprop('namevar', value)` to add it to `opt`. If the variable is already present, it will be overwritten.


Moreover you can 

*  add multiple variables by using `opt.newprops(struct('name1',value1,'name2', value2,...))`. The variables in opt with the same name of the ones in oldOpt will be overwritten.
* copy the content of a existent opt in a new one by using`opt.newprops(oldOpt)`. The variables in opt with the same name of the ones in oldOpt will be overwritten.
* create subvariables with `opt.newprop('name1.name2',value)`. It will be accessed by `opt.name1.name2`. Note that now `opt.name1` contains `struct('name2', value)`.
Thus the code

```matlab
opt.newprop('name1.name2',value1);
opt.newprop('name1.name3',value2);
```

is equivalent to write

```matlab
opt.newprop('name1', struct());
opt.name1.name2 = value1;
opt.name1.name3 = value2;
```
Note that the object `opt` is of class `GurlsOptions`, that is a subclass of `handle`, thus it is passed by reference and not by value as shown in the following example

```matlab
opt = defopt('foobar');
opt.newprop('v1', 5);
fun(opt);
disp(opt.v1); % it will display 6

function fun(opt)
  if isprop(opt, 'v1')
    opt.v1 = opt.v1 + 1;
  end
end
```



## Further Informations

- GURLS design is described here:
	https://github.com/CBCL/GURLS/wiki/3-User-Manual#wiki-Design

- Demos
	GURLS has exstensively commented demos in the "demo" subdirectory. 
	Have a look, and run gurls_helloworld.m. We feel this is the best way to learn how to use
	these packages.

- A Quick reference with several examples can be found here:
	https://github.com/LCSL/GURLS/blob/master/gurls/guide-train-test.pdf?raw=true

- A User manual can be found here:
	https://github.com/CBCL/GURLS/wiki/3-User-Manual#wiki-GURLS_Usage

- A collection of the most useful and common pipelines can be found here:
	https://github.com/CBCL/GURLS/wiki/3-User-Manual#wiki-Examples_in_GURLS
	
- Description of the available methods, demos and data for each package can be found at
	https://github.com/CBCL/GURLS/wiki/4-Code-Description#wiki-GURLS

- Developer's Guide
	A simple developer's guides is available in the gurls-manual.pdf file.
	GURLS is designed for easy expansion. Give it a try!




