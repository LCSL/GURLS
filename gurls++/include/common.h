 /*
  * The GURLS Package in C++
  *
  * Copyright (C) 2011, IIT@MIT Lab
  * All rights reserved.
  *
  * authors: M. Santoro
  * email:   msantoro@mit.edu
  * website: http://cbcl.mit.edu/IIT@MIT/IIT@MIT.html
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions
  * are met:
  *
  *     * Redistributions of source code must retain the above
  *       copyright notice, this list of conditions and the following
  *       disclaimer.
  *     * Redistributions in binary form must reproduce the above
  *       copyright notice, this list of conditions and the following
  *       disclaimer in the documentation and/or other materials
  *       provided with the distribution.
  *     * Neither the name(s) of the copyright holders nor the names
  *       of its contributors or of the Massacusetts Institute of
  *       Technology or of the Italian Institute of Technology may be
  *       used to endorse or promote products derived from this software
  *       without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGE.
  */


#ifndef _GURLS_COMMON_H_
#define _GURLS_COMMON_H_

/*!
\mainpage The GURLS Package in C++

<p>This C++ library provides an object oriented framework to solve
supervised mult-class learning problems using a regularized least
square approach.

<p>GURLS++ is the C++ implementation of the Matlab toolbox GURLS
(http://cbcl.mit.edu/gurls/).

<p>The ultimate goal of the authors and contributors to GURLS++ has
been to develop a comprehensive library to solve easily classication
and regression problems both in the setting of batch and online
learning from examples.

<p>Both GURLS and GURLS++ allows for training multiclass classifiers
based on the Regularized Least Squares (RLS) loss function.
The libraries encode ideas and experience from
classical and more recent developments in theory and practice of
multicategory classification. The software package includes efficient methods
to automatically select good model parameters, to use a variety of
output coding schemes and an easy interface for experiment management
and results visualization. The software package includes a number of
routines especially designed to deal with massive datasets,
furthermore, it is instrinsically modular and each method can be
used independently.


<img src="../pipeline.png"  title="" alt="pipeline" height="264" width="480" />

<p>GURLS and GURLS++ have been developed within the context of the
IIT\@MIT Lab, established
between the Massachusetts Institute of Technology (MIT) and the
Istituto Italiano di Tecnologia (IIT) to develop novel learning and
perception technologies/algorithms for learning, especially in the
visual perception domain, that are inspired by the neuroscience of
sensory systems and are developed within the rapidly growing theory
of computational learning.

<p><strong>Software Package</strong></p>

The GURLS package is more than a collection of optimization
routines, its design allows for quickly composing various experiments
by varying the parameter selection, performance evaluation and the
optimization procedures. The GURLS package has three main features:

<ol>
<li>A core library consisting of parameter selection and a wide range of optimization routines for building regularized least squares classifiers.
<li>Automated tools for building principled classifiers, where most decisions for building an optimal classifier are taken by the package based on theory.
<li>A simple and easy-to-use, and easy-to-extend benchmarking framework for rapidly prototyping learning algorithms with various components.
</ol>

<p><strong>Benchmarking Framework</strong></p>

The GURLS package allows the user to quickly design and compare
experiments, run with different pattern recognition pipelines
(for example, comparing the performance of linear vs. non-linear
classifiers on a given dataset). The output of each task in the
sequence is stored in a particular format by GURLS, which makes it
easy for the package to combine the results and present them as plots
or tables.


<p><strong>Current Development Status</strong></p>

The current version of this C++ library implements all the core
functionalities available in GURLS but still miss the functionalities
required to handle large scale machine learning. A major upgrade and
revision will be available soon. Check out the main GURLS repository
at (https://github.com/CBCL/GURLS) and keep updated for futher developments.


<p><strong>Author</strong>:</p>

<p><em>Matteo Santoro</em>. (<em>Contacts:</em> matteo.santoro@gmail.com)</p>


<p><strong>Copyright Notice</strong>:</p>


<p>Copyright (C) 2013, Matteo Santoro and IIT\@MIT Lab</p>

<p>All rights reserved.</p>


GURLS++ is released under the following licence:


  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions
  * are met:
  *
  *     * Redistributions of source code must retain the above
  *       copyright notice, this list of conditions and the following
  *       disclaimer.

  *     * Redistributions in binary form must reproduce the above
  *       copyright notice, this list of conditions and the following
  *       disclaimer in the documentation and/or other materials
  *       provided with the distribution.

  *     * Neither the name(s) of the copyright holders nor the names
  *       of its contributors or of the Massacusetts Institute of
  *       Technology or of the Italian Institute of Technology may be
  *       used to endorse or promote products derived from this software
  *       without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGE.


  <p><strong>History of revisions</strong></p>

<ol>


<li> Version 2.0.00, released 2013 January 24th.

<li> Version 1.0.00, released 2012 June 1th. GURLS++ 1.0.00 is the first
complete release.

<li> Version 0.1.00, released 2011 October 7th. GURLS++ 0.1.00 is an initial
release and should be used with caution. The library has already been
in use, but not every function/method has been thoroughly tested.
Only a limited subset of GURLS is interfaced. Commens, bugfixes,
enhancements, suggestions are wellcome.

</ol>
*/

/**
 * \defgroup LinearAlgebra LinearAlgebra
 * \brief Contains classes representing 1D/2D matrices and methods to compute basic linear algebra operations on them.
 */

/**
 * \defgroup Common Common
 * \brief Contains common constants, typedes, classes and methods used within the library.
 */

/**
 * \defgroup Kernels Kernels
 * \brief Contains classes that compute the kernel matrix.
 */

/**
 * \defgroup Optimization Optimization
 * \brief Contains classes and methods to implement the abstract concept of an <em>optimization algorithm</em>.
 */

/**
 * \defgroup ParameterSelection ParameterSelection
 * \brief Contains classes and methods to implement parameter selection.
 */

/**
 * \defgroup Prediction Prediction
 * \brief Contains classes that compute predictions.
 */

/**
 * \defgroup Performance Performance
 * \brief Contains classes and methods to evaluate predictions performances.
 */

/**
 * \defgroup Norms Norms
 * \brief Contains classes that spherify the data.
 */

/**
 * \defgroup Split Split
 * \brief Contains classes that split data into pairs of training and test samples.
 */

/**
 * \defgroup PredKernels PredKernels
 * \brief Contains classes that computes the kernel matrix for prediction.
 */

/**
 * \defgroup Confidence Confidence
 * \brief Contains classes that compute a confidence score for the predicted labels.
 */

/**
 * \defgroup Settings Settings
 * \brief Contains classes and methods useful to define a set of parameters used by any learning algorith.
 */

/**
 * \defgroup Tutorials Tutorials
 * \brief Contains example files that show how to use the GURLS++ API to build simple machine learning applications.
 */

/**
 * \defgroup Exceptions Exceptions
 * \brief Contains classes representing errors.
 */

/**
 * \ingroup Common Common
 * \file
 * \brief Header file containing a number of common definition.
 */


//! The main namespace of GURLS++
/*!
  \namespace gurls
  The namespace gurls comprises all the classes, constants, typedefs,
  and functions  currently available within the library.
  */
namespace gurls {

}


#endif // _GURLS_COMMON_H_
