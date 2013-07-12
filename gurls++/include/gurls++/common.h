/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-1013, IIT@MIT Lab
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


<p>GURLS++ is one of the 4 Grand Unified Laest Squares libraries (GURLS, GURLS++, bGURLS, bGURLS++),
a set of modular and easy-to-extend software libraries for efficient supervised learning,
which can be adopted effortlessly by both non-specialists and machine learning practitioners (http://lcsl.mit.edu/gurls.html).

<p>GURLS++, the C++ implementation of the Matlab toolbox GURLS, provides an object oriented framework to to solve supervised learning problems via
different state-of-the-art learning methods based on Regularized Least Squares.
It takes advantage of some favorable properties of regularized least squares algorithm.
is tailored to deal in particular with multi-category/multi-label problems.

<p>The library comprises useful routines to perform automatic parameter selection, and simple API's for specific
learning algorithms (currently only the Recursive RLS API is available but many others will be soon available).


The specification of the desired machine learning experiment in the library is very straightforward.
Basically, it is a formal description of a pipeline, i.e. an ordered sequence of steps.
Each step identifies an actual learning task, which can belong to a predefined category.
The core of the library is a class called GURLScore, which is responsible for processing the sequence of tasks
in the proper order and for linking the output of the former task to the input of the subsequent one.

\image html GURLS-design.png "GURLS++ pipeline"

A key role is played by the additional "options" structure, which we usually refer to as "opt".
It is used to store all configuration parameters required to customize the behaviour of single tasks in the pipeline.
Tasks receive configuration parameters from "opt" in read-only mode and - upon termination - the results are appended to the structure by GURLScore
in order to make them available to the subsequent tasks.
This allows the user to easily skip the execution of some tasks in a pipeline, by simply inserting the desired results directly into the options structure.
Currently, we identify six different task categories:
- automatic training
- validation
- test dataset splitting
- computation of the kernel matrix
- model selection
- optimization and training
- classifier prediction
- performance assessment.

All tasks in the same category can be interchanged with each other.

Further documentation:
- Installation instructions can be found here:
   https://github.com/CBCL/GURLS/wiki/2-Getting-Started#wiki-Installing_GURLSbGURLS
- Quick intructions on how to run the libraries for a default case can be found here:
   https://github.com/CBCL/GURLS/wiki/2-Getting-Started#wiki-Hello_World_in_GURLS-2
- A User manual with several examples can be found here:
   https://github.com/CBCL/GURLS/wiki/3-User-Manual#wiki-GURLS_Usage-2
- In gurls-manual.pdf you can find both the installation instructions and user manual,
  together with the Matlab and C++ Developer's Guide. GURLS is designed for easy expansion.
  Give it a try!


<p>The GURLS libraries have been developed within the context of the
IIT\@MIT Lab, established
between the Massachusetts Institute of Technology (MIT) and the
Istituto Italiano di Tecnologia (IIT) to develop novel learning and
perception technologies/algorithms for learning, especially in the
visual perception domain, that are inspired by the neuroscience of
sensory systems and are developed within the rapidly growing theory
of computational learning.

<p><strong>Current Development Status</strong></p>

The current version of this C++ library implements all the core
functionalities available in GURLS. Check out the main GURLS repository
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
* \defgroup Wrappers Wrappers
* \brief Contains gurls++ wrapper classes allowing users to easily bypass pipelines.
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
