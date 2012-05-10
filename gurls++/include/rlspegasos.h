/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, IIT@MIT Lab
 * All rights reserved.
 *
 * authors:  M. Santoro
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


#ifndef _GURLS_RLSPEGASOS_H_
#define _GURLS_RLSPEGASOS_H_

#include "optimization.h"
#include "utils.h"

#include <set>

namespace gurls {

template <typename T>
class RLSPegasos: public Optimizer<T>{

public:
   void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};


template <typename T>
void RLSPegasos<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt)
{
   //	lambda = opt.singlelambda(opt.paramsel.lambdas);
   std::vector<double> ll = OptNumberList::dynacast(opt.getOpt("lambdas"))->getValue();
   T lambda = static_cast<T>((OptFunction::dynacast(opt.getOpt("singlelambda")))->getValue(ll.data(), ll.size()));

   gMat2D<T> X(X_OMR.cols(), X_OMR.rows());
   X_OMR.transpose(X);

   gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
   Y_OMR.transpose(Y);


//   [n,d] = size(X);
   const unsigned long n = X_OMR.rows();
   const unsigned long d = X_OMR.cols();

//   T = size(bY,2);
   const unsigned long t = Y_OMR.cols();


//   opt.cfr.W = zeros(d,T);
   gMat2D<T>* W = new gMat2D<T>(d,t);
   set(W->getData(), (T)0.0, d*t);
   opt.removeOpt("W");
   opt.addOpt("W", new OptMatrix<gMat2D<T> >(*W));

//   opt.cfr.W_sum = zeros(d,T);
   gMat2D<T>* W_sum = new gMat2D<T>(d,t);
   copy(W_sum->getData(), W->getData(), d*t);
   opt.removeOpt("W_sum");
   opt.addOpt("W_sum", new OptMatrix<gMat2D<T> >(*W_sum));

//   opt.cfr.count = 0;
   opt.addOpt("count", new OptNumber(0.0));


//   opt.cfr.acc_last = [];
//   opt.cfr.acc_avg = [];

//           opt.cfr.t0 = ceil(norm(X(1,:))/sqrt(opt.singlelambda(opt.paramsel.lambdas)));
   T* row = new T[d];
   getRow(X.getData(), n, d, 0, row);
   opt.addOpt("t0", new OptNumber( ceil( nrm2(d, row, 1)/sqrt(lambda))));

   delete[] row;


//   % Run mulitple epochs
//   for i = 1:opt.epochs,
   int epochs = static_cast<int>(opt.getOptAsNumber("epochs"));
   for(int i=0; i<epochs; ++i)
   {
//       if opt.cfr.count == 0
//           opt.cfr.t0 = ceil(norm(X(1,:))/sqrt(opt.singlelambda(opt.paramsel.lambdas)));
//           fprintf('\n\tt0 is set to : %f\n', opt.cfr.t0);
//       end

//       opt.cfr = rls_pegasos_singlepass(X, bY, opt);
       rls_pegasos_driver(X.getData(), Y.getData(), opt, n, d, Y_OMR.rows(), t);

   }

//   cfr = opt.cfr;

//   cfr.W = opt.cfr.W_sum/opt.cfr.count;

    T count = static_cast<T>(opt.getOptAsNumber("count"));
    if(eq(count, (T)0.0))
        throw gException(Exception_Illegat_Argument_Value);

    set(W->getData(), (T)0.0, W->getSize());
    axpy(W->getSize(), (T)(1.0/count), W_sum->getData(), 1, W->getData(), 1);

}

}
#endif // _GURLS_RLSPEGASOS_H_
