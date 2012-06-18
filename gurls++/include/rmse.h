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


#ifndef _GURLS_RMSE_H_
#define _GURLS_RMSE_H_

#include "perf.h"

#include "utils.h"
#include "gvec.h"

#include "float.h"

namespace gurls {

    /**
     * \brief Rmse is the sub-class of Performance that evaluates prediction error
     */

template <typename T>
class Rmse: public Performance<T>{

public:
    /**
     * Evaluates the root mean square error of the predicted labels stored in the field pred of opt with respect to the true input labels Y. It computes it as the frobenius norm over the classes and the samples of the difference between the true and predicted labels matrices.
     * \param X not used
     * \param Y labels matrix
     * \param opt options with the following:
     *  - pred (settable with the class Prediction and its subclasses)
     *
     * \return adds to opt the field perf wich is a list containing the following fields (if opt already contains perf than adds to it the following fields):
     *  - rmse = root mean square error
     */
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
void Rmse<T>::execute(const gMat2D<T>& X_OMR, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt) throw(gException)
{
    const int rows = Y_OMR.rows();
    const int cols = Y_OMR.cols();


    gMat2D<T> Y(cols, rows);
    Y_OMR.transpose(Y);
    T* y_true = Y.getData();

//    if isfield (opt,'perf')
//        p = opt.perf; % lets not overwrite existing performance measures.
//                  % unless they have the same name
//    end

    if(!opt.hasOpt("perf"))
    {
        GurlsOptionsList* perf = new GurlsOptionsList("perf");
        opt.addOpt("perf", perf);
    }
    GurlsOptionsList* perf = static_cast<GurlsOptionsList*>(opt.getOpt("perf"));


    GurlsOption *pred_opt = opt.getOpt("pred");
    gMat2D<T> *pred_mat = &(OptMatrix<gMat2D<T> >::dynacast(pred_opt))->getValue();

    gMat2D<T> pred_t(cols, rows);
    pred_mat->transpose(pred_t);

    //pred = rows*cols
//    T* diff = new T[rows*cols];

//     n 	= size(X,1);
    const T n = static_cast<T>(X_OMR.rows());

//     diff 	= opt.pred - y;
//    copy(diff, pred_t.getData(), rows*cols);
    axpy(rows*cols, (T)-1.0, y_true, 1, pred_t.getData(), 1);

//  p.rmse = norm(diff,'fro') / sqrt(n);
    T rmse = nrm2<T>(rows*cols, pred_t.getData(), 1)/sqrt(n);

    perf->addOpt("rmse", new OptNumber(rmse));

}

}

#endif //_GURLS_RMSE_H_
