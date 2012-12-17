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
#include "optmatrix.h"

#include "float.h"

namespace gurls {

/**
 * \ingroup Performance
 * \brief PerfRmse is the sub-class of Performance that evaluates prediction error
 */

template <typename T>
class PerfRmse: public Performance<T>{

public:
    /**
     * Evaluates the root mean square error of the predicted labels stored in the field pred of opt with respect to the true input labels Y. It computes it as the frobenius norm over the classes and the samples of the difference between the true and predicted labels matrices.
     * \param X not used
     * \param Y labels matrix
     * \param opt options with the following:
     *  - pred (settable with the class Prediction and its subclasses)
     *
	 * \return perf, a GurslOptionList equal to the field pred of opt, with the following fields added or substituted:
     *  - rmse = root mean square error for each class/task
	 *  - forho = -rmse
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList* PerfRmse<T>::execute(const gMat2D<T>& /*X*/, const gMat2D<T>& Y, const GurlsOptionsList &opt) throw(gException)
{
    const unsigned long rows = Y.rows();
    const unsigned long cols = Y.cols();

    const T* y_true = Y.getData();

//    if isfield (opt,'perf')
//        p = opt.perf; % lets not overwrite existing performance measures.
//                  % unless they have the same name
//    end

    GurlsOptionsList* perf;

    if(opt.hasOpt("perf"))
    {
        GurlsOptionsList* tmp_opt = new GurlsOptionsList("tmp");
        tmp_opt->copyOpt("perf", opt);

        perf = GurlsOptionsList::dynacast(tmp_opt->getOpt("perf"));
        tmp_opt->removeOpt("perf", false);
        delete tmp_opt;

        perf->removeOpt("rmse");
        perf->removeOpt("forho");
//        perf->removeOpt("forplot");
    }
    else
        perf = new GurlsOptionsList("perf");


    const gMat2D<T> &pred = opt.getOptValue<OptMatrix<gMat2D<T> > >("pred");

    T *pred_t = new T[pred.getSize()];
    copy(pred_t, pred.getData(), pred.getSize());

    //pred = rows*cols

//     n 	= size(X,1);
    const T n = static_cast<T>(rows);

//     diff 	= opt.pred - y;
    axpy(rows*cols, (T)-1.0, y_true, 1, pred_t, 1);

//  p.rmse = norm(diff,'fro') / sqrt(n);
    T rmse = nrm2<T>(rows*cols, pred_t, 1)/sqrt(n);

    delete [] pred_t;

    perf->addOpt("rmse", new OptNumber(rmse));

//    p.forho 	= -p.rmse;
    perf->addOpt("forho", new OptNumber(-rmse));

//    p.forplot 	= p.rmse;
//    perf->addOpt("forplot", new OptNumber(rmse));

    return perf;
}

}

#endif //_GURLS_RMSE_H_
