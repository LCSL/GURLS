/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-1013, IIT@MIT Lab
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

#include "gurls++/perf.h"

#include "gurls++/utils.h"
#include "gurls++/gvec.h"
#include "gurls++/optmatrix.h"

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

    T *diff = new T[pred.getSize()];
    copy(diff, pred.getData(), pred.getSize());


//     diff 	= opt.pred - y;
    axpy(rows*cols, (T)-1.0, y_true, 1, diff, 1);

//    p.rmse = sqrt(sum(diff.^2,1)/n);
    gMat2D<T> *rmse = new gMat2D<T>(1, cols);
    T* r_it = rmse->getData();
    const T* diff_it = diff;

    for(unsigned long i=0; i< cols; ++i, ++r_it, diff_it+=rows)
	{
		*r_it = nrm2(rows, diff_it, 1);
		*r_it /= sqrt((T)rows);
	}
	
    delete [] diff;

    perf->addOpt("rmse", new OptMatrix<gMat2D<T> >(*rmse));

//    p.forho 	= -p.rmse;
    gMat2D<T> *forho = new gMat2D<T>(1, cols);
    set(forho->getData(), (T)0.0, cols);
    axpy(cols, (T)-1.0, rmse->getData(), 1, forho->getData(), 1);
    perf->addOpt("forho", new OptMatrix<gMat2D<T> >(*forho));


    return perf;
}

}

#endif //_GURLS_RMSE_H_
