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


#ifndef _GURLS_PRECISIONRECALL_H_
#define _GURLS_PRECISIONRECALL_H_

#include "gurls++/perf.h"

#include "gurls++/utils.h"
#include "gurls++/gvec.h"
#include "gurls++/optmatrix.h"

namespace gurls {

/**
 * \ingroup Performance
 * \brief PerfPrecRec is the sub-class of Performance that evaluates prediction precision
 */

template <typename T>
class PerfPrecRec: public Performance<T>{

public:
	///
	/// Default constructor
	///
	PerfPrecRec():Performance<T>("precrec"){}
	
	///
	/// Clone method
	///
	TaskBase *clone()
	{
		return new PerfPrecRec<T>();
	}

    /**
     * Evaluates the average precision per class through precision and recall.
     *
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt options with the following:
     *  - pred (settable with the class Prediction and its subclasses)
     *
     * \return perf, a GurslOptionList equal to the field pred of opt, with the following fields added or substituted:
     *  - ap = array of average precision for each class
     *  - forho = ap
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList* PerfPrecRec<T>::execute(const gMat2D<T>& /*X*/, const gMat2D<T>& Y, const GurlsOptionsList& opt) throw(gException)
{
    const int rows = Y.rows();
    const int cols = Y.cols();

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

        perf->removeOpt("ap");
        perf->removeOpt("forho");
//        perf->removeOpt("forplot");
    }
    else
        perf = new GurlsOptionsList("perf");

//    y_true = y;
    const T* y_true = Y.getData();

//    y_pred = opt.pred;
    const gMat2D<T> &y_pred = opt.getOptValue<OptMatrix<gMat2D<T> > >("pred");

    gMat2D<T>* ap_mat = new gMat2D<T>(1, cols);
    T* ap = ap_mat->getData();

    T* work = new T[4*rows];

//    T = size(y,2);
//    for t = 1:T,
    for(int i=0; i<cols; ++i)
    {
//        p.ap(t) = precrec_driver(y_pred(:,t), y_true(:,t),0);
//        p.forho(t) = p.ap(t);

        ap[i] = precrec_driver(y_pred.getData()+(i*rows), y_true+(i*rows), rows, work);

//        p.forplot(t) = p.ap(t);
    }

    delete [] work;

    OptMatrix<gMat2D<T> >* ap_opt = new OptMatrix<gMat2D<T> >(*ap_mat);
    perf->addOpt("ap", ap_opt);

    OptMatrix<gMat2D<T> >* forho_opt = new OptMatrix<gMat2D<T> >(*(new gMat2D<T>(*ap_mat)));
    perf->addOpt("forho", forho_opt);

//    OptMatrix<gMat2D<T> >* forplot_opt = new OptMatrix<gMat2D<T> >(*(new gMat2D<T>(*ap_mat)));
//    perf->addOpt("forplot", forplot_opt);

    return perf;
}

}

#endif //_GURLS_PRECISIONRECALL_H_
