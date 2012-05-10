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


#ifndef _GURLS_PRECISIONRECALL_H_
#define _GURLS_PRECISIONRECALL_H_

#include "perf.h"

#include "utils.h"
#include "gvec.h"

namespace gurls {

template <typename T>
class PrecisionRecall: public Performance<T>{

public:
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt)  throw(gException);
};

template<typename T>
void PrecisionRecall<T>::execute(const gMat2D<T>& /*X_OMR*/, const gMat2D<T>& Y_OMR, GurlsOptionsList& opt) throw(gException)
{
    const int rows = Y_OMR.rows();
    const int cols = Y_OMR.cols();

    gMat2D<T> Y(Y_OMR.cols(), Y_OMR.rows());
    Y_OMR.transpose(Y);

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

    if(perf->hasOpt("ap"))
        perf->removeOpt("ap");

    if(perf->hasOpt("forho"))
        perf->removeOpt("forho");

    T* ap = new T[cols];
    T* forho = new T[cols];

//    y_true = y;
    T* y_true = Y.getData();

//    y_pred = opt.pred;
    GurlsOption *pred_opt = opt.getOpt("pred");

//    if (pred_opt->getDataID() != typeid(T))
//            throw gException("Different types");

    gMat2D<T> *pred_mat = &(OptMatrix<gMat2D<T> >::dynacast(pred_opt))->getValue();
    gMat2D<float> y_pred(pred_mat->cols(), pred_mat->rows());
    pred_mat->transpose(y_pred);


//    T = size(y,2);
//    for t = 1:T,
    for(int i=0; i<cols; ++i)
    {
//        p.ap(t) = precrec_driver(y_pred(:,t), y_true(:,t),0);
//        p.forho(t) = p.ap(t);
        ap[i] = precrec_driver(y_pred.getData()+(i*rows), y_true+(i*rows), rows);

//        p.forplot(t) = p.ap(t);
        // TODO?
    }

    copy(forho, ap, cols);

    gMat2D<T>* ap_mat = new gMat2D<T>(ap, 1, cols, true);
    OptMatrix<gMat2D<T> >* ap_opt = new OptMatrix<gMat2D<T> >(*ap_mat);
    perf->addOpt("ap", ap_opt);

    gMat2D<T>* forho_mat = new gMat2D<T>(forho, 1, cols, true);
    OptMatrix<gMat2D<T> >* forho_opt = new OptMatrix<gMat2D<T> >(*forho_mat);
    perf->addOpt("forho", forho_opt);
}

}

#endif //_GURLS_PRECISIONRECALL_H_
