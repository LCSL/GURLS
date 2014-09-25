/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-1013, Matteo Santoro
 * All rights reserved.
 *
 * author:  M. Santoro
 * email:   matteo.santoro@gmail.com
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

#ifndef _GURLS_PREDRANDFEATS_H
#define _GURLS_PREDRANDFEATS_H

#include <cmath>

#include "gurls++/pred.h"

#include "gurls++/gmath.h"
#include "gurls++/gmat2d.h"
#include "gurls++/options.h"
#include "gurls++/optlist.h"
#include "gurls++/utils.h"


namespace gurls {

/**
 * \ingroup Prediction
 * \brief PredRandFeats is the sub-class of Prediction that computes the predictions
 *  of the linear classifier stored in opt.rls.W, and obtained the Random Features approach proposed in:
 *  Ali Rahimi, Ben Recht;
 *  Random Features for Large-Scale Kernel Machines;
 *  in Neural Information Processing Systems (NIPS) 2007.
 *  on the samples passed in the X matrix.
 */
template <typename T>
class PredRandFeats: public Prediction<T> {

public:
	///
	/// Default constructor
	///
	PredRandFeats():Prediction<T>("randfeats"){}
	
	///
	/// Clone method
	///
	TaskBase *clone()
	{
		return new PredRandFeats<T>();
	}

    /**
     * computes the predictions of the linear classifier stored in opt.rls.W
     * \param X input data matrix
     * \param Y labels matrix
     * \param opt structure of options with the following fields (and subfields):
     *  - optimizer.W (set by the optimizer tasks)
     *  - optimizer.proj (set by the optimizer tasks)
     *
     * \return matrix of predicted labels
     */
    GurlsOptionsList *execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt);
};

template <typename T>
GurlsOptionsList *PredRandFeats<T>::execute(const gMat2D<T>& X, const gMat2D<T>& /*Y*/, const GurlsOptionsList &opt)
{

//    G = rp_apply_real(X, opt.rls.proj);
    const gMat2D<T>& proj = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.proj");
    gMat2D<T> *G = rp_apply_real(X, proj);

//    scores = G*opt.rls.W;
    const gMat2D<T>& W = opt.getOptValue<OptMatrix<gMat2D<T> > >("optimizer.W");

    gMat2D<T> *scores_mat = new gMat2D<T>(G->rows(), W.cols());
    dot(G->getData(), W.getData(), scores_mat->getData(), G->rows(), G->cols(), W.rows(), W.cols(), G->rows(), W.cols(), CblasNoTrans, CblasNoTrans, CblasColMajor);


    GurlsOptionsList* pred = new GurlsOptionsList("pred");
    pred->addOpt("scores", new OptMatrix<gMat2D<T> >(*scores_mat));

    return pred;
}

}

#endif // _GURLS_PREDRANDFEATS_H
