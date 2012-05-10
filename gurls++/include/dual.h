/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, Matteo Santoro
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


#ifndef _GURLS_DUAL_H
#define _GURLS_DUAL_H

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "gmath.h"
#include "options.h"
#include "optlist.h"

#include "pred.h"
#include "primal.h"

using namespace std;



namespace gurls {

template <typename T>
class PredDual: public Prediction<T> {

public:
    void execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt);
};

template <typename T>
void PredDual<T>::execute(const gMat2D<T>& X, const gMat2D<T>& Y, GurlsOptionsList& opt)
{
    try
    {

        GurlsOptionsList* kernel = static_cast<GurlsOptionsList*>(opt.getOpt("kernel"));

        if(kernel->getOptAsString("type") == "linear")
        {
            PredPrimal<T> pred;
            return pred.execute(X, Y, opt);
        }
        else
        {
            GurlsOptionsList* predkernel = static_cast<GurlsOptionsList*>(opt.getOpt("predkernel"));
            GurlsOption *K_opt = predkernel->getOpt("K");
            gMat2D<T> *K = &(OptMatrix<gMat2D<T> >::dynacast(K_opt))->getValue();


            GurlsOptionsList* rls = static_cast<GurlsOptionsList*>(opt.getOpt("rls"));
            GurlsOption *C_opt = rls->getOpt("C");
            gMat2D<T> *C = &(OptMatrix<gMat2D<T> >::dynacast(C_opt))->getValue();


            gMat2D<T>* Z = new gMat2D<T>(K->rows(), C->cols());
            set(Z->getData(), (T)0.0, Z->getSize());


            dot(K->getData(), C->getData(), Z->getData(), K->rows(), K->cols(), C->rows(), C->cols(), Z->rows(), Z->cols(), CblasNoTrans, CblasNoTrans, CblasRowMajor);


            if(opt.hasOpt("pred"))
                opt.removeOpt("pred");

             opt.addOpt("pred", new OptMatrix<gMat2D<T> >(*Z));
        }
    }
    catch (gException& ex)
    {
        throw ex;
    }

}

}



#endif // _GURLS_DUAL_H
