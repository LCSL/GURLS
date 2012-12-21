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


#ifndef _GURLS_BIGMACROAVG_H_
#define _GURLS_BIGMACROAVG_H_

#include "optlist.h"
#include "bigperf.h"
#include "utils.h"
#include "optmatrix.h"
#include "bigarray.h"

namespace gurls {

/**
 * \ingroup Performance
 * \brief BigPerfMacroAvg is the sub-class of BigPerformance that evaluates prediction accuracy
 */

template <typename T>
class BigPerfMacroAvg: public BigPerformance<T>{

public:
    /**
     * Evaluates the average accuracy per class.
     *
     * \param X input data bigarray
     * \param Y labels bigarray
     * \param opt options with the following:
     *  - pred (settable with the class Prediction and its subclasses)
     *
     * \return perf, a GurslOptionList equal to the field pred of opt, with the following fields added or substituted:
     *  - acc = array of prediction accuracy for each class
     *  - forho = acc
     */
    GurlsOptionsList* execute(const BigArray<T>& X, const BigArray<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList* BigPerfMacroAvg<T>::execute(const BigArray<T>& /*X*/, const BigArray<T>& /*Y*/, const GurlsOptionsList& /*opt*/) throw(gException)
{
    // TODO

    return new GurlsOptionsList("perf");
}

}

#endif //_GURLS_BIGMACROAVG_H_
