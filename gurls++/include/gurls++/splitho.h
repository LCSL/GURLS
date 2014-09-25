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


#ifndef _GURLS_SPLITHO_H_
#define _GURLS_SPLITHO_H_


#include "gurls++/split.h"
#include "gurls++/gmath.h"

namespace gurls {

/**
 * \ingroup Split
 * \brief SplitHoMulti is the sub-class of Split that splits data into one or more pairs of training and test samples.
 */

template <typename T>
class SplitHo: public Split<T>
{
public:
	///
	/// Default constructor
	///
	SplitHo():Split<T>("ho"){}
	
	///
	/// Clone method
	///
	TaskBase *clone()
	{
		return new SplitHo<T>();
	}

    /**
     * Splits data into one or more pairs of training and test samples, to be used for cross-validation. The fraction of samples for the validation set is specified in the field hoproportion of opt, and the number of pairs is specified in the field nholdouts of opt
     * \param X not used
     * \param Y labels matrix
     * \param opt options with the following field
     *   - hoproportion (default)
     *   - nholdouts (default)
     *
     * \return adds to opt the field split, which is a list containing the following fields:
     *  - indices = n x nholdouts matrix, each column contains the indices of training and validation samples
     *  - lasts = nholdout x 1 array, each row contains the number of elements of training set, which will be built taking the samples corresponding to the first lasts+1 elements of indices, the remainder indices will be used for validation.
     */

    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList *SplitHo<T>::execute(const gMat2D<T>& /*X*/, const gMat2D<T>& Y, const GurlsOptionsList &opt) throw(gException)
{
//    nSplits = opt.nholdouts;
    const int nSplits = static_cast<int>(opt.getOptAsNumber("nholdouts"));

//    fraction = opt.hoproportion;
    const double fraction = opt.getOptAsNumber("hoproportion");

//    [n,T] = size(y);
    const int n = Y.rows();
    const int t = Y.cols();
	
    gMat2D<unsigned long>* m_indices = new gMat2D<unsigned long>(n, nSplits);
    unsigned long* indices = m_indices->getData();
	
    gMat2D<unsigned long>* m_lasts = new gMat2D<unsigned long>(nSplits, 1);

	//check if Y is in gurls classification format
	bool islabel = isLabel(Y);

	if(islabel)
	{

	int nclass=t;
	if(t==1)
		nclass=2;

//    [dummy, y] = max(y,[],2);
    T* work = new T[Y.getSize()];
    unsigned long* y = new unsigned long[n];

//    maxPerRow(Y.getData(), n, t, y);
    indicesOfMax(Y.getData(), n, t, y, work, 2);

    delete[] work;

//    for t = 1:T,
//        classes{t}.idx = find(y == t);
//    end

    int* nSamples= new int[nclass];
    int nva = 0;
    unsigned long* idx = new unsigned long[n];
    unsigned long* it_idx = idx;

//    for t = 1:T,
    for(int i=0; i<t; ++i)
    {
//        nSamples(t) = numel(classes{t}.idx);
//        nva = nva + floor(fraction*nSamples(t));
		if(t>1){
			indicesOfEqualsTo<unsigned long>(y, n, i, it_idx, nSamples[i]);
			nva += static_cast<int>(std::floor(fraction*nSamples[i]));
			it_idx += nSamples[i];
		}
		else{
			indicesOfEqualsTo(Y.getData(), n, Y.min(), it_idx, nSamples[0]);
			indicesOfEqualsTo(Y.getData(), n, Y.max(), it_idx+nSamples[0], nSamples[1]);
			nva =  static_cast<int>(std::floor(fraction*(nSamples[0]))+std::floor(fraction*(nSamples[1])));
			}
    }

	set(m_lasts->getData(), (unsigned long) (n-nva), nSplits);
    delete[] y;


    int count_tr;
    int count_va;

//    for state = 1:nSplits,
    for(int state=0; state<nSplits; ++state)
    {
//    %% Shuffle each class
        count_tr = 0;
        count_va = n-nva;

        it_idx = idx;

//        for t = 1:T,
        for(int i=0; i<nclass; ++i)
        {
            const int nsamples = nSamples[i];

            if(nsamples == 0)
                continue;

            randperm(nsamples, it_idx, false);

            int last = nsamples - static_cast<int>(std::floor(nsamples*fraction));

            copy(indices+count_tr+state*n, it_idx, last, 1, 1);

            copy(indices+count_va+state*n, it_idx+last, nsamples-last, 1, 1);
			
            count_tr += last;
            count_va += nsamples-last;
            it_idx += nsamples;

        }

    }
	delete[] idx;
	delete[] nSamples;
	}
	else
	{
		int nva =  static_cast<int>(std::floor(fraction*(n)));
		set(m_lasts->getData(), (unsigned long) (n-nva), nSplits);
		
		unsigned long* idx = new unsigned long[n];
		unsigned long* it_idx = idx;

		for(int state=0; state<nSplits; ++state)
			{
            randperm(n, it_idx, true, 0);
            copy(indices+state*n, it_idx, n, 1, 1);
			}
		delete[] idx;
	}
    GurlsOptionsList* split = new GurlsOptionsList("split");
    split->addOpt("indices", new OptMatrix<gMat2D<unsigned long> >(*m_indices));
    split->addOpt("lasts", new OptMatrix<gMat2D<unsigned long> >(*m_lasts));

    return split;
}
}

#endif //_GURLS_SPLITHO_H_