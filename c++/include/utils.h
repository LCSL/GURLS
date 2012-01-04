/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011, IIT@MIT Lab
 * All rights reserved.
 *
 * author:  M. Santoro
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


#ifndef _GURLS_UTILS_H_
#define _GURLS_UTILS_H_

#include <map>
#include <algorithm>

#include "gmath.h"
#include "gmat2d.h"

using namespace std;

namespace gurls {

/*
 * compute the average precision through precision and recall
 */
template <typename T>
double precrec_driver(T* out, T* gt, unsigned long N){

	map<T, T> data;
	for (unsigned long i = 0; i < N; i++){
		data[-out[i]] = gt[i];
	}

	T tp[N], fp[N], prec[N], rec[N];
	unsigned long tpcumsum = 0;
	unsigned long fpcumsum = 0;

	typename map<T, T>::iterator it = data.begin();
	typename map<T, T>::iterator end = data.end();

	unsigned long idx = 0;
	while (it != end) {
		tp[idx] = it->second > 0 ? tpcumsum+=1 : tpcumsum;
		fp[idx] = it->second < 0 ? fpcumsum+=1 : fpcumsum;
		it++;
		idx++;
	}

	for (idx = 0; idx < N; idx++){
		rec[idx] = tp[idx]/tpcumsum;
		prec[idx] = tp[idx]/(tp[idx]+fp[idx]);
	}

	// compute average precision
	T ap = 0;
	T p = 0.0;
	T t = 0.f;
	int steps = 0;
	while ( t <= 1.0f ) {
		p = 0;
		for (idx = 0; idx < N; idx++){
			if (rec[idx] >= t){
				p = std::max(p, prec[idx]);
			}
		}
		ap+=p;
		steps++;
		t+=.1;
	}
	return ap/=steps;
}

}


#endif // _GURLS_UTILS_H_
