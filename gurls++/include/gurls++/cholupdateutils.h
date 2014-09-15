/*
 * The GURLS Package in C++
 *
 * Copyright (C) 2011-1013, IIT@MIT Lab
 * All rights reserved.
 *
 * author:  Raffaello Camoriano, Arjan Gijsberts
 * email:   raffaello.camoriano@iit.it
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

// Content: Utility functions used by the primalrecupdatecholesky class,
// performing the rank-k update of the Cholesky factor R

#include "gurls++/gmat2d.h"

namespace gurls {

void dchud(double* r, int ldr, int p, double* x, double* z, int ldz, int nz,
           double* y, double* rho, double* c, double* s,
           unsigned char rtrans, unsigned char ztrans) {
  
    // WARNING: This function only updates R. The update of z and rho is not implemented.
    //
    // Performs the rank-1 update of the Cholesky factor r with the input vector x by applying p Givens rotations.
    //
    // Arguments:
    //
    // r:           Upper Cholesky factor buffer
    // ldr:         Unknown, not used (probably it's the order of r, which must be == p)
    // p:           Dimension of vector x and order of matrix r
    // x:           Update sample buffer
    // z:           Known terms matrix
    // ldz:         Number of columns of z
    // nz:          Number of rows of z
    // y:           Output labels vector
    // rho:         On entry, the norms of the residual vectors that are to be updated. On exit, RHO has been updated. If RHO(i) is negative on entry then it is not changed.
    // c:           Cosines of the transforming rotations
    // s:           Sines of the transforming rotations
    // rtrans:      Bool to transpose r
    // ztrans:      Bool to transpose z    
    
    unsigned int i;
    unsigned int stp;   // Step to apply to r to get the next row element
    double* work = (double*) 0x0;   // Pointer to the elements of the update sample
    double* tbuff = (double*) 0x0;  // Pointer to the elements of the diagonal of R

    // create working copy of x
    work = (double*) malloc(p * sizeof(double));
    int one = 1;
    dcopy_(&p, x, &one, work, &one);

    stp = (rtrans == 1) ? p : 1;

    // Update r and fill c and s on the go
    for(i = 0, tbuff = r; (int)i < p; tbuff += (p+1), i++) {
        
        // Compute i-th Givens rotation
        drotg_(tbuff, work+i, c+i, s+i);   

        // force positive values on the diagonal (not strictly necessary)
        //if(*tbuff < 0) {
        //    *tbuff = -(*tbuff);
        //    c[i] = -c[i];
        //    s[i] = -s[i];
        //}

        // Apply Givens rotation
        if((int)i < p - 1) {
            
            int n = p - i - 1;
            int one = 1;
            int step = stp;
            double* xcoord = tbuff + stp;
            double* ycoord = work + i + 1;
            double* costh = c + i;
            double* sinth = s + i;
            
            drot_(&n, xcoord, &step, ycoord, &one, costh, sinth);
            
        }
    }
    
    free(work);
        
    // NOTE: The update of z and rho is not implemented
// update z and rho if applicable
//     if(nz > 0) {
//         work = (double*) malloc(nz * sizeof(double));
//         dcopy_(&nz, y, &one, work, &one); 
// 
//         stp = (ztrans == 1) ? 1 : nz;
//         stp2 = (ztrans == 1) ? ldz : 1;
// 
//         // update z
//         for(i = 0; (int)i < p; i++) {
//             //drot_(nz, z+i*stp, stp2, work, 1, c[i], s[i]);
//         }
// 
//         // update rho
//         for(i = 0; (int)i < nz; i++) {
//             if(work[i] < 0) {
//                 work[i] = -work[i];
//             }
//             if(rho[i] > 0) {
//                 scale = work[i] + rho[i];
//                 workscale = work[i] / scale;
//                 rhoscale = rho[i] / scale;
//                 rho[i] = scale * sqrt((workscale)*(workscale) + (rhoscale)*(rhoscale));
//                 // above is equal to line below, but less sensitive to numerical instability
//                 //rho[i] = sqrt(work[i]*work[i] + rho[i]*rho[i]);
//             }
//         }
//         free(work);
//     }
}

// WARNING: cholupdate only works with double parameters, not generic typename T
void cholupdate(gMat2D<double>& R, double& x, bool rtrans) {  
   
    // rPtr points to the vectorized R matrix
    double* rPtr = R.getData();
    
    int d = (int) R.cols();
    
    double* xPtr = &x;          // Initialize pointer to the x vector
    double* c = new double[d];  
    double* s = new double[d];  

    // dchud should be called in place of gsl_linalg_cholesky_update, since GSL is not available in GURLS
    dchud(rPtr, 0, d, xPtr, 0, 0, 0,
          0, 0, c, s,
          (unsigned char) rtrans, 0);
    
    // Store the updated matrix in R    
    R = gMat2D<double>(rPtr, (unsigned long) d, (unsigned long) d, 1);

    // Free memory
    delete[] c;
    delete[] s;
}
}