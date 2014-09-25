#ifndef _GURLS_ABSERR_H_
#define _GURLS_ABSERR_H_

#include "gurls++/perf.h"

#include "gurls++/utils.h"
#include "gurls++/gvec.h"
#include "gurls++/optmatrix.h"

#include "float.h"

namespace gurls {

/**
 * \ingroup Performance
 * \brief PerfAbsErr is the sub-class of Performance that evaluates prediction error
 */

template <typename T>
class PerfAbsErr: public Performance<T>{

public:
	///
	/// Default constructor
	///
	PerfAbsErr():Performance<T>("abserr"){}
	
	///
	/// Clone method
	///
	TaskBase *clone()
	{
		return new PerfAbsErr<T>();
	}

    /**
     * Evaluates the absolute mean error of the predicted labels stored in the field pred of opt with respect to the true input labels Y. 
     * \param X not used
     * \param Y labels matrix
     * \param opt options with the following:
     *  - pred (settable with the class Prediction and its subclasses)
     *
     * \return perf, a GurslOptionList equal to the field pred of opt, with the following fields added or substituted:
     *  - abserr = absolute mean error for each class/task
     *  - forho = -abserr
     */
    GurlsOptionsList* execute(const gMat2D<T>& X, const gMat2D<T>& Y, const GurlsOptionsList& opt) throw(gException);
};

template<typename T>
GurlsOptionsList* PerfAbsErr<T>::execute(const gMat2D<T>& /*X*/, const gMat2D<T>& Y, const GurlsOptionsList &opt) throw(gException)
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

        perf->removeOpt("abserr");
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

//    p.abserr = sum(abs(diff),1);
	gMat2D<T> *abserr = new gMat2D<T>(1, cols);
	T* r_it = abserr->getData();
	const T* diff_it = diff;

    for(unsigned long i=0; i< cols; ++i, ++r_it, diff_it+=rows){
		*r_it = static_cast<T>(0.0);
		for (unsigned long j=0; j < rows; ++j)
			*r_it += (*(diff_it+j)) * ((*(diff_it+j) >= 0.0) ? 1.0 : -1.0);
		*r_it /= rows;
	}
		
    delete [] diff;

    perf->addOpt("abserr", new OptMatrix<gMat2D<T> >(*abserr));

//    p.forho 	= -p.abserr;
    gMat2D<T> *forho = new gMat2D<T>(1, cols);
    set(forho->getData(), (T)0.0, cols);
    axpy(cols, (T)-1.0, abserr->getData(), 1, forho->getData(), 1);
    perf->addOpt("forho", new OptMatrix<gMat2D<T> >(*forho));


    return perf;
}

}

#endif //_GURLS_ABSERR_H_
