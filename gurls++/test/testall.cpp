#include "gmat2d.h"
#include "options.h"
#include "optlist.h"

#include "linearkernel.h"
#include "rbfkernel.h"
#include "chisquaredkernel.h"

#include "predkerneltraintest.h"

#include "precisionrecall.h"
#include "macroavg.h"
#include "rmse.h"

#include "rlsauto.h"
#include "rlsprimal.h"
#include "rlsprimalr.h"
#include "rlsdual.h"
#include "rlsdualr.h"
#include "rlspegasos.h"

#include "loocvprimal.h"
#include "loocvdual.h"
#include "fixlambda.h"
#include "fixsiglam.h"
#include "siglam.h"
#include "siglamho.h"
#include "hoprimal.h"
#include "hodual.h"

#include "primal.h"
#include "dual.h"

#include "splitho.h"

#include "boltzman.h"
#include "boltzmangap.h"
#include "gap.h"
#include "maxscore.h"

#include <cstdlib>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE yeast

#include <boost/test/unit_test.hpp>

#include "test.h"

using namespace gurls::test;

typedef double T;


//BOOST_AUTO_TEST_SUITE(yeast)

const path yeastDataPath( path(std::string(GURLS_DATA_DIR))/path("yeast"));


//BOOST_AUTO_TEST_CASE(TestSplitHo)
//{
//    Data data(yeastDataPath, "splitho", true);

//    data.loadDefaults();

//    Fixture<T, gurls::SplitHo<T> >fixture(yeastDataPath, "splitho", data);

//    fixture.runTask();

//    fixture.checkResults("split");
//}

BOOST_AUTO_TEST_CASE(TestParamSelSiglam)
{
    Data data(yeastDataPath, "paramselsiglam", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelSiglam<T> > fixture(yeastDataPath, "paramselsiglam", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

BOOST_AUTO_TEST_CASE(TestParamSelSiglamHo)
{
    Data data(yeastDataPath, "paramselsiglamho", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelSiglamHo<T> > fixture(yeastDataPath, "paramselsiglamho", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

BOOST_AUTO_TEST_CASE(TestKernelLinear)
{
    Data data(yeastDataPath, "kernellinear", true);

    data.loadDefaults();

    Fixture<T, gurls::KernelLinear<T> > fixture(yeastDataPath, "kernellinear", data);

    fixture.runTask();

    fixture.checkResults("kernel");
}

BOOST_AUTO_TEST_CASE(TestKernelRBF)
{
    Data data(yeastDataPath, "kernelgauss", true);

    data.loadDefaults();

    Fixture<T, gurls::KernelRBF<T> > fixture(yeastDataPath, "kernelgauss", data);

    fixture.runTask();

    fixture.checkResults("kernel");
}

BOOST_AUTO_TEST_CASE(TestKernelChisquared)
{
    Data data(yeastDataPath, "kernelchisquared", true);

    data.loadDefaults();

    Fixture<T, gurls::KernelChisquared<T> > fixture(yeastDataPath, "kernelchisquared", data);

    fixture.runTask();

    fixture.checkResults("kernel");
}

BOOST_AUTO_TEST_CASE(TestParamSelLoocvPrimal)
{
    Data data(yeastDataPath, "paramselloocvprimal", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelLoocvPrimal<T> > fixture(yeastDataPath, "paramselloocvprimal", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

BOOST_AUTO_TEST_CASE(TestParamSelHoPrimal)
{
    Data data(yeastDataPath, "paramselhoprimal", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelHoPrimal<T> > fixture(yeastDataPath, "paramselhoprimal", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

BOOST_AUTO_TEST_CASE(TestParamSelLoocvDual_linearkernel)
{
    Data data(yeastDataPath, "paramselloocvdual_linear", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelLoocvDual<T> > fixture(yeastDataPath, "paramselloocvdual_linear", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

BOOST_AUTO_TEST_CASE(TestParamSelHoDual_linearkernel)
{
    Data data(yeastDataPath, "paramselhodual_linear", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelHoDual<T> > fixture(yeastDataPath, "paramselhodual_linear", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

BOOST_AUTO_TEST_CASE(TestParamSelLoocvDual_gausskernel)
{
    Data data(yeastDataPath, "paramselloocvdual_gauss", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelLoocvDual<T> > fixture(yeastDataPath, "paramselloocvdual_gauss", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

BOOST_AUTO_TEST_CASE(TestParamSelHoDual_gausskernel)
{
    Data data(yeastDataPath, "paramselhodual_gauss", true);

    data.loadDefaults();

    Fixture<T, gurls::ParamSelHoDual<T> > fixture(yeastDataPath, "paramselhodual_gauss", data);

    fixture.runTask();

    fixture.checkResults("paramsel");
}

//BOOST_AUTO_TEST_CASE(TestParamSelHoPrimalr)
//{
//    Data data(yeastDataPath, "paramselhoprimalr", true);

//    data.loadDefaults();

//    Fixture<T, gurls::ParamSelHoPrimalr<T> > fixture(yeastDataPath, "paramselhoprimalr", data);

//    fixture.runTask();

//    fixture.checkResults("paramsel");
//}

//BOOST_AUTO_TEST_CASE(TestParamSelHoDualr_linearkernel)
//{
//    Data data(yeastDataPath, "paramselhodualr_linear", true);

//    data.loadDefaults();

//    Fixture<T, gurls::ParamSelHoDualr<T> > fixture(yeastDataPath, "paramselhodualr_linear", data);

//    fixture.runTask();

//    fixture.checkResults("paramsel");
//}

//BOOST_AUTO_TEST_CASE(TestParamSelHoDualr_gausskernel)
//{
//    Data data(yeastDataPath, "paramselhodualr_gauss", true);

//    data.loadDefaults();

//    Fixture<T, gurls::ParamSelHoDualr<T> > fixture(yeastDataPath, "paramselhodualr_gauss", data);

//    fixture.runTask();

//    fixture.checkResults("paramsel");
//}


BOOST_AUTO_TEST_CASE(TestRLSPrimal)
{
    Data data(yeastDataPath, "rlsprimal", true);

    data.loadDefaults();

    Fixture<T, gurls::RLSPrimal<T> >fixture(yeastDataPath, "rlsprimal", data);

    fixture.runTask();

    fixture.checkResults("optimizer");
}

BOOST_AUTO_TEST_CASE(TestRLSDual_linearkernel)
{
    Data data(yeastDataPath, "rlsdual_linear", true);

    data.loadDefaults();

    Fixture<T, gurls::RLSDual<T> >fixture(yeastDataPath, "rlsdual_linear", data);

    fixture.runTask();

    fixture.checkResults("optimizer");
}

BOOST_AUTO_TEST_CASE(TestRLSDual_gausskernel)
{
    Data data(yeastDataPath, "rlsdual_gauss", true);

    data.loadDefaults();

    Fixture<T, gurls::RLSDual<T> >fixture(yeastDataPath, "rlsdual_gauss", data);

    fixture.runTask();

    fixture.checkResults("optimizer");
}


//BOOST_AUTO_TEST_CASE(TestRLSPrimalr)
//{
//    Data data(yeastDataPath, "rlsprimalr", true);

//    data.loadDefaults();

//    Fixture<T, gurls::RLSPrimalr<T> >fixture(yeastDataPath, "rlsprimalr", data);

//    fixture.runTask();

//    fixture.checkResults("optimizer");
//}

//BOOST_AUTO_TEST_CASE(TestRLSDualr_linearkernel)
//{
//    Data data(yeastDataPath, "rlsdualr_linear", true);

//    data.loadDefaults();

//    Fixture<T, gurls::RLSDualr<T> >fixture(yeastDataPath, "rlsdualr_linear", data);

//    fixture.runTask();

//    fixture.checkResults("optimizer");
//}

//BOOST_AUTO_TEST_CASE(TestRLSDualr_gausskernel)
//{
//    Data data(yeastDataPath, "rlsdualr_gauss", true);

//    data.loadDefaults();

//    Fixture<T, gurls::RLSDualr<T> > fixture(yeastDataPath, "rlsdualr_gauss", data);

//    fixture.runTask();

//    fixture.checkResults("optimizer");
//}


BOOST_AUTO_TEST_CASE(TestPredKernelTrainTest_chisquaredkernel)
{
    Data data(yeastDataPath, "predkernelchisquared", true);

    data.loadDefaults();

    Fixture<T, gurls::PredKernelTrainTest<T> > fixture(yeastDataPath, "predkernelchisquared", data);

    fixture.runTask();

    fixture.checkResults("predkernel");
}

BOOST_AUTO_TEST_CASE(TestPredKernelTrainTest_gausskernel)
{
    Data data(yeastDataPath, "predkernel_gauss", true);

    data.loadDefaults();

    Fixture<T, gurls::PredKernelTrainTest<T> > fixture(yeastDataPath, "predkernel_gauss", data);

    fixture.runTask();

    fixture.checkResults("predkernel");
}

BOOST_AUTO_TEST_CASE(TestPredPrimal)
{
    Data data(yeastDataPath, "predprimal", true);

    data.loadDefaults();

    Fixture<T, gurls::PredPrimal<T> > fixture(yeastDataPath, "predprimal", data);

    fixture.runTask();

    fixture.checkResults("pred");
}

BOOST_AUTO_TEST_CASE(TestPredDual_linearkernel)
{
    Data data(yeastDataPath, "preddual_linear", true);

    data.loadDefaults();

    Fixture<T, gurls::PredDual<T> > fixture(yeastDataPath, "preddual_linear", data);

    fixture.runTask();

    fixture.checkResults("pred");
}

BOOST_AUTO_TEST_CASE(TestPredDual_gausskernel)
{
    Data data(yeastDataPath, "preddual_gauss", true);

    data.loadDefaults();

    Fixture<T, gurls::PredDual<T> > fixture(yeastDataPath, "preddual_gauss", data);

    fixture.runTask();

    fixture.checkResults("pred");
}

BOOST_AUTO_TEST_CASE(TestPerfRmse)
{
    Data data(yeastDataPath, "perfrmse", true);

    data.loadDefaults();

    Fixture<T, gurls::PerfRmse<T> >fixture(yeastDataPath, "perfrmse", data);

    fixture.runTask();

    fixture.checkResults("perf");
}

BOOST_AUTO_TEST_CASE(TestPerfMacroAvg)
{
    Data data(yeastDataPath, "perfmacroavg", true);

    data.loadDefaults();

    Fixture<T, gurls::PerfMacroAvg<T> >fixture(yeastDataPath, "perfmacroavg", data);

    fixture.runTask();

    fixture.checkResults("perf");
}

BOOST_AUTO_TEST_CASE(TestPerfPrecRec)
{
    Data data(yeastDataPath, "perfprecrec", true);

    data.loadDefaults();

    Fixture<T, gurls::PerfPrecRec<T> >fixture(yeastDataPath, "perfprecrec", data);

    fixture.runTask();

    fixture.checkResults("perf");
}

BOOST_AUTO_TEST_CASE(TestConfGap)
{
    Data data(yeastDataPath, "confgap", true);

    data.loadDefaults();

    Fixture<T, gurls::ConfGap<T> >fixture(yeastDataPath, "confgap", data);

    fixture.runTask();

    fixture.checkResults("conf");
}

BOOST_AUTO_TEST_CASE(TestConfMaxScore)
{
    Data data(yeastDataPath, "confmaxscore", true);

    data.loadDefaults();

    Fixture<T, gurls::ConfMaxScore<T> >fixture(yeastDataPath, "confmaxscore", data);

    fixture.runTask();

    fixture.checkResults("conf");
}

BOOST_AUTO_TEST_CASE(TestConfBoltzman)
{
    Data data(yeastDataPath, "confboltzman", true);

    data.loadDefaults();

    Fixture<T, gurls::ConfBoltzman<T> >fixture(yeastDataPath, "confboltzman", data);

    fixture.runTask();

    fixture.checkResults("conf");
}

BOOST_AUTO_TEST_CASE(TestConfBoltzmanGap)
{
    Data data(yeastDataPath, "confboltzmangap", true);

    data.loadDefaults();

    Fixture<T, gurls::ConfBoltzmanGap<T> >fixture(yeastDataPath, "confboltzmangap", data);

    fixture.runTask();

    fixture.checkResults("conf");
}

//BOOST_AUTO_TEST_SUITE_END()
