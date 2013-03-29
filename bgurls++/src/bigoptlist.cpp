#include "bigoptlist.h"

namespace gurls
{

BGurlsOptionsList::BGurlsOptionsList(const std::string &ExpName, const std::string &sharedDir, bool usedefopt): GurlsOptionsList(ExpName, usedefopt)
{
    init(sharedDir, usedefopt);
}

BGurlsOptionsList::BGurlsOptionsList(const std::string &ExpName, const std::wstring &sharedDir, bool usedefopt): GurlsOptionsList(ExpName, usedefopt)
{
    std::string dir = std::string(sharedDir.begin(), sharedDir.end());

    init(dir, usedefopt);
}

void BGurlsOptionsList::init(const std::string& sharedDir, bool usedefopt)
{
    using namespace boost::filesystem;

    if(usedefopt)
    {
        (*table)["nb_pred"] = new OptNumber(1);
        (*table)["memlimit"] = new OptNumber(std::pow(2.0, 30)); // default 1 GB

        (*table)["shared_dir"] = new OptString(sharedDir);

        path sharedDirPath(sharedDir);

        (*table)["tmpfile"] = new OptString(path(sharedDirPath / "tmp").native());

        GurlsOptionsList* files = new GurlsOptionsList("files");
        files->addOpt("Xva_filename", path(sharedDirPath / "Xva.h5").native());
        files->addOpt("Yva_filename", path(sharedDirPath / "Yva.h5").native());
        files->addOpt("XvatXva_filename", path(sharedDirPath / "XvatXva.h5").native());
        files->addOpt("XvatYva_filename", path(sharedDirPath / "XvatYva.h5").native());
        files->addOpt("XtX_filename", path(sharedDirPath / "XtX.h5").native());
        files->addOpt("Xty_filename", path(sharedDirPath / "Xty.h5").native());
        files->addOpt("pred_filename", path(sharedDirPath / "pred.h5").native());
        files->addOpt("optimizer_W_filename", path(sharedDirPath / "optimizer_W.h5").native());
        files->addOpt("optimizer_C_filename", path(sharedDirPath / "optimizer_C.h5").native());
        files->addOpt("optimizer_X_filename", path(sharedDirPath / "optimizer_X.h5").native());

        (*table)["files"] = files;

    }
}

}
