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

#ifndef _GURLS_BIGOPTLIST_H_
#define _GURLS_BIGOPTLIST_H_


#include "optlist.h"

#include <boost/filesystem/path.hpp>

namespace gurls
{


#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable : 4251)
#endif

/**
  * \ingroup Settings
  * \brief BGurlsOptionsList is an option containing a list of options
  * mapped by name.
  */
class GURLS_EXPORT BGurlsOptionsList: public GurlsOptionsList
{
public:
    /**
      * Constructor. Builds an optionlist with a name and optionally a set of default options
      *
      * \param ExpName name of the options list
      * \param usedefopt if \a true the list is filled with a set of default options, if \a false the list is left empty
      */
    BGurlsOptionsList(const std::string& ExpName, const std::string& sharedDir, bool usedefopt = false) : GurlsOptionsList(ExpName, usedefopt)
    {
        using namespace boost::filesystem3;

        if(usedefopt)
        {
            (*table)["nb_pred"] = new OptNumber(1);
            (*table)["memlimit"] = new OptNumber(std::pow(2, 30)); // default 1 GB

            (*table)["shared_dir"] = new OptString(sharedDir);

            path sharedDirPath(sharedDir);

            (*table)["tmpfile"] = new OptString(path(sharedDirPath / "tmp").native());

            GurlsOptionsList* files = new GurlsOptionsList("files");
            files->addOpt("Xva_filename", path(sharedDirPath / "Xva.nc").native());
            files->addOpt("Yva_filename", path(sharedDirPath / "Yva.nc").native());
            files->addOpt("XvatXva_filename", path(sharedDirPath / "XvatXva.nc").native());
            files->addOpt("XvatYva_filename", path(sharedDirPath / "XvatYva.nc").native());
            files->addOpt("XtX_filename", path(sharedDirPath / "XtX.nc").native());
            files->addOpt("Xty_filename", path(sharedDirPath / "Xty.nc").native());
            files->addOpt("pred_filename", path(sharedDirPath / "pred.nc").native());
            files->addOpt("optimizer_W_filename", path(sharedDirPath / "optimizzer_W.nc").native());
            files->addOpt("optimizer_C_filename", path(sharedDirPath / "optimizzer_C.nc").native());
            files->addOpt("optimizer_X_filename", path(sharedDirPath / "optimizzer_X.nc").native());

            (*table)["files"] = files;

        }
    }

};

#ifdef _WIN32
#pragma warning(pop)
#endif

}

#endif // _GURLS_BIGOPTLIST_H_
