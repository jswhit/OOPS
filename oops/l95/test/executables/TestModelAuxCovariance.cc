/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/L95Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/ModelAuxCovariance.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::ModelAuxCovariance<lorenz95::L95Traits> tests;
  run.execute(tests);
  return 0;
};

