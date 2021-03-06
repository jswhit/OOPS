/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_QGIDTLM_H_
#define QG_MODEL_QGIDTLM_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "oops/interface/LinearModelBase.h"

#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "model/QgTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
// -----------------------------------------------------------------------------
/// QG linear identity model definition.
/*!
 *  QG linear identity model definition and configuration parameters.
 */

class QgIdTLM: public oops::LinearModelBase<QgTraits>,
              private util::ObjectCounter<QgIdTLM> {
 public:
  static const std::string classname() {return "qg::QgIdTLM";}

  QgIdTLM(const QgGeometry &, const eckit::Configuration &);
  ~QgIdTLM();

/// Model trajectory computation
  void setTrajectory(const QgState &, QgState &, const ModelBias &) override;

/// Run TLM and its adjoint
  void initializeTL(QgIncrement &) const override;
  void stepTL(QgIncrement &, const ModelBiasIncrement &) const override;
  void finalizeTL(QgIncrement &) const override;

  void initializeAD(QgIncrement &) const override;
  void stepAD(QgIncrement &, ModelBiasIncrement &) const override;
  void finalizeAD(QgIncrement &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const QgGeometry & resolution() const {return resol_;}

 private:
  void print(std::ostream &) const override;

// Data
  int keyConfig_;
  util::Duration tstep_;
  const QgGeometry resol_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_QGIDTLM_H_
