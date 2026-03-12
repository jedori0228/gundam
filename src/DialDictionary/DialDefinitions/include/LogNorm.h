//
// Created by Adrien Blanchet on 29/11/2022.
//

#ifndef GUNDAM_LOGNORM_H
#define GUNDAM_LOGNORM_H

#include "DialBase.h"

// for log-norm dial, no cache is needed


class LogNorm : public DialBase {

public:
  LogNorm() = default;

  [[nodiscard]] std::unique_ptr<DialBase> clone() const override { return std::make_unique<LogNorm>(*this); }
  [[nodiscard]] std::string getDialTypeName() const override { return {"lnN"}; }
  [[nodiscard]] double evalResponse(const DialInputBuffer& input_) const override { return exp(input_.getInputBuffer()[0]); }

  /// Build the dial with no input arguments.  This is here for completeness,
  /// but could eventually do... something.
  virtual void buildDial(const std::string& option_="") override {}

};


#endif //GUNDAM_LOGNORM_H
