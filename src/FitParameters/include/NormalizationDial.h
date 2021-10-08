//
// Created by Nadrino on 26/05/2021.
//

#ifndef XSLLHFITTER_NORMALIZATIONDIAL_H
#define XSLLHFITTER_NORMALIZATIONDIAL_H

#include "Dial.h"

class NormalizationDial : public Dial {

public:
  NormalizationDial();

  void reset() override;
  void initialize() override;

  double evalResponse(const double& parameterValue_) override;

protected:
  void fillResponseCache() override;

};


#endif //XSLLHFITTER_NORMALIZATIONDIAL_H