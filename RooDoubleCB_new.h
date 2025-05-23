/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef RooDoubleCB_new_h
#define RooDoubleCB_new_h

#include <RooAbsCategory.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooCategoryProxy.h>
#include <RooRealProxy.h>
#include <TMath.h>

#include <complex>

class RooDoubleCB_new : public RooAbsPdf {
public:
  RooDoubleCB_new() {}
  RooDoubleCB_new(const char *name, const char *title, RooAbsReal &_m,
                  RooAbsReal &_m0, RooAbsReal &_sigma, RooAbsReal &_alpha1,
                  RooAbsReal &_n1, RooAbsReal &_alpha2, RooAbsReal &_n2);
  RooDoubleCB_new(RooDoubleCB_new const &other, const char *name = nullptr);
  TObject *clone(const char *newname) const override {
    return new RooDoubleCB_new(*this, newname);
  }

protected:
  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;

  double evaluate() const override;
  void doEval(RooFit::EvalContext &) const override;
  void translate(RooFit::Detail::CodeSquashContext &ctx) const override;

private:
  ClassDefOverride(RooDoubleCB_new, 1) // Your description goes here...
};
inline double RooDoubleCB_new_evaluate(double m, double m0, double sigma,
                                       double alpha1, double n1, double alpha2,
                                       double n2) {
  Double_t t = (m - m0) / sigma;
  if (alpha1 < 0)
    t = -t;

  Double_t absAlpha = fabs((Double_t)alpha1); //fabs ok
  Double_t absAlpha2 = fabs((Double_t)alpha2);

  if (t >= -absAlpha && t < absAlpha2) {
    return exp(-0.5 * t * t); //exp ok
  }
  if (t < -absAlpha) {
    Double_t a =
        TMath::Power(n1 / absAlpha, n1) * exp(-0.5 * absAlpha * absAlpha);
    Double_t b = n1 / absAlpha - absAlpha;
    return a / TMath::Power(b - t, n1);
  }
  if (t >= absAlpha2) {
    Double_t c =
        TMath::Power(n2 / absAlpha2, n2) * exp(-0.5 * absAlpha2 * absAlpha2);
    Double_t d = n2 / absAlpha2 - absAlpha2;
    return c / TMath::Power(d + t, n2);
  }
  return 0.0;
}

#endif // RooDoubleCB_new_h