/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "RooDoubleCB.h"

#include <RooAbsReal.h>
#include <RooAbsCategory.h>

#include <Riostream.h>
#include <TMath.h>

#include <cmath>

ClassImp(RooDoubleCB);

RooDoubleCB::RooDoubleCB(const char *name, const char *title,
                        RooAbsReal& _m,
                        RooAbsReal& _m0,
                        RooAbsReal& _sigma,
                        RooAbsReal& _alpha1,
                        RooAbsReal& _n1,
                        RooAbsReal& _alpha2,
                        RooAbsReal& _n2)
   : RooAbsPdf(name,title),
   m("m","m",this,_m),
   m0("m0","m0",this,_m0),
   sigma("sigma","sigma",this,_sigma),
   alpha1("alpha1","alpha1",this,_alpha1),
   n1("n1","n1",this,_n1),
   alpha2("alpha2","alpha2",this,_alpha2),
   n2("n2","n2",this,_n2)
{
}

RooDoubleCB::RooDoubleCB(RooDoubleCB const &other, const char *name)
   : RooAbsPdf(other,name),
   m("m",this,other.m),
   m0("m0",this,other.m0),
   sigma("sigma",this,other.sigma),
   alpha1("alpha1",this,other.alpha1),
   n1("n1",this,other.n1),
   alpha2("alpha2",this,other.alpha2),
   n2("n2",this,other.n2)
{
}


double RooDoubleCB::evaluate() const 
{
   return RooDoubleCB_evaluate(m, m0, sigma, alpha1, n1, alpha2, n2); 
}

//void RooDoubleCB::computeBatch(double *output, std::size_t size, RooFit::Detail::DataMap const &dataMap) const
void RooDoubleCB::doEval(RooFit::EvalContext &ctx) const 
{ 
   /*std::span<const double> mSpan = dataMap.at(m);
   std::span<const double> m0Span = dataMap.at(m0);
   std::span<const double> sigmaSpan = dataMap.at(sigma);
   std::span<const double> alpha1Span = dataMap.at(alpha1);
   std::span<const double> n1Span = dataMap.at(n1);
   std::span<const double> alpha2Span = dataMap.at(alpha2);
   std::span<const double> n2Span = dataMap.at(n2);*/

   std::span<const double> mSpan = ctx.at(m);
   std::span<const double> m0Span = ctx.at(m0);
   std::span<const double> sigmaSpan = ctx.at(sigma);
   std::span<const double> alpha1Span = ctx.at(alpha1);
   std::span<const double> n1Span = ctx.at(n1);
   std::span<const double> alpha2Span = ctx.at(alpha2);
   std::span<const double> n2Span = ctx.at(n2);

   std::size_t n = ctx.output().size();

   //for (std::size_t i = 0; i < size; ++i) {
   for (std::size_t i = 0; i < n; ++i) {
      //output[i] = RooDoubleCB_evaluate(mSpan.size() > 1 ? mSpan[i] : mSpan[0],
      ctx.output()[i] = RooDoubleCB_evaluate(mSpan.size() > 1 ? mSpan[i] : mSpan[0],
                               m0Span.size() > 1 ? m0Span[i] : m0Span[0],
                               sigmaSpan.size() > 1 ? sigmaSpan[i] : sigmaSpan[0],
                               alpha1Span.size() > 1 ? alpha1Span[i] : alpha1Span[0],
                               n1Span.size() > 1 ? n1Span[i] : n1Span[0],
                               alpha2Span.size() > 1 ? alpha2Span[i] : alpha2Span[0],
                               n2Span.size() > 1 ? n2Span[i] : n2Span[0]);
   }
} 
void RooDoubleCB::translate(RooFit::Detail::CodeSquashContext &ctx) const
{
   ctx.addResult(this, ctx.buildCall("RooDoubleCB_evaluate", m, m0, sigma, alpha1, n1, alpha2, n2));
}
