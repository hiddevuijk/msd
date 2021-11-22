#ifndef GUARD_POTENTIAL_H
#define GUARD_POTENTIAL_H

#include "vec2.h"

#include <math.h>

class WCApotential {
 public:
  WCApotential(double sigma, double epsilon, double alpha)
    : sigma_(sigma), epsilon_(epsilon), alpha_(alpha) {}

  Vec2 force(const Vec2& r1, const Vec2& r2) const;
  Vec2 force(const Vec2& r1, const Vec2& r2, double L) const;
  double getSigmaCutOff() const { return sigmaCutOff_; }
 private:
  double sigma_;
  double epsilon_;
  double alpha_;

  double sigmaCutOff_;

};

Vec2 WCApotential::force(const Vec2& r1, const Vec2& r2) const
{
  Vec2 f(0,0);
  Vec2 d = r1 - r1;

  double l = sqrt(d.x * d.x + d.y * d.y);

  if ( l < sigmaCutOff_ ) {
    f = -epsilon_*pow(sigma_/l, alpha_) * d;
  }
  
  return f;
}


Vec2 WCApotential::force(const Vec2& r1, const Vec2& r2, double L) const
{
  // TO DO: add pbc
  Vec2 f(0,0);
  Vec2 d = r1 - r1;

  double l = sqrt(d.x * d.x + d.y * d.y);

  if ( l < sigmaCutOff_ ) {
    f = -epsilon_*pow(sigma_/l, alpha_) * d;
  }
  
  return f;
}




#endif
