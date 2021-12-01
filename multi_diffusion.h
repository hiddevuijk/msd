#ifndef GUARD_MULTI_DIFFUSION_H
#define GUARD_MULTI_DIFFUSION_H

#include "diffusion.h"

#include <vector>
#include <string>
#include <fstream>

class MDiffusion {
 public:
  MDiffusion(int N, int ntime, double dt, int sampleEvery = 1);

  void sample(int i, double x);

  std::vector<double> getR2T(int i) const;

  void save(std::string outname);
 private:
  int N_;
  int ntime_;
  std::vector<Diffusion> diffs_; 
};


MDiffusion::MDiffusion(int N, int ntime, double dt, int sampleEvery)
  : N_(N), ntime_(ntime), diffs_(N, Diffusion(ntime, dt, sampleEvery))
{ }

void MDiffusion::sample(int i, double x)
{
  diffs_[i].sample(x);
}

std::vector<double> MDiffusion::getR2T(int i) const
{
  return diffs_[i].getR2T();
}

void MDiffusion::save(std::string outname) 
{
  std::ofstream out;
  out.open(outname);

  
  for (int ti = 0; ti < ntime_; ++ti) {
    out << diffs_[0].time_[ti];
    if (ti < ntime_ - 1) out << "\t";
  }
  for (int i = 0; i < N_; ++i) {
    out << "\n";
    std::vector<double> r2t = diffs_[i].getR2T();    
    for (int ti = 0; ti < ntime_; ++ti) {
      out << r2t[ti];
      if (ti < ntime_ - 1) out << "\t";
    }
  }
}

#endif
