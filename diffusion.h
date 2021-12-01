#ifndef GUARD_DIFFUSION_H
#define GUARD_DIFFUSION_H

#include <vector>
#include <string>

class MDiffusion;

class Diffusion {
 friend class MDiffusion;
 public:
  Diffusion(int ntime, double dt, int sampleEvery = 1);

  void sample(double x);

  std::vector<double> getR2T() const;

  void save(std::string outname) const;

 private:
  void calculateR2T();

  int sampleEvery_;
  int sampleInterval_;
  int nsamp_;
  int ntime_;
  int nfill_;
  int lastAdded_;
  std::vector<double> time_;
  std::vector<double> r2t_;

  std::vector<double> x_;
};

Diffusion::Diffusion(int ntime, double dt, int sampleEvery)
  : sampleEvery_(sampleEvery), sampleInterval_(0), nsamp_(0), ntime_(ntime),
    nfill_(0), lastAdded_(0), time_(ntime), r2t_(ntime, 0), x_(ntime)
{
  for(int i = 0; i < ntime_; ++i) {
    time_[i] = i * dt;
  }
};

void Diffusion::sample(double x)
{
  if(nfill_ < ntime_ - 1) {
    x_[nfill_] = x;
    lastAdded_ = nfill_;
    ++nfill_;
  } else {
    lastAdded_ = (lastAdded_ + 1 ) % ntime_;
    x_[lastAdded_] = x;
    ++sampleInterval_;
    if (sampleInterval_ == sampleEvery_) {
      sampleInterval_ = 0;
      calculateR2T();
    }
  }
}


void Diffusion::calculateR2T()
{

  ++nsamp_;

  int j;
  double d;
  for(int i = 0; i < ntime_; ++i) {
    j = lastAdded_ - i;
    if( j < 0 ) j += ntime_;
    d = x_[lastAdded_] - x_[j];
    r2t_[i] += d*d;
  }
  
}

std::vector<double> Diffusion::getR2T() const
{
  std::vector<double> r2tAverage(ntime_, 0);
  for(int i=0; i < ntime_; ++i) r2tAverage[i] = r2t_[i] / nsamp_;
  return r2tAverage;
}

void Diffusion::save(std::string outname) const
{
  std::ofstream out;
  out.open(outname);
  for (int i = 0; i < ntime_; ++i) {
    out << time_[i] << "\t" << r2t_[i]/nsamp_;
    if (i < ntime_ - 1 ) out << "\n";
  }
  out.close();
}

#endif
