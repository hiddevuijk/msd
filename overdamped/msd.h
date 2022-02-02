#ifndef  GUARD_MSD_H
#define  GUARD_MSD_H

#include "vec2.h"
#include "system.h"

#include <fstream>

class MSD {
 public:
  MSD(std::vector<Vec2> positions)
    : N_(positions.size()),initPositions_(positions) {}; 

  void sample(double time, std::vector<Vec2> positions);
  void save(std::string outname) const;

 private:

  int N_;
  std::vector<Vec2> initPositions_;
  std::vector<double> distance_;
  std::vector<double> time_;
};

void MSD::sample(double time, std::vector<Vec2> positions)
{
  time_.push_back(time); 
  distance_.push_back(0);
  for (int i = 0; i < N_; ++i) {
    double d = positions[i].x - initPositions_[i].x;
    distance_.back() +=  d*d/N_;
    d = positions[i].y - initPositions_[i].y;
    distance_.back() +=  d*d/N_;
  }

}

void MSD::save(std::string outname) const
{
  std::ofstream out;
  out.open(outname);
  for (unsigned int i = 0; i < distance_.size(); ++i) {
    out << time_[i] << "\t" << distance_[i];
    if (i < distance_.size() - 1 ) out << "\n";
  }
  out.close();
}
#endif
