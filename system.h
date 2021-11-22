#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include "vec2.h"
#include "potential.h"

#include <vector>
#include <algorithm>

#include <iostream>

template <class Potential>
class System {
 public:
  System(int N, double L, double dt, double rVerlet, Potential potential)
    : N_(N), L_(L), dt_(dt), rVerlet_(rVerlet), positions_(N), velocities_(N_),
      forces_(N), potential_(potential),
      verletList_(N, std::vector<int>(N, -1)),
      distanceSinceUpdate_(N)
    { 
      initOnGrid();
      resetVelocities();
      updateVerletList();
    };

  // evolve in time
  void integrate(double t);

  int getNumberOfParticles() const { return N_; }
  int getSystemSize() const { return L_; }

 private:
  // set all positions of the particles to 
  void initOnGrid();

  // set all velocities to 0
  void resetVelocities();

  // calculate forces
  void setForces();

  // update Verlet list
  void updateVerletList();

  // single time step
  void step(double dt);

  // single step without magnetic field
  void initStep(double dt);


  int N_;    // number of particles
  double L_; // systems size
  double dt_; // time step
  double rVerlet_; // verlet radius

  std::vector<Vec2> positions_;  
  std::vector<Vec2> velocities_;
  std::vector<Vec2> forces_;

  Potential potential_;

  std::vector<std::vector<int> > verletList_;

  std::vector<double> distanceSinceUpdate_;

  double gamma_;
};

template <class Potential>
void System<Potential>::integrate(double t)
{
  while (t - dt_ > 0) {
    step(dt_); 
    t -= dt_;
  }
  step(t - dt_);
}

template <class Potential>
void System<Potential>::step( double dt)
{
  setForces();
  Vec2 d(0, 0);

  double maxD = 0;  
  for (int i = 0; i < N_; ++i) {
    d = velocities_[i] * dt;

    velocities_[i] += forces_[i] * dt / gamma_;
  
    positions_[i] += d;

    // CHECK: TO DO
    distanceSinceUpdate_[i] += std::sqrt(d.x * d.x + d.y * d.y);
  }

}

template <class Potential>
void System<Potential>::initStep( double dt)
{
  setForces();
  //TO DO
}



template <class Potential>
void System<Potential>::setForces()
{
  // check if verletList needs to be updated
  double maxD = 0;
  for (int i = 0; i < N_; ++i) {
    if ( distanceSinceUpdate_[i] > maxD ) maxD = distanceSinceUpdate_[i];
  }

 if( 2 * maxD > rVerlet_ - potential_.getSigmaCutOff() ) {
    updateVerletList();
 }

  Vec2 f;
  int k;
  for (int i = 0; i < N_; ++i) {
    for (int j = 0; j < N_; ++j ) {
      k = verletList_[i][j];

      if (k == -1) break;

      // add PBC !!!!!!!
      f = potential_.force(positions_[i], positions_[k] );

      forces_[i] += f;
      forces_[k] -= f;
      
    }
  }

}

template <class Potential>
void System<Potential>::updateVerletList()
{
  Vec2 d;
  for (int i = 0; i < N_; ++i) {
    std::fill(verletList_[i].begin(), verletList_[i].end(), -1);
    int k = 0;
    for (int j = 0; j < i; ++j) {
      d = positions_[i] - positions_[j]; 
      if (d.x*d.x + d.y*d.y < rVerlet_ ) {
        verletList_[i][k] = j;
        ++k;
      }
    }
  } 
}

template <class Potential>
void System<Potential>::initOnGrid()
{

  // n is the number of grid points in one dimension
  int n = std::ceil( std::sqrt(N_) );

  // dl is the spacing between the grid points
  double dl = L_ / n;

  std::vector<int> indices(n*n);
  for( int i = 0; i < n*n; ++i) indices[i] = i;

  // shuffle indices;
  // TO DO

  for (int i = 0; i < N_; ++i) {
    // x and y indices
    // TO DO
    int ix = i % n;
    int iy = (i - ix) / n;

    positions_[i].x = ix * dl;
    positions_[i].y = iy * dl;
  }

}

template <class Potential>
void System<Potential>::resetVelocities()
{
  for (int i = 0; i < N_; ++i) {
    velocities_[i].x = 0;
    velocities_[i].y = 0;
  }
}

#endif //GUARD_SYSTEM_H
