#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

//
// TO TO
//

#include "vec2.h"
#include "potential.h"

#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <boost/random.hpp>

#include <iostream>



template <class Potential>
class System {
 public:
  System(int N, double L, double dt, double rVerlet, Potential potential,
         double T, double gamma, int seed);
   
  // evolve in time
  void integrate(double t);

  void setPosition(int i, double x, double y) {
      positions_[i].x = x; positions_[i].y = y; }

  double getPositionX(int i) const { return positions_[i].x; };
  double getPositionY(int i) const { return positions_[i].y; };

  int getNumberOfParticles() const { return N_; }
  int getSystemSize() const { return L_; }
  std::vector<Vec2> getPositions() const { return positions_; }

  void setKappa(double kappa1, int N1, double kappa2);
  void setT(double T1, int N1, double T2);
  void setGamma(double gamma1, int N1, double gamma2);

  void savePositions(std::string outname) const;


  void backInBox();

 private:
  // set all positions of the particles to 
  void initOnGrid();

  // calculate D from T_ and gamma_
  void resetD();

  // calculate forces
  void setForces();

  // update Verlet list
  void updateVerletList();

  // single time step
  void step(double dt);
  // single time setp without interactions
  void stepFree(double dt);



  int N_;    // number of particles
  double L_; // systems size
  double dt_; // time step
  double rVerlet_; // verlet radius

  std::vector<Vec2> positions_;  
  std::vector<Vec2> forces_;

  Potential potential_;

  std::vector<std::vector<int> > verletList_;

  std::vector<Vec2> positionsAtUpdate_;

  std::vector<double> kappa_;
  std::vector<double> T_;
  std::vector<double> gamma_;
  std::vector<double> D_;

	// random number generator

	int seed_;
	const boost::normal_distribution<double> ndist;
	boost::mt19937 rng;		
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > rndist;
};

template <class Potential>
System<Potential>::System(int N, 
                          double L,
                          double dt,
                          double rVerlet,
                          Potential potential,
                          double T,
                          double gamma,
                          int seed)
    : N_(N),
      L_(L),
      dt_(dt),
      rVerlet_(rVerlet),
      positions_(N),
      forces_(N),
      potential_(potential),
      verletList_(N, std::vector<int>(N, -1)),
      positionsAtUpdate_(N),
      kappa_(N, 0.0),
      T_(N, T),
      gamma_(N, gamma),
      D_(N, 0.0),
      seed_(seed),
      ndist(0., 1.),
      rng(seed),
      rndist(rng, ndist) 
{ 
  initOnGrid();
  updateVerletList();
  resetD();

};


template <class Potential>
void System<Potential>::integrate(double t)
{
  if( potential_.getEpsilon() > 0 ) {
    while (t - dt_ > 0) {
      step(dt_); 
      t -= dt_;
    }
    step(t);
  } else {
    while (t - dt_ > 0) {
      stepFree(dt_); 
      t -= dt_;
    }
    stepFree(t);
  }
}

template <class Potential>
void System<Potential>::stepFree(double dt)
{
  Vec2 dr(0, 0);
  
  for (int i = 0; i < N_; ++i) {
    dr.x = sqrt( 2*D_[i] * dt) * rndist();
    dr.y = sqrt( 2*D_[i] * dt) * rndist();
    positions_[i]  += dr;
  }

}

template <class Potential>
void System<Potential>::step(double dt)
{
  setForces();
  Vec2 dr(0, 0);
  
  bool update = false;
  for (int i = 0; i < N_; ++i) {

    dr.x = -1*forces_[i].x - kappa_[i] * forces_[i].y;
    dr.y = -1*forces_[i].y + kappa_[i] * forces_[i].x;
    dr *= dt/(gamma_[i]*(1+kappa_[i]*kappa_[i]));

    dr.x += sqrt( 2*D_[i] * dt) * rndist();
    dr.y += sqrt( 2*D_[i] * dt) * rndist();

    positions_[i]  += dr;

    // CHECK: TO DO
    dr = positionsAtUpdate_[i] - positions_[i];
    if ( 2 * std::sqrt( dr.x*dr.x + dr.y*dr.y)
                  > rVerlet_ - potential_.getSigmaCutOff() ) {
      update = true;
    }

  }

  if (update) { updateVerletList(); } 

}




template <class Potential>
void System<Potential>::setForces()
{
  for (int i = 0; i < N_; ++i) { forces_[i] *= 0; }
  
  Vec2 f;
  int j;
  for (int i = 0; i < N_; ++i) {
    for (int k = 0; k < i; ++k ) {
      j = verletList_[i][k]; 
      if( j == -1 ) break;
      
      
      f = potential_.force(positions_[i], positions_[j], L_);

      forces_[i] += f;
      forces_[j] -= f;
      
    }
  }

}

template <class Potential>
void System<Potential>::updateVerletList()
{
  positionsAtUpdate_ = positions_; 

  Vec2 d;
  for (int i = 0; i < N_; ++i) {
    std::fill(verletList_[i].begin(), verletList_[i].end(), -1);

    int k = 0;
    for (int j = 0; j < i; ++j) {
      d = positions_[i] - positions_[j]; 
      d.x -= L_ * std::round( d.x / L_);
      d.y -= L_ * std::round( d.y / L_);

      // if neighbour j is in verlet radius, add to list
      if (d.x*d.x + d.y*d.y < rVerlet_*rVerlet_) {
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

  boost::random::uniform_int_distribution<> ridist(0, n*n - 1);

  // dl is the spacing between the grid points
  double dl = L_ / n;

  std::vector<int> indices(n*n);
  for (int i = 0; i < n*n; ++i) indices[i] = i;
  // shuffle indices
  for (int i = 0; i < n*n; ++i) {
    int j = ridist(rng);    
    int temp = indices[i]; 
    indices[i] = indices[j];
    indices[j] = temp;
  }

  for (int i = 0; i < N_; ++i) {
    int I = indices[i]; // index of the node
    int ix = I % n; // x index of the node
    int iy = (I - ix) / n; // y index of the node

    positions_[i].x = ix * dl;
    positions_[i].y = iy * dl;
  }

}

template <class Potential>
void System<Potential>::setKappa(double kappa1, int N1, double kappa2)
{
  for (int i = 0; i < N_; ++i) {
    if (i < N1) {
      kappa_[i] = kappa1;
    } else {
      kappa_[i] = kappa2;
    }
  }

}

template <class Potential>
void System<Potential>::setT(double T1, int N1, double T2)
{
  for (int i = 0; i < N_; ++i) {
    if (i < N1) {
      T_[i] = T1;
    } else {
      T_[i] = T2;
    }
  }
  resetD();
}

template <class Potential>
void System<Potential>::setGamma(double gamma1, int N1, double gamma2)
{
  for (int i = 0; i < N_; ++i) {
    if (i < N1) {
      gamma_[i] = gamma1;
    } else {
      gamma_[i] = gamma2;
    }
  }
  resetD();
}

template <class Potential>
void System<Potential>::resetD()
{
  for (int i = 0; i < N_; ++i) {
    D_[i] = T_[i]/gamma_[i];
  }
}

template <class Potential>
void System<Potential>::savePositions(std::string outname) const
{

	std::ofstream out;
	out.open(outname);
  for (int i = 0; i < N_; ++i) {
    out << positions_[i].x << "\t"
        << positions_[i].y;
    if (i < N_ - 1) out << "\n";
  }

	out.close();
}

template <class Potential>
void System<Potential>::backInBox()
{
  for (int i = 0; i < N_; ++i) {
    positions_[i].x -= L_ * std::floor( positions_[i].x / L_ );
    positions_[i].y -= L_ * std::floor( positions_[i].y / L_ );
  }

}


#endif //GUARD_SYSTEM_H
