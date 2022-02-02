
#include "system.h"
#include "msd.h"
#include "diffusion.h"
#include "multi_diffusion.h"

#include "ConfigFile.h"

#include <iostream>
#include <vector>

using namespace std;

void sampleAll(const System<WCApotential>& sys, MDiffusion& md)
{
  for (int i = 0; i < sys.getNumberOfParticles(); ++i) {
    md.sample(2 * i    , sys.getPositionX(i));
    md.sample(2 * i + 1, sys.getPositionY(i));
  }
}

int main() {

  ConfigFile config("input.txt");
  int N          = config.read<int>("N");
  double L       = config.read<double>("L");
  double dt      = config.read<double>("dt");
  double T1      = config.read<double>("T1");
  double T2      = config.read<double>("T2");
  double gamma1  = config.read<double>("gamma1");
  double gamma2  = config.read<double>("gamma2");
  int seed       = config.read<int>("seed");
  double rVerlet = config.read<double>("rVerlet");
  int N1         = config.read<int>("N1");
  double kappa1  = config.read<double>("kappa1");
  double kappa2  = config.read<double>("kappa2");

  double sigma   = config.read<double>("sigma");
  double sigmaCO = config.read<double>("sigmaCO");
  double epsilon = config.read<double>("epsilon");
  double alpha   = config.read<double>("alpha");
 
  double equilibrationTime = config.read<double>("equilibrationTime");
  double NsampleTime       = config.read<double>("NsampleTime");
  double sampleDTime       = config.read<double>("sampleDTime");
  double totalTime         = config.read<double>("totalTime");

  WCApotential potential(sigma, sigmaCO, epsilon, alpha);

  System<WCApotential> system(N, L, dt, rVerlet, potential,T1,gamma1, seed);

  system.setKappa(kappa1, N1, kappa2);
  system.setT(T1, N1, T2);
  system.setGamma(gamma1, N1, gamma2);

  system.savePositions("r0.dat");

  double time = 0;

  system.integrate(equilibrationTime);

  //MSD msd(system.getPositions());
  MDiffusion diff(2*N, NsampleTime, sampleDTime);

  double Temperature = 0;
  int i = 0;

  while (time < totalTime) {
    if( i % 100 == 0  and i > 0){
    cout << time << endl;
    }
    i += 1;
    
    system.integrate(sampleDTime);
    sampleAll(system, diff); 
    //msd.sample(time, system.getPositions());
    time += sampleDTime;
  }
  cout << Temperature/(i*N) << endl;
  system.backInBox();
  system.savePositions("r.dat");
  diff.save("msd.dat");  
  //msd.save("msd.dat");

  return 0;
}
