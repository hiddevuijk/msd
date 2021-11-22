
#include "system.h"

#include <iostream>
#include <vector>

using namespace std;

int main() {
  int N = 11;
  double L = 10;
  double dt = 1.;
  double rVerlet = 1.;
  WCApotential potential(0, 0, 0);

  System<WCApotential> system(N, L, dt, rVerlet, potential);

  return 0;
}
