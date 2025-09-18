#include "s-curve.hpp"

int main() {
  motion_planner mp(1000,10000,1000,10000); // Amax=1000, Jmax=10000, Vmax=100

  // // L = 1 mm
  auto plan = mp.plan_single_linear_asym(100, /*v0=*/0.0, /*v1=*/0.0,/*vp=*/100);

  mp.sample_profile_asym(2000.0, plan.t, plan.v0, plan.vp, plan.v1,
                    mp.acc, mp.vel, mp.pos);

  mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  mp.print_profile(plan.t);


  return 0;
}