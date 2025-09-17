#include "s-curve.hpp"

int main() {
  motion_planner mp(1000,10000,100); // Amax=1000, Jmax=10000, Vmax=100

  // // L = 10 mm
  // mp.plan_SingleLinearMotion(10,0,0,0,0,0);
  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile();

  // // L = 20 mm
  // mp.plan_SingleLinearMotion(20,0,0,0,0,0);
  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile();

  // // L = 40 mm
  // mp.plan_SingleLinearMotion(40,0,0,0,0,0);
  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile();

  // L = 100 mm
  auto plan = mp.plan_single_linear_asym(1, /*v0=*/0.0, /*v1=*/0.0,/*vp=*/100,
                                      /*Jup=*/10000, /*Aup=*/1000,
                                      /*Jdn=*/10000, /*Adn=*/1000);

  mp.sample_profile_asym(2000.0, plan.t, plan.v0, plan.vp, plan.v1,
                    10000, 1000, 10000, 1000,
                    mp.acc, mp.vel, mp.pos);

  mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  mp.print_profile(plan.t);
  // mp.update(1000,200,10000);
  // mp.plan_SingleLinearMotion(100,0,0,0,0,0);
  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile();

  // mp.plan_SingleLinearMotion(30,0,0,0,0,0);
  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile();

  return 0;
}