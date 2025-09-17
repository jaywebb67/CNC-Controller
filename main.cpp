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
  mp.plan_SingleLinearMotion(60,80,0,0,0,0);
  mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  mp.print_profile();

  mp.update(1000,200,10000);
  mp.plan_SingleLinearMotion(100,0,0,0,0,0);
  mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  mp.print_profile();

  // mp.plan_SingleLinearMotion(30,0,0,0,0,0);
  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile();

  return 0;
}