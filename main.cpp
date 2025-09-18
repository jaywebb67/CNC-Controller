#include "s-curve.hpp"

int main() {
  motion_planner mp(1000,10000,1000,10000); // Amax=1000, Jmax=10000, Vmax=100

  // // // L = 1 mm
  // auto plan = mp.plan_single_linear_asym(100, /*v0=*/0.0, /*v1=*/0.0,/*vp=*/100);

  // mp.sample_profile_asym(2000.0, plan.t, plan.v0, plan.vp, plan.v1,
  //                   mp.acc, mp.vel, mp.pos);

  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile(plan.t);

  std::vector<Vec3> pts = { {0,0,0}, {100,0,0}, {200,0,0}, {300,0,0} };
  std::vector<double> feeds = { 100, 100, 100 }; // mm/s caps per block

  double eps   = 0.1;    // mm path tolerance
  double anmax = 2000.0; // mm/s^2 centripetal cap at blends
  double v0=0, vN=0;

  auto mbp = mp.plan_path_with_lookahead(pts, feeds, eps, anmax,
                                        /*Jup,Aup*/ 10000,1000,
                                        /*Jdn,Adn*/ 10000,1000,
                                        v0, vN);

  for (size_t i = 0; i < mbp.blocks.size(); ++i) {
    const auto& b = mbp.blocks[i];
    mp.sample_profile_asym(2000.0, b.t, b.v0, b.vp, b.v1, mp.acc, mp.vel, mp.pos);
    mp.print_profile(b.t, "block" + std::to_string(i) + "_times.csv");
    mp.save_timeSteps_csv("block" + std::to_string(i) + "_timeseries.csv", mp.acc, mp.vel, mp.pos, 2000.0);
  }


  return 0;
}