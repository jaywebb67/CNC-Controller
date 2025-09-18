#include "s-curve.hpp"

int main() {
  motion_planner mp(1000,10000,1000,10000); // Amax=1000, Jmax=10000, Vmax=100

  // // // L = 1 mm
  // auto plan = mp.plan_single_linear_asym(100, /*v0=*/0.0, /*v1=*/0.0,/*vp=*/100);

  // mp.sample_profile_asym(2000.0, plan.t, plan.v0, plan.vp, plan.v1,
  //                   mp.acc, mp.vel, mp.pos);

  // mp.print_timeSteps(mp.acc,mp.vel,mp.pos,2000);
  // mp.print_profile(plan.t);

  std::vector<Vec3> pts = { {0,0,0}, {100,0,0}, {150,50,0}, {100,100,0}, {0,100,0}, {0,0,0} };
  std::vector<double> feeds = { 100, 100, 100, 100, 100}; // mm/s caps per block

  double eps   = 0.01;    // mm path tolerance
  double anmax = 2000.0; // mm/s^2 centripetal cap at blends
  double v0=0, vN=0;

  auto mbp = mp.plan_path_with_lookahead(pts, feeds, eps, anmax,
                                        /*Jup,Aup*/ 10000,1000,
                                        /*Jdn,Adn*/ 10000,1000,
                                        v0, vN);

  for (size_t i = 0; i < mbp.blocks.size(); ++i) {
      const auto& b = mbp.blocks[i];
      std::cout << "Block" << std::to_string(i)<< " L(raw)= " << mbp.L_raw[i];
      std::cout << "\tL(trimmed)" << mbp.L_trim[i] << "\n";
      std::cout << "Block" << std::to_string(i)<< " D(Left)= " << mbp.d_left[i];
      std::cout << "\tD(right)" << mbp.d_right[i] << "\n";
      // For each block i
      BlockGeom G = motion_planner::make_block_geom(pts[i], pts[i+1], mbp.d_left[i], mbp.d_right[i]);

      std::vector<double> a,v,s;
      mp.sample_profile_asym(2000.0, b.t, b.v0, b.vp, b.v1, a, v, s);

      // make sure last s is exactly the trimmed length:
      if (!s.empty()) s.back() = G.L_trim;

      // map to XYZ
      std::vector<double> X,Y,Z;
      motion_planner::map_block_positions(G, s, X, Y, Z);

      // save
      mp.save_timeSteps_csv_extended("block" + std::to_string(i) + "_timeseries.csv",
                                  2000.0, G, a, v, s, X, Y, Z);

  }


  return 0;
}