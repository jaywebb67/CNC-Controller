#include <cstdint>
#include <iostream>
#include <iomanip>
#include <vector>
#include "utilities.cpp"
#include <fstream>
#include <string>
#include <limits> // you already use std::numeric_limits



struct HalfTimes { double t1, t2; }; // t3 == t1
struct BlockTimes { double t1u,t2u,t3u,t4,t5d,t6d,t7d; };

struct BlockPlan {
  BlockTimes t;         // t1u,t2u,t3u,t4,t5d,t6d,t7d
  double v0, vp, v1;    // entry, peak/cruise, exit speeds
};

// Results for the whole polyline
struct MultiBlockPlan {
  std::vector<double> L_raw, L_trim;      // lengths before/after trimming
  std::vector<double> v_entry, v_exit;    // per-block entry/exit speeds
  std::vector<double> R_junc, vcap_junc;  // per junction (size N-1)
  std::vector<double> d_left, d_right;    // trim applied to each block
  std::vector<BlockPlan> blocks;          // full S-curve plan per block
};

class motion_planner {
    private:
        double accel_max, Jup, deccel_max, Jdn;
    
    
        public:
        std::vector<double> acc, vel, pos;


        motion_planner(double accel, double jerkU, double deccel, double jerkD)
        : accel_max(accel), Jup(jerkU), deccel_max(deccel), Jdn(jerkD) {}

        BlockPlan plan_single_linear_asym(double L, double v0, double v1,
                                        double Vmax,double Aup,double Adn)
        {
            BlockPlan bp{};
            bp.v0 = v0; 
            bp.v1 = v1;

            // Peak speed bracket: vp ∈ [lo, hi]
            double lo = std::max(v0, v1);
            double hi = std::max(Vmax, lo);   // ensure hi >= lo

            // Try cruise at hi
            auto up_try    = solve_half(std::max(0.0, hi - v0), Jup, Aup);
            auto s_up_try  = s_half(v0, up_try, Jup, Aup);

            auto dn_try    = solve_half(std::max(0.0, hi - v1), Jdn, Adn);
            auto s_dn_try  = s_half(v1, dn_try, Jdn, Adn);   // <-- FIX: v1, not hi/vp

            double need = s_up_try + s_dn_try;
            if (L >= need) {
                bp.t.t1u = up_try.t1; bp.t.t2u = up_try.t2; bp.t.t3u = up_try.t1;
                bp.t.t5d = dn_try.t1; bp.t.t6d = dn_try.t2; bp.t.t7d = dn_try.t1;
                bp.t.t4  = (hi > 0.0) ? (L - need) / hi : 0.0;
                if (bp.t.t4 < 0) bp.t.t4 = 0.0;
                bp.vp = hi;
                return bp;
            }

            // No-cruise: solve vp by bisection
            for (int it = 0; it < 60; ++it) {
                double vp    = 0.5 * (lo + hi);
                auto upM     = solve_half(std::max(0.0, vp - v0), Jup, Aup);
                auto s_upM   = s_half(v0, upM, Jup, Aup);

                auto dnM     = solve_half(std::max(0.0, vp - v1), Jdn, Adn);
                auto s_dnM   = s_half(v1, dnM, Jdn, Adn);     // <-- FIX: v1 here too

                double f = s_upM + s_dnM - L;
                if (f > 0) hi = vp; else lo = vp;
                if (std::abs(f) < 1e-12) { lo = hi = vp; break; }
            }
            bp.vp = 0.5 * (lo + hi);

            auto upF = solve_half(std::max(0.0, bp.vp - v0), Jup, Aup);
            auto dnF = solve_half(std::max(0.0, bp.vp - v1), Jdn, Adn);

            bp.t.t1u = upF.t1; bp.t.t2u = upF.t2; bp.t.t3u = upF.t1;
            bp.t.t4  = 0.0;
            bp.t.t5d = dnF.t1; bp.t.t6d = dnF.t2; bp.t.t7d = dnF.t1;

            // Clamp tiny negatives to zero
            auto clamp0 = [](double& x){ if (x < 1e-12) x = 0.0; };
            clamp0(bp.t.t1u); clamp0(bp.t.t2u); clamp0(bp.t.t3u);
            clamp0(bp.t.t4);  clamp0(bp.t.t5d); clamp0(bp.t.t6d); clamp0(bp.t.t7d);

            return bp;
        }


        BlockPlan plan_single_linear_asym(double L, double v0, double v1,
                                        double Vmax) {
            return plan_single_linear_asym(L, v0, v1, Vmax, accel_max, deccel_max);
        }
        

            // Look-ahead over a polyline: pts.size()==N+1, feeds.size()==N
            // an_max: centripetal accel cap for junctions; eps: blend tolerance
            MultiBlockPlan plan_path_with_lookahead(const std::vector<Vec3>& pts,
                                                    const std::vector<double>& feeds,
                                                    double eps, double an_max,
                                                    double Jup, double Aup,
                                                    double Jdn, double Adn,
                                                    double v_start, double v_end)
            {
                const std::size_t N = feeds.size(); // # of linear blocks
                MultiBlockPlan out{};
                out.L_raw.resize(N);
                out.L_trim.assign(N, 0.0);
                out.v_entry.assign(N, 0.0);
                out.v_exit.assign(N, 0.0);
                out.blocks.resize(N);

                // geometry
                std::vector<Vec3> dir(N);
                for (std::size_t i=0;i<N;++i) {
                    Vec3 d = pts[i+1]-pts[i];
                    out.L_raw[i] = norm(d);
                    dir[i] = unit(d);
                }

                // junctions (N-1 of them)
                out.R_junc.assign(N?N-1:0, 0.0);
                out.vcap_junc.assign(N?N-1:0, std::numeric_limits<double>::infinity());

                // trims applied to each block’s entry/exit
                out.d_left.assign(N, 0.0);
                out.d_right.assign(N, 0.0);

                const double ANG_EPS = 1e-7;

                for (std::size_t i=0; i+1<N; ++i) {
                    const double c = clamp01(dot(dir[i], dir[i+1]));
                    const double theta = std::acos(c);            // turn angle
                    // Straight-ish: no blend, no cap
                    if (theta < 1e-6) { out.R_junc[i]=INFINITY; out.vcap_junc[i]=feeds[i+1]; continue; }
                    // Near 180° corner: force stop
                    if (theta > M_PI - ANG_EPS) { out.R_junc[i]=0.0; out.vcap_junc[i]=0.0; continue; }

                    const double half = 0.5*theta;
                    const double cosH = std::cos(half);
                    const double tanH = std::tan(half);

                    // Ideal radius from tolerance
                    double R = eps / std::max(1e-12, (1.0 - cosH));
                    // Limit radius so trims fit inside adjacent blocks
                    const double L_i   = out.L_raw[i];
                    const double L_ip1 = out.L_raw[i+1];
                    const double R_by_len = std::min(0.5*L_i, 0.5*L_ip1) / std::max(1e-12, tanH);
                    if (std::isfinite(R)) R = std::min(R, R_by_len); else R = R_by_len;

                    const double d = R * tanH;              // trim on each adjacent line
                    out.d_right[i]   = std::min(0.5*L_i,   d);
                    out.d_left[i+1]  = std::min(0.5*L_ip1, d);

                    out.R_junc[i] = R;
                    out.vcap_junc[i] = std::sqrt(std::max(0.0, an_max * R));
                    // also cap by nominal feeds on either side
                    out.vcap_junc[i] = std::min(out.vcap_junc[i], std::min(feeds[i], feeds[i+1]));
                }

                // Trimmed lengths
                for (std::size_t i=0;i<N;++i) {
                    double Ltr = out.L_raw[i] - out.d_left[i] - out.d_right[i];
                    out.L_trim[i] = std::max(1e-9, Ltr); // never zero/negative
                }

                // Per-block entry/exit caps from junction caps and feeds
                auto entry_cap = [&](std::size_t i){
                    double c = feeds[i];
                    if (i>0)  c = std::min(c, out.vcap_junc[i-1]);
                    if (i==0) c = std::min(c, v_start);
                    return c;
                };
                auto exit_cap = [&](std::size_t i){
                    double c = feeds[i];
                    if (i+1<N) c = std::min(c, out.vcap_junc[i]);
                    if (i+1==N) c = std::min(c, v_end);
                    return c;
                };

                // Forward pass
                out.v_entry[0] = std::min(v_start, entry_cap(0));
                for (std::size_t i=0;i<N;++i) {
                    double v_cap_exit = exit_cap(i);
                    double L          = out.L_trim[i];
                    double v_in       = out.v_entry[i];

                    out.v_exit[i] = solve_vexit_forward(v_in, L, v_cap_exit, Jup, Aup, Jdn, Adn);
                    if (i+1<N) {
                        // next entry is current exit, but obey its entry cap
                        out.v_entry[i+1] = std::min(out.v_exit[i], entry_cap(i+1));
                    }
                }

                // Backward pass
                if (N) out.v_exit[N-1] = std::min(out.v_exit[N-1], v_end);
                for (int i=int(N)-1; i>=0; --i) {
                    const double L = out.L_trim[std::size_t(i)];
                    const double v_in  = out.v_entry[std::size_t(i)];
                    const double v_out = out.v_exit[std::size_t(i)];

                    if (Dmin_change(v_in, v_out, Jup, Aup, Jdn, Adn) > L + 1e-12) {
                        // reduce entry to the largest value that still fits
                        double cap = entry_cap(std::size_t(i));
                        double v_in_new = solve_ventry_backward(v_out, L, std::min(cap, v_in), Jup, Aup, Jdn, Adn);
                        out.v_entry[std::size_t(i)] = v_in_new;
                        if (i-1 >= 0) out.v_exit[std::size_t(i-1)] = std::min(out.v_exit[std::size_t(i-1)], v_in_new);
                    }
                }

                // Final per-block S-curve plans (on trimmed lengths)
                for (std::size_t i=0;i<N;++i) {
                    out.blocks[i] = plan_single_linear_asym(out.L_trim[i],
                                                            out.v_entry[i], out.v_exit[i],
                                                            /*Vmax*/ feeds[i],
                                                            /*up*/Aup, /*dn*/Adn);
                }
                return out;
            }


        void sample_profile_asym(double fs,const BlockTimes& bt,
            double v0, double vp, double v1,double Aup, double Adn,
            std::vector<double>& a_out,std::vector<double>& v_out,
            std::vector<double>& s_out)
        {
            const double t1u=bt.t1u, t2u=bt.t2u, t3u=bt.t3u;
            const double t4 =bt.t4;
            const double t5d=bt.t5d, t6d=bt.t6d, t7d=bt.t7d;

            const double dt = 1.0 / fs;
            const double T  = t1u+t2u+t3u+t4+t5d+t6d+t7d;
            if (T <= 0) { a_out.clear(); v_out.clear(); s_out.clear(); return; }

            const std::size_t N = static_cast<std::size_t>(std::floor(T/dt)) + 1;
            a_out.assign(N, 0.0); v_out.assign(N, 0.0); s_out.assign(N, 0.0);

            // Use A_eff for triangular halves
            const double Aup_eff = (t2u > 0) ? Aup : (Jup * t1u);
            const double Adn_eff = (t6d > 0) ? Adn : (Jdn * t5d);

            // Boundaries
            const double b1=t1u, b2=b1+t2u, b3=b2+t3u, b4=b3+t4, b5=b4+t5d, b6=b5+t6d, b7=b6+t7d;

            // Acc half end states (with Aup_eff)
            const double v1e = v0 + 0.5*Jup*t1u*t1u;
            const double s1e = v0*t1u + (1.0/6.0)*Jup*t1u*t1u*t1u;

            const double v2e = v1e + Aup_eff*t2u;
            const double s2e = s1e + v1e*t2u + 0.5*Aup_eff*t2u*t2u;

            const double v3e = v2e + Aup_eff*t3u - 0.5*Jup*t3u*t3u; // == vp
            const double s3e = s2e + v2e*t3u + 0.5*Aup_eff*t3u*t3u - (1.0/6.0)*Jup*t3u*t3u*t3u;

            // Cruise
            const double s4e = s3e + vp*t4;

            // Dec half end states (with Adn_eff)
            const double v5e = vp - 0.5*Jdn*t5d*t5d;
            const double s5e = s4e + vp*t5d - (1.0/6.0)*Jdn*t5d*t5d*t5d;

            const double v6e = v5e - Adn_eff*t6d;
            const double s6e = s5e + v5e*t6d - 0.5*Adn_eff*t6d*t6d;

            for (std::size_t i=0; i<N; ++i) {
                const double t = i*dt;
                double a, v, s;

                if (t < b1) { // seg1 +Jup
                    const double tau=t;
                    a =  Jup*tau;
                    v =  v0 + 0.5*Jup*tau*tau;
                    s =  v0*tau + (1.0/6.0)*Jup*tau*tau*tau;
                }
                else if (t < b2) { // seg2 +Aup_eff
                    const double tau=t-b1;
                    a =  Aup_eff;
                    v =  v1e + Aup_eff*tau;
                    s =  s1e + v1e*tau + 0.5*Aup_eff*tau*tau;
                }
                else if (t < b3) { // seg3 -Jup from Aup_eff
                    const double tau=t-b2;
                    a =  Aup_eff - Jup*tau;
                    v =  v2e + Aup_eff*tau - 0.5*Jup*tau*tau;
                    s =  s2e + v2e*tau + 0.5*Aup_eff*tau*tau - (1.0/6.0)*Jup*tau*tau*tau;
                }
                else if (t < b4) { // cruise
                    const double tau=t-b3;
                    a = 0.0; v = vp; s = s3e + vp*tau;
                }
                else if (t < b5) { // seg5 -Jdn
                    const double tau=t-b4;
                    a = -Jdn*tau;
                    v =  vp - 0.5*Jdn*tau*tau;
                    s =  s4e + vp*tau - (1.0/6.0)*Jdn*tau*tau*tau;
                }
                else if (t < b6) { // seg6 -Adn_eff
                    const double tau=t-b5;
                    a = -Adn_eff;
                    v =  v5e - Adn_eff*tau;
                    s =  s5e + v5e*tau - 0.5*Adn_eff*tau*tau;
                }
                else {             // seg7 +Jdn back to zero accel
                    const double tau=t-b6;
                    a = -Adn_eff + Jdn*tau;
                    v =  v6e - Adn_eff*tau + 0.5*Jdn*tau*tau;
                    s =  s6e + v6e*tau - 0.5*Adn_eff*tau*tau + (1.0/6.0)*Jdn*tau*tau*tau;
                }

                a_out[i]=a; v_out[i]=v; s_out[i]=s;
            }
        }

        void sample_profile_asym(double fs,const BlockTimes& bt,
            double v0, double vp, double v1,std::vector<double>& a_out,
            std::vector<double>& v_out,std::vector<double>& s_out)
        {
            return sample_profile_asym(fs,bt,v0,vp,v1,accel_max,deccel_max,a_out,v_out,s_out);
        }

        // Save the block times (t1..t7) to CSV, while still printing to stdout.
        void print_profile(const BlockTimes& t, const std::string& csv_path = "") const {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "t1="<<t.t1u<<" t2="<<t.t2u<<" t3="<<t.t3u<<" t4="<<t.t4
                    <<" t5="<<t.t5d<<" t6="<<t.t6d<<" t7="<<t.t7d<<"\n";

            if (!csv_path.empty()) {
                std::ofstream ofs(csv_path, std::ios::trunc);
                if (!ofs) {
                    std::cerr << "Failed to open CSV file: " << csv_path << "\n";
                    return;
                }
                ofs << std::fixed << std::setprecision(9);
                ofs << "t1,t2,t3,t4,t5,t6,t7\n";
                ofs << t.t1u << ',' << t.t2u << ',' << t.t3u << ','
                    << t.t4  << ',' << t.t5d << ',' << t.t6d << ',' << t.t7d << '\n';
            }
        }

        // Save the sampled profile (time series) to CSV for charting in Excel.
        void save_timeSteps_csv(const std::string& csv_path,
                                const std::vector<double>& acc,
                                const std::vector<double>& vel,
                                const std::vector<double>& pos,
                                double fs) const
        {
            const std::size_t N = std::min({acc.size(), vel.size(), pos.size()});
            std::ofstream ofs(csv_path, std::ios::trunc);
            if (!ofs) {
                std::cerr << "Failed to open CSV file: " << csv_path << "\n";
                return;
            }
            ofs << std::fixed << std::setprecision(9);

            // Header
            ofs << "index,time_s,acc,vel,pos\n";

            const double dt = (fs > 0.0) ? 1.0 / fs : 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                const double t = i * dt;
                ofs << i << ',' << t << ',' << acc[i] << ',' << vel[i] << ',' << pos[i] << '\n';
            }
        }


        void update(double newAccel, double newJup, double newJdn,double newDeccel){

            accel_max = newAccel;
            Jup = newJup;
            Jdn = newJdn;
            deccel_max = newDeccel;

        }

    private:

        inline HalfTimes solve_half(double dv, double J, double A) {
            HalfTimes h{};
            double t1 = A/J;
            if (dv >= A*t1) {             // with plateau
                h.t1 = t1;
                h.t2 = dv/A - t1;
            } else {                      // triangular: reduce peak accel
                h.t1 = std::sqrt(dv / J);
                h.t2 = 0.0;
            }
            return h;
        }

        inline double s_half(double v_start, const HalfTimes& h, double J, double A) {
            // If triangular, A means J*h.t1
            const double Aeff = (h.t2>0) ? A : (J*h.t1);
            return v_start * (2*h.t1 + h.t2)
                + (Aeff*h.t1*h.t1 + 1.5*Aeff*h.t1*h.t2 + 0.5*Aeff*h.t2*h.t2);
        }


            // Distance needed to change speed vA -> vB with zero accel endpoints
        // using jerk/accel limits (up vs down).
        inline double Dmin_change(double vA, double vB,
                                double Jup, double Aup,
                                double Jdn, double Adn)
        {
            if (vB >= vA) {
                const double dv = vB - vA;
                auto h = solve_half(dv, Jup, Aup);
                return s_half(vA, h, Jup, Aup);
            } else {
                const double dv = vA - vB;
                auto h = solve_half(dv, Jdn, Adn);
                // time-reversal symmetry: distance equals accel from vB -> vA
                return s_half(vB, h, Jdn, Adn);
            }
        }

        // Solve for maximal reachable v_exit given (v_entry, L, cap_hi),
        // using bisection on a monotone Dmin(v_entry->x)
        inline double solve_vexit_forward(double v_entry, double L,
                                        double cap_hi,
                                        double Jup, double Aup,
                                        double Jdn, double Adn)
        {
            // If we must decel (cap below entry), we'll try to decel as much as we can
            if (cap_hi < v_entry) {
                // find minimal v_exit >= 0 s.t. Dmin(v_entry->v_exit)==L
                double lo = 0.0, hi = v_entry;               // decel range
                // If we can reach the cap, clamp to it; else the solver returns the
                // largest decel possible within L (higher than cap), to be fixed in backward pass.
                if (Dmin_change(v_entry, cap_hi, Jup, Aup, Jdn, Adn) <= L)
                    return cap_hi;
                for (int it=0; it<50; ++it) {
                    double mid = 0.5*(lo+hi);
                    double D = Dmin_change(v_entry, mid, Jup, Aup, Jdn, Adn);
                    (D > L ? lo : hi) = mid;
                }
                return 0.5*(lo+hi);
            }

            // Accelerate up toward cap_hi
            if (Dmin_change(v_entry, cap_hi, Jup, Aup, Jdn, Adn) <= L)
                return cap_hi;

            double lo = v_entry, hi = cap_hi;
            for (int it=0; it<50; ++it) {
                double mid = 0.5*(lo+hi);
                double D = Dmin_change(v_entry, mid, Jup, Aup, Jdn, Adn);
                (D > L ? hi : lo) = mid;
            }
            return 0.5*(lo+hi);
        }

        // Backward: given fixed v_exit and L, find the largest v_entry allowed.
        inline double solve_ventry_backward(double v_exit, double L, double cap_hi,
                                            double Jup, double Aup,
                                            double Jdn, double Adn)
        {
            // If even from cap_hi we can decel to v_exit within L, keep cap
            if (Dmin_change(cap_hi, v_exit, Jup, Aup, Jdn, Adn) <= L)
                return cap_hi;

            double lo = 0.0, hi = cap_hi;
            for (int it=0; it<50; ++it) {
                double mid = 0.5*(lo+hi);
                double D = Dmin_change(mid, v_exit, Jup, Aup, Jdn, Adn);
                (D > L ? hi : lo) = mid;
            }
            return 0.5*(lo+hi);
        }

};