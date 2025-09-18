#include <cstdint>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

struct HalfTimes { double t1, t2; }; // t3 == t1
struct BlockTimes { double t1u,t2u,t3u,t4,t5d,t6d,t7d; };

struct BlockPlan {
  BlockTimes t;         // t1u,t2u,t3u,t4,t5d,t6d,t7d
  double v0, vp, v1;    // entry, peak/cruise, exit speeds
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

        void print_timeSteps(const std::vector<double>& acc,
                            const std::vector<double>& vel,
                            const std::vector<double>& pos,
                            double fs)
        {
            const std::size_t N = std::min({acc.size(), vel.size(), pos.size()});
            const double dt = 1.0 / fs;
            std::cout << "Time Step\tt[s]\tacc\tvel\tpos\n";
            std::cout << std::fixed << std::setprecision(6);
            for (std::size_t i = 0; i < N; ++i) {
                const double t = i * dt;
                std::cout << i << '\t' << t << '\t' << acc[i] << '\t' << vel[i] << '\t' << pos[i] << '\n';
            }
        }

        void print_profile(BlockTimes t) const {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "t1="<<t.t1u<<" t2="<<t.t2u<<" t3="<<t.t3u<<" t4="<<t.t4
                    <<" t5="<<t.t5d<<" t6="<<t.t6d<<" t7="<<t.t7d<<"\n";
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


};