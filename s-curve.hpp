#include <cstdint>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

class motion_planner {
    private:
        double accel_max, jerk_max, max_feedrate;
        // segment times
        double t1=0,t2=0,t3=0,t4=0,t5=0,t6=0,t7=0;
        int case_id = 0;
    public:
        std::vector<double> acc, vel, pos;

    public:
        motion_planner(uint16_t accel, uint16_t jerk, uint16_t feedrate)
        : accel_max(accel), jerk_max(jerk), max_feedrate(feedrate) {}


        void plan_SingleLinearMotion(float X_target, float Y_target, float Z_target,
                        float X_pos,    float Y_pos,    float Z_pos) {
            const double dX = X_target - X_pos;
            const double dY = Y_target - Y_pos;
            const double dZ = Z_target - Z_pos;
            const double L  = std::sqrt(dX*dX + dY*dY + dZ*dZ);
            if (L <= 0.0) { 
                t1=t2=t3=t4=t5=t6=t7=0; 
                case_id=0; 
                return; 
            }

            const double ux = dX/L;
            const double uy = dY/L;
            const double uz = dZ/L;


            t1 = accel_max / jerk_max;
            const double t2_cap = std::max(0.0, max_feedrate/accel_max - t1);
            const double s_ad   = 2*(accel_max*t1*t1 + 1.5*accel_max*t1*t2_cap + 0.5*accel_max*t2_cap*t2_cap);
            constexpr double EPS = 1e-6;

            if (std::abs(L - s_ad) <= EPS) {
                t2=t2_cap; t3=t1; t4=0; t5=t1; t6=t2_cap; t7=t1;
                case_id = (t2_cap==0.0) ? 4 : 3; // edge collapses to case 4 if no plateau
                //return;
            }
            else if (L > s_ad) {
                const double v_acc = accel_max*(t1+t2_cap);
                t2=t2_cap; t3=t1; t4=(L-s_ad)/v_acc; t5=t1; t6=t2_cap; t7=t1;
                case_id = (t2_cap==0.0) ? 1 : 2;
                //return;
            }
            else {
                const double a=1.0, b=3.0*t1, c=2.0*t1*t1 - L/accel_max;
                double t2_sol = posRootQuad(a,b,c);

                if (t2_sol >= -EPS) {                 // got a non-negative solution (within tol)
                    t2_sol = std::max(0.0, t2_sol);   // clamp tiny negatives to 0
                    const double v_acc = accel_max * (t1 + t2_sol);
                    if (t2_sol <= t2_cap + EPS && v_acc <= max_feedrate + 1e-12) {
                        if (t2_sol <= EPS) {
                            // Case 4: no cruise, t2 == 0 (edge for this L)
                            t2=0; t3=t1; t4=0; t5=t1; t6=0; t7=t1;
                            case_id = 4;
                            //return;
                        } else {
                            // Case 3: no cruise, t2 > 0
                            t2=t2_sol; t3=t1; t4=0; t5=t1; t6=t2_sol; t7=t1;
                            case_id = 3;
                            //return;
                        }
                    }
                }

                // Case 5: no-cruise, reduce t1 (pure jerk-bounded)
                const double t1p = std::cbrt(L / (2.0 * jerk_max));
                t1=t1p; t2=0; t3=t1p; t4=0; t5=t1p; t6=0; t7=t1p;
                case_id = 5;
            }
            sample_profile(2000,acc,vel,pos,t1,t2,t3,t4,t5,t6,t7,accel_max,jerk_max);
        }

        

        void sample_profile(double fs,
                    std::vector<double>& a_out,
                    std::vector<double>& v_out,
                    std::vector<double>& s_out,
                    double t1, double t2, double t3, double t4,
                    double t5, double t6, double t7,
                    double A,  double J)
        {
            const double dt = 1.0 / fs;
            const double T  = t1+t2+t3+t4+t5+t6+t7;

            const size_t N = static_cast<size_t>(std::floor(T/dt)) + 1;
            a_out.assign(N, 0.0);
            v_out.assign(N, 0.0);
            s_out.assign(N, 0.0);

            // Segment boundaries (cumulative)
            const double b1 = t1;
            const double b2 = b1 + t2;
            const double b3 = b2 + t3;                 // = t1+t2+t1
            const double b4 = b3 + t4;
            const double b5 = b4 + t5;                 // + t1
            const double b6 = b5 + t6;                 // + t2
            const double b7 = b6 + t7;                 // + t1 == T

            // Useful constants at segment boundaries
            const double v1 = 0.5*J*t1*t1;
            const double s1 = (1.0/6.0)*J*t1*t1*t1;

            const double v2 = v1 + A*t2;
            const double s2 = s1 + v1*t2 + 0.5*A*t2*t2;

            // End of seg3 (accel phase finished, accel back to 0)
            const double v3 = v2 + (A*t1 - 0.5*J*t1*t1);        // = v2 + v1
            const double s3 = s2 + v2*t1 + 0.5*A*t1*t1 - (1.0/6.0)*J*t1*t1*t1;

            const double v_acc = v3;                             // cruise speed
            const double s4 = s3 + v_acc*t4;

            // Start of seg6 (after seg5)
            const double v5e = v_acc - 0.5*J*t1*t1;              // end of seg5
            const double s5e = s4 + v_acc*t1 - (1.0/6.0)*J*t1*t1*t1;

            // Start of seg7 (after seg6)
            const double v6e = v5e - A*t2;                       // = v1
            const double s6e = s5e + v5e*t2 - 0.5*A*t2*t2;

            // Sample
            for (size_t i = 0; i < N; ++i) {
                const double t = i*dt;
                double a, v, s;

                if (t < b1) {                     // seg1: +J
                    const double tau = t;
                    a = J*tau;
                    v = 0.5*J*tau*tau;
                    s = (1.0/6.0)*J*tau*tau*tau;
                }
                else if (t < b2) {                // seg2: +A
                    const double tau = t - b1;
                    a = A;
                    v = v1 + A*tau;
                    s = s1 + v1*tau + 0.5*A*tau*tau;
                }
                else if (t < b3) {                // seg3: -J
                    const double tau = t - b2;
                    a = A - J*tau;
                    v = v2 + A*tau - 0.5*J*tau*tau;
                    s = s2 + v2*tau + 0.5*A*tau*tau - (1.0/6.0)*J*tau*tau*tau;
                }
                else if (t < b4) {                // seg4: cruise
                    const double tau = t - b3;
                    a = 0.0;
                    v = v_acc;
                    s = s3 + v_acc*tau;
                }
                else if (t < b5) {                // seg5: -J
                    const double tau = t - b4;
                    a = -J*tau;
                    v = v_acc - 0.5*J*tau*tau;
                    s = s4 + v_acc*tau - (1.0/6.0)*J*tau*tau*tau;
                }
                else if (t < b6) {                // seg6: -A
                    const double tau = t - b5;
                    a = -A;
                    v = v5e - A*tau;
                    s = s5e + v5e*tau - 0.5*A*tau*tau;
                }
                else {                             // seg7: +J
                    const double tau = t - b6;
                    a = -A + J*tau;
                    v = v6e - A*tau + 0.5*J*tau*tau;
                    s = s6e + v6e*tau - 0.5*A*tau*tau + (1.0/6.0)*J*tau*tau*tau;
                }

                a_out[i] = a;
                v_out[i] = v;
                s_out[i] = s;
            }
            std::cout << "Computed acc, vel and pos\n";
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

        void print_profile() const {
            std::cout << "Case " << case_id << " selected\n";
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "t1="<<t1<<" t2="<<t2<<" t3="<<t3<<" t4="<<t4
                    <<" t5="<<t5<<" t6="<<t6<<" t7="<<t7<<"\n";
        }

        void update(double newAccel, double newFeed, double newJerk){

            accel_max = newAccel;
            max_feedrate = newFeed;
            jerk_max = newJerk;
            

        }

    private:
        static inline double posRootQuad(double a,double b,double c, double eps=1e-9){
            double D = b*b - 4*a*c;
            if (D < -eps) return -1.0;                 // no real roots
            double r = (std::abs(D) <= eps) ? 0.0 : std::sqrt(D);
            double x1 = (-b + r)/(2*a);
            double x2 = (-b - r)/(2*a);

            // Accept non-negative within tolerance; clamp tiny negatives to 0
            double best = 1e300;
            if (x1 >= -eps) best = std::min(best, std::max(0.0, x1));
            if (x2 >= -eps) best = std::min(best, std::max(0.0, x2));

            return (best==1e300) ? -1.0 : best;
        }

        static inline double Sad(double A,double t1,double t2){
            return 2.0 * (A*t1*t1 + 1.5*A*t1*t2 + 0.5*A*t2*t2);
        }

};