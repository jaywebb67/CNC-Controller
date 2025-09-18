#include <cmath>
#include <algorithm>

static constexpr double M_PI = 3.141592653589793238462643383279502884;

struct Vec3 { double x,y,z; };

static inline Vec3 operator-(const Vec3& a, const Vec3& b){
     return {a.x-b.x,a.y-b.y,a.z-b.z}; 
}

static inline Vec3 operator+(const Vec3& a, const Vec3& b){
     return {a.x+b.x,a.y+b.y,a.z+b.z}; 
}


inline Vec3 operator*(const Vec3& v, double s) { 
    return { v.x*s, v.y*s, v.z*s }; 
}

inline Vec3 operator/(const Vec3& v, double s) { 
    return { v.x/s, v.y/s, v.z/s }; 
}


static inline double dot(const Vec3& a, const Vec3& b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline double norm(const Vec3& a){ 
    return std::sqrt(dot(a,a)); 
}

static inline Vec3   unit(const Vec3& a){ 
    double n=norm(a); 
    return (n>0)? Vec3{a.x/n,a.y/n,a.z/n} : Vec3{0,0,0}; 
}

static inline double clamp01(double x){ 
    return std::max(-1.0, std::min(1.0, x)); 
}
