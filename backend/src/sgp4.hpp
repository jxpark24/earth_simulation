// SGP4 orbital propagator (WGS-72), near-Earth model.
// Follows the Vallado/Kelso reference formulation of the SGP4 theory used to
// propagate NORAD two-line element sets.
#pragma once

#include <string>

namespace sgp4 {

constexpr double kPi     = 3.14159265358979323846;
constexpr double kTwoPi  = 2.0 * kPi;
constexpr double kDeg2Rad = kPi / 180.0;

// WGS-72 — the gravity model the TLEs themselves are generated against.
constexpr double kRadiusEarthKm = 6378.135;
constexpr double kMu            = 398600.8;
constexpr double kJ2            = 0.001082616;
constexpr double kJ3            = -0.00000253881;
constexpr double kJ4            = -0.00000165597;
constexpr double kJ3OverJ2      = kJ3 / kJ2;

// Orbital elements parsed straight out of a TLE, plus every derived constant
// SGP4 precomputes once at initialisation so propagation stays cheap.
struct Satellite {
    std::string name;
    int    norad_id   = 0;
    bool   valid      = false;
    bool   deep_space = false;   // period >= 225 min: reduced-accuracy model
    int    error      = 0;

    // Epoch.
    double jdsatepoch = 0.0;

    // Mean elements at epoch.
    double bstar = 0.0, ecco = 0.0, inclo = 0.0, nodeo = 0.0, argpo = 0.0;
    double mo = 0.0, no_kozai = 0.0, no_unkozai = 0.0;

    // Derived / precomputed.
    int    isimp = 0;
    double aycof = 0, con41 = 0, cc1 = 0, cc4 = 0, cc5 = 0, d2 = 0, d3 = 0, d4 = 0;
    double delmo = 0, eta = 0, argpdot = 0, omgcof = 0, sinmao = 0, t2cof = 0;
    double t3cof = 0, t4cof = 0, t5cof = 0, x1mth2 = 0, x7thm1 = 0, mdot = 0;
    double nodedot = 0, xlcof = 0, xmcof = 0, nodecf = 0;

    // --- deep-space (SDP4) terms, populated only when deep_space is true ---
    double gsto = 0;                 // sidereal time at epoch
    int    irez = 0;                 // 0 none, 1 synchronous, 2 half-day resonance
    // Lunar-solar long-period coefficients.
    double e3=0, ee2=0, peo=0, pgho=0, pho=0, pinco=0, plo=0;
    double se2=0, se3=0, sgh2=0, sgh3=0, sgh4=0, sh2=0, sh3=0;
    double si2=0, si3=0, sl2=0, sl3=0, sl4=0;
    double xgh2=0, xgh3=0, xgh4=0, xh2=0, xh3=0, xi2=0, xi3=0;
    double xl2=0, xl3=0, xl4=0, zmol=0, zmos=0;
    // Secular rates from lunar-solar attraction.
    double dedt=0, didt=0, dmdt=0, dnodt=0, domdt=0;
    // Geopotential resonance coefficients.
    double d2201=0, d2211=0, d3210=0, d3222=0, d4410=0, d4422=0;
    double d5220=0, d5232=0, d5421=0, d5433=0;
    double del1=0, del2=0, del3=0, xfact=0, xlamo=0;
};

struct StateVector {
    double r[3];   // position, km, TEME frame
    double v[3];   // velocity, km/s, TEME frame
    int    error;  // 0 == ok
};

// Parse a TLE line pair (plus optional name) into an initialised Satellite.
Satellite tle_to_satellite(const std::string& name,
                           const std::string& line1,
                           const std::string& line2);

// Propagate to `tsince` minutes relative to the satellite's own epoch.
StateVector propagate(const Satellite& s, double tsince);

// --- Time utilities ---------------------------------------------------------

double julian_date(int year, int month, int day, int hour, int minute, double second);

// Greenwich Mean Sidereal Time (radians) for a given Julian date.
double gmst(double jd);

}  // namespace sgp4
