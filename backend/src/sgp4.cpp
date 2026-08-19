#include "sgp4.hpp"
#include "deepspace.hpp"

#include <cmath>
#include <cstdlib>

namespace sgp4 {
namespace {

// xke = 60 / sqrt(Re^3 / mu) — converts between the canonical time unit SGP4
// works in and minutes.
const double kXke = 60.0 / std::sqrt(kRadiusEarthKm * kRadiusEarthKm * kRadiusEarthKm / kMu);
constexpr double kX2o3 = 2.0 / 3.0;

double fmod2p(double x) {
    double r = std::fmod(x, kTwoPi);
    if (r < 0.0) r += kTwoPi;
    return r;
}

// Pull a fixed-width field out of a TLE line without tripping over short lines.
double field(const std::string& s, size_t start, size_t len) {
    if (start >= s.size()) return 0.0;
    return std::atof(s.substr(start, std::min(len, s.size() - start)).c_str());
}

// TLE exponent fields look like " 14375-3" meaning 0.14375e-3.
double decimal_exponent_field(const std::string& s, size_t start, size_t len) {
    if (start >= s.size()) return 0.0;
    std::string f = s.substr(start, std::min(len, s.size() - start));
    std::string mantissa, exponent;
    size_t i = 0;
    while (i < f.size() && f[i] == ' ') ++i;
    if (i < f.size() && (f[i] == '-' || f[i] == '+')) mantissa += f[i++];
    while (i < f.size() && (isdigit(f[i]) || f[i] == '.')) mantissa += f[i++];
    if (i < f.size() && (f[i] == '-' || f[i] == '+')) exponent += f[i++];
    while (i < f.size() && isdigit(f[i])) exponent += f[i++];
    if (mantissa.empty() || mantissa == "-" || mantissa == "+") return 0.0;
    double m = std::atof(("0." + (mantissa[0] == '-' || mantissa[0] == '+'
                                      ? mantissa.substr(1) : mantissa)).c_str());
    if (mantissa[0] == '-') m = -m;
    int e = exponent.empty() ? 0 : std::atoi(exponent.c_str());
    return m * std::pow(10.0, e);
}

// One-time setup: un-Kozai the mean motion and precompute the secular rates and
// drag coefficients that propagation reuses on every step.
void sgp4init(Satellite& s) {
    const double ss      = 78.0 / kRadiusEarthKm + 1.0;
    const double qzms2t  = std::pow((120.0 - 78.0) / kRadiusEarthKm, 4);

    const double eccsq  = s.ecco * s.ecco;
    const double omeosq = 1.0 - eccsq;
    const double rteosq = std::sqrt(omeosq);
    const double cosio  = std::cos(s.inclo);
    const double cosio2 = cosio * cosio;

    // Recover the unperturbed ("un-Kozai'd") mean motion from the TLE value.
    const double ak  = std::pow(kXke / s.no_kozai, kX2o3);
    const double d1  = 0.75 * kJ2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    double del       = d1 / (ak * ak);
    const double adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
    del              = d1 / (adel * adel);
    s.no_unkozai     = s.no_kozai / (1.0 + del);

    const double ao    = std::pow(kXke / s.no_unkozai, kX2o3);
    const double sinio = std::sin(s.inclo);
    const double po    = ao * omeosq;
    const double con42 = 1.0 - 5.0 * cosio2;
    s.con41            = -con42 - cosio2 - cosio2;
    const double posq  = po * po;
    const double rp    = ao * (1.0 - s.ecco);

    s.deep_space = (kTwoPi / s.no_unkozai) >= 225.0;

    if (omeosq < 0.0 || s.no_unkozai < 0.0) { s.error = 4; return; }

    // A *mean* perigee below the reference sphere is not grounds for rejection:
    // highly eccentric orbits (the Cluster probes, e ~ 0.91) sit slightly under
    // it yet propagate fine, spending almost all their time far outside. Decay
    // is caught per-position by the mrt < 1 check during propagation instead.

    s.isimp = 0;
    // Very low perigee, or deep space: skip the higher-order drag expansion.
    if (rp < (220.0 / kRadiusEarthKm + 1.0) || s.deep_space) s.isimp = 1;

    double sfour  = ss;
    double qzms24 = qzms2t;
    const double perige = (rp - 1.0) * kRadiusEarthKm;

    // The atmospheric drag model uses a different reference altitude down low.
    if (perige < 156.0) {
        sfour = perige - 78.0;
        if (perige < 98.0) sfour = 20.0;
        qzms24 = std::pow((120.0 - sfour) / kRadiusEarthKm, 4);
        sfour  = sfour / kRadiusEarthKm + 1.0;
    }
    const double pinvsq = 1.0 / posq;

    const double tsi   = 1.0 / (ao - sfour);
    s.eta              = ao * s.ecco * tsi;
    const double etasq = s.eta * s.eta;
    const double eeta  = s.ecco * s.eta;
    const double psisq = std::fabs(1.0 - etasq);
    const double coef  = qzms24 * std::pow(tsi, 4);
    const double coef1 = coef / std::pow(psisq, 3.5);

    const double cc2 = coef1 * s.no_unkozai *
        (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
         0.375 * kJ2 * tsi / psisq * s.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    s.cc1 = s.bstar * cc2;

    double cc3 = 0.0;
    if (s.ecco > 1.0e-4)
        cc3 = -2.0 * coef * tsi * kJ3OverJ2 * s.no_unkozai * sinio / s.ecco;

    s.x1mth2 = 1.0 - cosio2;
    s.cc4 = 2.0 * s.no_unkozai * coef1 * ao * omeosq *
        (s.eta * (2.0 + 0.5 * etasq) + s.ecco * (0.5 + 2.0 * etasq) -
         kJ2 * tsi / (ao * psisq) *
             (-3.0 * s.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
              0.75 * s.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) *
                  std::cos(2.0 * s.argpo)));
    s.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

    const double cosio4 = cosio2 * cosio2;
    const double temp1  = 1.5 * kJ2 * pinvsq * s.no_unkozai;
    const double temp2  = 0.5 * temp1 * kJ2 * pinvsq;
    const double temp3  = -0.46875 * kJ4 * pinvsq * pinvsq * s.no_unkozai;

    // Secular rates of mean anomaly, argument of perigee and RAAN under J2/J4.
    s.mdot    = s.no_unkozai + 0.5 * temp1 * rteosq * s.con41 +
                0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
    s.argpdot = -0.5 * temp1 * con42 +
                0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
    const double xhdot1 = -temp1 * cosio;
    s.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                          2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;

    s.omgcof = s.bstar * cc3 * std::cos(s.argpo);
    s.xmcof  = 0.0;
    if (s.ecco > 1.0e-4) s.xmcof = -kX2o3 * coef * s.bstar / eeta;
    s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1;
    s.t2cof  = 1.5 * s.cc1;

    // Guard the near-polar singularity at cos(i) = -1.
    if (std::fabs(cosio + 1.0) > 1.5e-12)
        s.xlcof = -0.25 * kJ3OverJ2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    else
        s.xlcof = -0.25 * kJ3OverJ2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;

    s.aycof  = -0.5 * kJ3OverJ2 * sinio;
    const double delmotemp = 1.0 + s.eta * std::cos(s.mo);
    s.delmo  = delmotemp * delmotemp * delmotemp;
    s.sinmao = std::sin(s.mo);
    s.x7thm1 = 7.0 * cosio2 - 1.0;

    if (s.isimp != 1) {
        const double cc1sq = s.cc1 * s.cc1;
        s.d2 = 4.0 * ao * tsi * cc1sq;
        const double temp = s.d2 * tsi * s.cc1 / 3.0;
        s.d3 = (17.0 * ao + sfour) * temp;
        s.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * s.cc1;
        s.t3cof = s.d2 + 2.0 * cc1sq;
        s.t4cof = 0.25 * (3.0 * s.d3 + s.cc1 * (12.0 * s.d2 + 10.0 * cc1sq));
        s.t5cof = 0.2 * (3.0 * s.d4 + 12.0 * s.cc1 * s.d3 + 6.0 * s.d2 * s.d2 +
                         15.0 * cc1sq * (2.0 * s.d2 + cc1sq));
    }
    // Sidereal time at epoch: the phase reference for resonance terms.
    s.gsto = gmst(s.jdsatepoch);

    if (s.deep_space) {
        const double tc = 0.0;
        ds::Common c = ds::dscom(s.jdsatepoch - 2433281.5, s.ecco, s.argpo, tc,
                                 s.inclo, s.nodeo, s.no_unkozai, s);
        double em = c.em, inclm = s.inclo, nm = c.nm;
        double argpm = 0.0, nodem = 0.0, mm = 0.0;
        ds::dsinit(s, c, tc, s.argpdot + s.nodedot, em, argpm, inclm, mm, nm, nodem);
    }

    s.valid = true;
}

}  // namespace

StateVector propagate(const Satellite& s, double tsince) {
    StateVector out{};
    out.error = 0;
    if (!s.valid) { out.error = s.error ? s.error : 1; return out; }

    // --- Secular update of the mean elements ---
    const double xmdf   = s.mo + s.mdot * tsince;
    const double argpdf = s.argpo + s.argpdot * tsince;
    const double nodedf = s.nodeo + s.nodedot * tsince;
    double argpm = argpdf;
    double mm    = xmdf;
    const double t2 = tsince * tsince;
    double nodem = nodedf + s.nodecf * t2;

    double tempa = 1.0 - s.cc1 * tsince;
    double tempe = s.bstar * s.cc4 * tsince;
    double templ = s.t2cof * t2;

    // Full drag expansion (skipped for low-perigee and deep-space objects).
    if (s.isimp != 1) {
        const double delomg   = s.omgcof * tsince;
        const double delmtemp = 1.0 + s.eta * std::cos(xmdf);
        const double delm     = s.xmcof * (delmtemp * delmtemp * delmtemp - s.delmo);
        const double temp     = delomg + delm;
        mm    = xmdf + temp;
        argpm = argpdf - temp;
        const double t3 = t2 * tsince;
        const double t4 = t3 * tsince;
        tempa = tempa - s.d2 * t2 - s.d3 * t3 - s.d4 * t4;
        tempe = tempe + s.bstar * s.cc5 * (std::sin(mm) - s.sinmao);
        templ = templ + s.t3cof * t3 + t4 * (s.t4cof + tsince * s.t5cof);
    }

    double nm    = s.no_unkozai;
    double em    = s.ecco;
    double inclm = s.inclo;

    // Lunar-solar secular drift plus resonance integration.
    if (s.deep_space)
        ds::dspace(s, tsince, tsince, em, argpm, inclm, mm, nodem, nm);

    if (nm <= 0.0) { out.error = 2; return out; }
    const double am = std::pow(kXke / nm, kX2o3) * tempa * tempa;
    nm = kXke / std::pow(am, 1.5);
    em = em - tempe;

    if (em >= 1.0 || em < -0.001) { out.error = 1; return out; }
    if (em < 1.0e-6) em = 1.0e-6;

    mm += s.no_unkozai * templ;
    double xlm = mm + argpm + nodem;

    nodem = fmod2p(nodem);
    argpm = fmod2p(argpm);
    xlm   = fmod2p(xlm);
    mm    = fmod2p(xlm - argpm - nodem);

    // These stay local: propagate() must not mutate the shared Satellite, or
    // the OpenMP batch loop would race.
    double ep = em, xincp = inclm, argpp = argpm, nodep = nodem, mp = mm;
    double sinip = std::sin(inclm), cosip = std::cos(inclm);
    double aycof = s.aycof, xlcof = s.xlcof;
    double con41 = s.con41, x1mth2 = s.x1mth2, x7thm1 = s.x7thm1;

    if (s.deep_space) {
        ds::dpper(s, tsince, ep, xincp, nodep, argpp, mp);
        if (xincp < 0.0) { xincp = -xincp; nodep += kPi; argpp -= kPi; }
        if (ep < 0.0 || ep > 1.0) { out.error = 3; return out; }

        // Inclination has moved, so the inclination-dependent coefficients
        // have to be rebuilt rather than reused from epoch.
        sinip = std::sin(xincp);
        cosip = std::cos(xincp);
        aycof = -0.5 * kJ3OverJ2 * sinip;
        if (std::fabs(cosip + 1.0) > 1.5e-12)
            xlcof = -0.25 * kJ3OverJ2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
        else
            xlcof = -0.25 * kJ3OverJ2 * sinip * (3.0 + 5.0 * cosip) / 1.5e-12;
        const double cosisq = cosip * cosip;
        con41  = 3.0 * cosisq - 1.0;
        x1mth2 = 1.0 - cosisq;
        x7thm1 = 7.0 * cosisq - 1.0;
    }

    // --- Long-period periodics ---
    const double axnl = ep * std::cos(argpp);
    double temp = 1.0 / (am * (1.0 - ep * ep));
    const double aynl = ep * std::sin(argpp) + temp * aycof;
    const double xl   = mp + argpp + nodep + temp * xlcof * axnl;

    // --- Kepler's equation, solved by Newton-Raphson ---
    const double u = fmod2p(xl - nodep);
    double eo1  = u;
    double tem5 = 9999.9;
    int ktr = 1;
    double sineo1 = 0.0, coseo1 = 0.0;
    while (std::fabs(tem5) >= 1.0e-12 && ktr <= 10) {
        sineo1 = std::sin(eo1);
        coseo1 = std::cos(eo1);
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if (std::fabs(tem5) >= 0.95) tem5 = tem5 > 0.0 ? 0.95 : -0.95;
        eo1 += tem5;
        ++ktr;
    }

    // --- Short-period periodics ---
    const double ecose = axnl * coseo1 + aynl * sineo1;
    const double esine = axnl * sineo1 - aynl * coseo1;
    const double el2   = axnl * axnl + aynl * aynl;
    const double pl    = am * (1.0 - el2);
    if (pl < 0.0) { out.error = 4; return out; }

    const double rl     = am * (1.0 - ecose);
    const double rdotl  = std::sqrt(am) * esine / rl;
    const double rvdotl = std::sqrt(pl) / rl;
    const double betal  = std::sqrt(1.0 - el2);
    temp = esine / (1.0 + betal);
    const double sinu = am / rl * (sineo1 - aynl - axnl * temp);
    const double cosu = am / rl * (coseo1 - axnl + aynl * temp);
    double su = std::atan2(sinu, cosu);
    const double sin2u = (cosu + cosu) * sinu;
    const double cos2u = 1.0 - 2.0 * sinu * sinu;

    temp = 1.0 / pl;
    const double temp1 = 0.5 * kJ2 * temp;
    const double temp2 = temp1 * temp;

    const double mrt = rl * (1.0 - 1.5 * temp2 * betal * con41) +
                       0.5 * temp1 * x1mth2 * cos2u;
    su    = su - 0.25 * temp2 * x7thm1 * sin2u;
    const double xnode = nodep + 1.5 * temp2 * cosip * sin2u;
    const double xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
    const double mvt   = rdotl - nm * temp1 * x1mth2 * sin2u / kXke;
    const double rvdot = rvdotl + nm * temp1 * (x1mth2 * cos2u + 1.5 * con41) / kXke;

    if (mrt < 1.0) { out.error = 6; return out; }   // decayed

    // --- Orientation vectors -> TEME position and velocity ---
    const double sinsu = std::sin(su),    cossu = std::cos(su);
    const double snod  = std::sin(xnode), cnod  = std::cos(xnode);
    const double sini  = std::sin(xinc),  cosi  = std::cos(xinc);
    const double xmx = -snod * cosi;
    const double xmy =  cnod * cosi;

    const double ux = xmx * sinsu + cnod * cossu;
    const double uy = xmy * sinsu + snod * cossu;
    const double uz = sini * sinsu;
    const double vx = xmx * cossu - cnod * sinsu;
    const double vy = xmy * cossu - snod * sinsu;
    const double vz = sini * cossu;

    out.r[0] = mrt * ux * kRadiusEarthKm;
    out.r[1] = mrt * uy * kRadiusEarthKm;
    out.r[2] = mrt * uz * kRadiusEarthKm;
    const double vkmps = kRadiusEarthKm * kXke / 60.0;
    out.v[0] = (mvt * ux + rvdot * vx) * vkmps;
    out.v[1] = (mvt * uy + rvdot * vy) * vkmps;
    out.v[2] = (mvt * uz + rvdot * vz) * vkmps;
    return out;
}

double julian_date(int year, int month, int day, int hour, int minute, double second) {
    return 367.0 * year
         - std::floor(7.0 * (year + std::floor((month + 9.0) / 12.0)) * 0.25)
         + std::floor(275.0 * month / 9.0)
         + day + 1721013.5
         + ((second / 60.0 + minute) / 60.0 + hour) / 24.0;
}

double gmst(double jd) {
    const double tut1 = (jd - 2451545.0) / 36525.0;
    double t = -6.2e-6 * tut1 * tut1 * tut1
             + 0.093104 * tut1 * tut1
             + (876600.0 * 3600.0 + 8640184.812866) * tut1
             + 67310.54841;
    t = std::fmod(t * kDeg2Rad / 240.0, kTwoPi);   // seconds -> radians
    if (t < 0.0) t += kTwoPi;
    return t;
}

Satellite tle_to_satellite(const std::string& name,
                           const std::string& line1,
                           const std::string& line2) {
    Satellite s;
    s.name = name;
    if (line1.size() < 63 || line2.size() < 63) { s.error = 2; return s; }

    s.norad_id = static_cast<int>(field(line2, 2, 5));
    s.bstar    = decimal_exponent_field(line1, 53, 8);

    // Epoch: two-digit year + fractional day-of-year.
    int    epochyr  = static_cast<int>(field(line1, 18, 2));
    double epochdays = field(line1, 20, 12);
    int year = epochyr < 57 ? epochyr + 2000 : epochyr + 1900;

    // Day-of-year -> calendar date, then Julian date.
    const bool leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    int mdays[12] = {31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int dayofyr = static_cast<int>(epochdays);
    if (dayofyr < 1) dayofyr = 1;
    int month = 0, remaining = dayofyr;
    while (month < 11 && remaining > mdays[month]) { remaining -= mdays[month]; ++month; }
    const double frac = (epochdays - static_cast<double>(dayofyr)) * 24.0;
    const int    hour = static_cast<int>(frac);
    const double fmin = (frac - hour) * 60.0;
    const int    minute = static_cast<int>(fmin);
    const double second = (fmin - minute) * 60.0;
    s.jdsatepoch = julian_date(year, month + 1, remaining, hour, minute, second);

    s.inclo = field(line2,  8, 8) * kDeg2Rad;
    s.nodeo = field(line2, 17, 8) * kDeg2Rad;
    s.ecco  = field(line2, 26, 7) * 1.0e-7;
    s.argpo = field(line2, 34, 8) * kDeg2Rad;
    s.mo    = field(line2, 43, 8) * kDeg2Rad;
    s.no_kozai = field(line2, 52, 11) * kTwoPi / 1440.0;   // rev/day -> rad/min

    if (s.no_kozai <= 0.0 || s.ecco < 0.0 || s.ecco >= 1.0) { s.error = 3; return s; }

    sgp4init(s);
    return s;
}

}  // namespace sgp4
