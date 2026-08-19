#include "deepspace.hpp"

#include <cmath>

namespace sgp4 {
namespace ds {
namespace {

constexpr double kZes = 0.01675;        // solar eccentricity
constexpr double kZel = 0.05490;        // lunar eccentricity
constexpr double kZns = 1.19459e-5;     // solar mean motion, rad/min
constexpr double kZnl = 1.5835218e-4;   // lunar mean motion, rad/min
constexpr double kRptim = 4.37526908801129966e-3;   // Earth rotation, rad/min

double fmod2p(double x) {
    double r = std::fmod(x, kTwoPi);
    if (r < 0.0) r += kTwoPi;
    return r;
}

}  // namespace

Common dscom(double epoch, double ep, double argpp, double tc,
             double inclp, double nodep, double np, Satellite& s) {
    constexpr double c1ss = 2.9864797e-6;
    constexpr double c1l  = 4.7968065e-7;
    constexpr double zsinis = 0.39785416, zcosis = 0.91744867;
    constexpr double zcosgs = 0.1945905,  zsings = -0.98088458;

    Common c{};
    c.nm     = np;
    c.em     = ep;
    c.snodm  = std::sin(nodep);
    c.cnodm  = std::cos(nodep);
    c.sinomm = std::sin(argpp);
    c.cosomm = std::cos(argpp);
    c.sinim  = std::sin(inclp);
    c.cosim  = std::cos(inclp);
    c.emsq   = c.em * c.em;
    const double betasq = 1.0 - c.emsq;
    c.rtemsq = std::sqrt(betasq);

    s.peo = s.pinco = s.plo = s.pgho = s.pho = 0.0;

    // Days since 1900 Jan 0.5, which the lunar ephemeris below is phased to.
    c.day = epoch + 18261.5 + tc / 1440.0;
    const double xnodce = std::fmod(4.5236020 - 9.2422029e-4 * c.day, kTwoPi);
    const double stem = std::sin(xnodce), ctem = std::cos(xnodce);

    // Orientation of the lunar orbit plane.
    const double zcosil = 0.91375164 - 0.03568096 * ctem;
    const double zsinil = std::sqrt(1.0 - zcosil * zcosil);
    const double zsinhl = 0.089683511 * stem / zsinil;
    const double zcoshl = std::sqrt(1.0 - zsinhl * zsinhl);
    c.gam = 5.8351514 + 0.0019443680 * c.day;
    double zx = 0.39785416 * stem / zsinil;
    const double zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx = std::atan2(zx, zy);
    zx = c.gam + zx - xnodce;
    const double zcosgl = std::cos(zx), zsingl = std::sin(zx);

    // Pass 1 accumulates solar terms, pass 2 lunar terms.
    double zcosg = zcosgs, zsing = zsings;
    double zcosi = zcosis, zsini = zsinis;
    double zcosh = c.cnodm, zsinh = c.snodm;
    double cc = c1ss;
    const double xnoi = 1.0 / c.nm;

    for (int lsflg = 1; lsflg <= 2; ++lsflg) {
        const double a1  =  zcosg * zcosh + zsing * zcosi * zsinh;
        const double a3  = -zsing * zcosh + zcosg * zcosi * zsinh;
        const double a7  = -zcosg * zsinh + zsing * zcosi * zcosh;
        const double a8  =  zsing * zsini;
        const double a9  =  zsing * zsinh + zcosg * zcosi * zcosh;
        const double a10 =  zcosg * zsini;
        const double a2  =  c.cosim * a7 + c.sinim * a8;
        const double a4  =  c.cosim * a9 + c.sinim * a10;
        const double a5  = -c.sinim * a7 + c.cosim * a8;
        const double a6  = -c.sinim * a9 + c.cosim * a10;

        const double x1 =  a1 * c.cosomm + a2 * c.sinomm;
        const double x2 =  a3 * c.cosomm + a4 * c.sinomm;
        const double x3 = -a1 * c.sinomm + a2 * c.cosomm;
        const double x4 = -a3 * c.sinomm + a4 * c.cosomm;
        const double x5 =  a5 * c.sinomm;
        const double x6 =  a6 * c.sinomm;
        const double x7 =  a5 * c.cosomm;
        const double x8 =  a6 * c.cosomm;

        const double z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        const double z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        const double z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        double z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * c.emsq;
        double z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * c.emsq;
        double z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * c.emsq;
        const double z11 = -6.0 * a1 * a5 + c.emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        const double z12 = -6.0 * (a1 * a6 + a3 * a5) +
                           c.emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        const double z13 = -6.0 * a3 * a6 + c.emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        const double z21 = 6.0 * a2 * a5 + c.emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        const double z22 = 6.0 * (a4 * a5 + a2 * a6) +
                           c.emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        const double z23 = 6.0 * a4 * a6 + c.emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
        z1 = z1 + z1 + betasq * z31;
        z2 = z2 + z2 + betasq * z32;
        z3 = z3 + z3 + betasq * z33;

        const double s3 = cc * xnoi;
        const double s2 = -0.5 * s3 / c.rtemsq;
        const double s4 = s3 * c.rtemsq;
        const double s1 = -15.0 * c.em * s4;
        const double s5 = x1 * x3 + x2 * x4;
        const double s6 = x2 * x3 + x1 * x4;
        const double s7 = x2 * x4 - x1 * x3;

        if (lsflg == 1) {
            c.ss1 = s1; c.ss2 = s2; c.ss3 = s3; c.ss4 = s4;
            c.ss5 = s5; c.ss6 = s6; c.ss7 = s7;
            c.sz1 = z1; c.sz2 = z2; c.sz3 = z3;
            c.sz11 = z11; c.sz12 = z12; c.sz13 = z13;
            c.sz21 = z21; c.sz22 = z22; c.sz23 = z23;
            c.sz31 = z31; c.sz32 = z32; c.sz33 = z33;
            // Switch the driving body to the Moon for the second pass.
            zcosg = zcosgl; zsing = zsingl;
            zcosi = zcosil; zsini = zsinil;
            zcosh = zcoshl * c.cnodm + zsinhl * c.snodm;
            zsinh = c.snodm * zcoshl - c.cnodm * zsinhl;
            cc = c1l;
        } else {
            c.s1 = s1; c.s2 = s2; c.s3 = s3; c.s4 = s4;
            c.s5 = s5; c.s6 = s6; c.s7 = s7;
            c.z1 = z1; c.z2 = z2; c.z3 = z3;
            c.z11 = z11; c.z12 = z12; c.z13 = z13;
            c.z21 = z21; c.z22 = z22; c.z23 = z23;
            c.z31 = z31; c.z32 = z32; c.z33 = z33;
        }
    }

    s.zmol = fmod2p(4.7199672 + 0.22997150 * c.day - c.gam);
    s.zmos = fmod2p(6.2565837 + 0.017201977 * c.day);

    // Solar coefficients.
    s.se2  =  2.0 * c.ss1 * c.ss6;
    s.se3  =  2.0 * c.ss1 * c.ss7;
    s.si2  =  2.0 * c.ss2 * c.sz12;
    s.si3  =  2.0 * c.ss2 * (c.sz13 - c.sz11);
    s.sl2  = -2.0 * c.ss3 * c.sz2;
    s.sl3  = -2.0 * c.ss3 * (c.sz3 - c.sz1);
    s.sl4  = -2.0 * c.ss3 * (-21.0 - 9.0 * c.emsq) * kZes;
    s.sgh2 =  2.0 * c.ss4 * c.sz32;
    s.sgh3 =  2.0 * c.ss4 * (c.sz33 - c.sz31);
    s.sgh4 = -18.0 * c.ss4 * kZes;
    s.sh2  = -2.0 * c.ss2 * c.sz22;
    s.sh3  = -2.0 * c.ss2 * (c.sz23 - c.sz21);

    // Lunar coefficients.
    s.ee2  =  2.0 * c.s1 * c.s6;
    s.e3   =  2.0 * c.s1 * c.s7;
    s.xi2  =  2.0 * c.s2 * c.z12;
    s.xi3  =  2.0 * c.s2 * (c.z13 - c.z11);
    s.xl2  = -2.0 * c.s3 * c.z2;
    s.xl3  = -2.0 * c.s3 * (c.z3 - c.z1);
    s.xl4  = -2.0 * c.s3 * (-21.0 - 9.0 * c.emsq) * kZel;
    s.xgh2 =  2.0 * c.s4 * c.z32;
    s.xgh3 =  2.0 * c.s4 * (c.z33 - c.z31);
    s.xgh4 = -18.0 * c.s4 * kZel;
    s.xh2  = -2.0 * c.s2 * c.z22;
    s.xh3  = -2.0 * c.s2 * (c.z23 - c.z21);

    return c;
}

void dsinit(Satellite& s, const Common& c, double tc, double xpidot,
            double& em, double& argpm, double& inclm, double& mm,
            double& nm, double& nodem) {
    constexpr double q22 = 1.7891679e-6,  q31 = 2.1460748e-6,  q33 = 2.2123015e-7;
    constexpr double root22 = 1.7891679e-6, root44 = 7.3636953e-9, root54 = 2.1765803e-9;
    constexpr double root32 = 3.7393792e-7, root52 = 1.1428639e-7;
    constexpr double x2o3 = 2.0 / 3.0;

    s.irez = 0;
    // One-day (synchronous) resonance...
    if (nm < 0.0052359877 && nm > 0.0034906585) s.irez = 1;
    // ...and the half-day resonance of eccentric 12-hour orbits.
    if (nm >= 8.26e-3 && nm <= 9.24e-3 && em >= 0.5) s.irez = 2;

    // Solar secular rates.
    const double ses  =  c.ss1 * kZns * c.ss5;
    const double sis  =  c.ss2 * kZns * (c.sz11 + c.sz13);
    const double sls  = -kZns * c.ss3 * (c.sz1 + c.sz3 - 14.0 - 6.0 * c.emsq);
    const double sghs =  c.ss4 * kZns * (c.sz31 + c.sz33 - 6.0);
    double shs        = -kZns * c.ss2 * (c.sz21 + c.sz23);
    // Near-equatorial orbits have an ill-defined node; drop the node term.
    if (inclm < 5.2359877e-2 || inclm > kPi - 5.2359877e-2) shs = 0.0;
    if (c.sinim != 0.0) shs /= c.sinim;
    const double sgs = sghs - c.cosim * shs;

    // Lunar secular rates, added on top.
    s.dedt  = ses + c.s1 * kZnl * c.s5;
    s.didt  = sis + c.s2 * kZnl * (c.z11 + c.z13);
    s.dmdt  = sls - kZnl * c.s3 * (c.z1 + c.z3 - 14.0 - 6.0 * c.emsq);
    const double sghl = c.s4 * kZnl * (c.z31 + c.z33 - 6.0);
    double shll       = -kZnl * c.s2 * (c.z21 + c.z23);
    if (inclm < 5.2359877e-2 || inclm > kPi - 5.2359877e-2) shll = 0.0;
    s.domdt = sgs + sghl;
    s.dnodt = shs;
    if (c.sinim != 0.0) {
        s.domdt -= c.cosim / c.sinim * shll;
        s.dnodt += shll / c.sinim;
    }

    const double theta = fmod2p(s.gsto + tc * kRptim);

    if (s.irez == 0) return;

    const double aonv = std::pow(nm / (60.0 / std::sqrt(kRadiusEarthKm * kRadiusEarthKm *
                                                        kRadiusEarthKm / kMu)), x2o3);

    if (s.irez == 2) {
        // 12-hour geopotential resonance (Molniya-class orbits).
        const double cosisq = c.cosim * c.cosim;
        const double emo = em;
        em = s.ecco;
        const double emsq = s.ecco * s.ecco;
        const double eoc = em * emsq;
        const double g201 = -0.306 - (em - 0.64) * 0.440;
        double g211, g310, g322, g410, g422, g520, g533, g521, g532;

        if (em <= 0.65) {
            g211 =    3.616  -   13.2470 * em +   16.2900 * emsq;
            g310 =  -19.302  +  117.3900 * em -  228.4190 * emsq +  156.591  * eoc;
            g322 =  -18.9068 +  109.7927 * em -  214.6334 * emsq +  146.5816 * eoc;
            g410 =  -41.122  +  242.6940 * em -  471.0940 * emsq +  313.953  * eoc;
            g422 = -146.407  +  841.8800 * em - 1629.014  * emsq + 1083.435  * eoc;
            g520 = -532.114  + 3017.977  * em - 5740.032  * emsq + 3708.276  * eoc;
        } else {
            g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc;
            g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc;
            g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc;
            g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc;
            g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
            if (em > 0.715)
                g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
            else
                g520 =  1464.74 -  4664.75 * em +  3763.64 * emsq;
        }
        if (em < 0.7) {
            g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc;
            g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
            g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4   * eoc;
        } else {
            g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
            g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
            g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
        }

        const double sini2 = c.sinim * c.sinim;
        const double f220 = 0.75 * (1.0 + 2.0 * c.cosim + cosisq);
        const double f221 = 1.5 * sini2;
        const double f321 =  1.875 * c.sinim * (1.0 - 2.0 * c.cosim - 3.0 * cosisq);
        const double f322 = -1.875 * c.sinim * (1.0 + 2.0 * c.cosim - 3.0 * cosisq);
        const double f441 = 35.0 * sini2 * f220;
        const double f442 = 39.3750 * sini2 * sini2;
        const double f522 = 9.84375 * c.sinim *
            (sini2 * (1.0 - 2.0 * c.cosim - 5.0 * cosisq) +
             0.33333333 * (-2.0 + 4.0 * c.cosim + 6.0 * cosisq));
        const double f523 = c.sinim *
            (4.92187512 * sini2 * (-2.0 - 4.0 * c.cosim + 10.0 * cosisq) +
             6.56250012 * (1.0 + 2.0 * c.cosim - 3.0 * cosisq));
        const double f542 = 29.53125 * c.sinim *
            (2.0 - 8.0 * c.cosim + cosisq * (-12.0 + 8.0 * c.cosim + 10.0 * cosisq));
        const double f543 = 29.53125 * c.sinim *
            (-2.0 - 8.0 * c.cosim + cosisq * (12.0 + 8.0 * c.cosim - 10.0 * cosisq));

        const double xno2  = nm * nm;
        const double ainv2 = aonv * aonv;
        double temp1 = 3.0 * xno2 * ainv2;
        double temp  = temp1 * root22;
        s.d2201 = temp * f220 * g201;
        s.d2211 = temp * f221 * g211;
        temp1 *= aonv;
        temp = temp1 * root32;
        s.d3210 = temp * f321 * g310;
        s.d3222 = temp * f322 * g322;
        temp1 *= aonv;
        temp = 2.0 * temp1 * root44;
        s.d4410 = temp * f441 * g410;
        s.d4422 = temp * f442 * g422;
        temp1 *= aonv;
        temp = temp1 * root52;
        s.d5220 = temp * f522 * g520;
        s.d5232 = temp * f523 * g532;
        temp = 2.0 * temp1 * root54;
        s.d5421 = temp * f542 * g521;
        s.d5433 = temp * f543 * g533;

        s.xlamo = fmod2p(s.mo + s.nodeo + s.nodeo - theta - theta);
        s.xfact = s.mdot + s.dmdt + 2.0 * (s.nodedot + s.dnodt - kRptim) - s.no_unkozai;
        em = emo;
    }

    if (s.irez == 1) {
        // Synchronous (geostationary) resonance.
        const double g200 = 1.0 + c.emsq * (-2.5 + 0.8125 * c.emsq);
        const double g310 = 1.0 + 2.0 * c.emsq;
        const double g300 = 1.0 + c.emsq * (-6.0 + 6.60937 * c.emsq);
        const double f220 = 0.75 * (1.0 + c.cosim) * (1.0 + c.cosim);
        const double f311 = 0.9375 * c.sinim * c.sinim * (1.0 + 3.0 * c.cosim)
                          - 0.75 * (1.0 + c.cosim);
        double f330 = 1.0 + c.cosim;
        f330 = 1.875 * f330 * f330 * f330;

        double del1 = 3.0 * nm * nm * aonv * aonv;
        s.del2 = 2.0 * del1 * f220 * g200 * q22;
        s.del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
        s.del1 = del1 * f311 * g310 * q31 * aonv;

        s.xlamo = fmod2p(s.mo + s.nodeo + s.argpo - theta);
        s.xfact = s.mdot + xpidot - kRptim + s.dmdt + s.domdt + s.dnodt - s.no_unkozai;
    }
}

void dpper(const Satellite& s, double t,
           double& ep, double& inclp, double& nodep, double& argpp, double& mp) {
    // Solar contribution.
    double zm = s.zmos + kZns * t;
    double zf = zm + 2.0 * kZes * std::sin(zm);
    double sinzf = std::sin(zf);
    double f2 =  0.5 * sinzf * sinzf - 0.25;
    double f3 = -0.5 * sinzf * std::cos(zf);
    const double ses  = s.se2 * f2 + s.se3 * f3;
    const double sis  = s.si2 * f2 + s.si3 * f3;
    const double sls  = s.sl2 * f2 + s.sl3 * f3 + s.sl4 * sinzf;
    const double sghs = s.sgh2 * f2 + s.sgh3 * f3 + s.sgh4 * sinzf;
    const double shs  = s.sh2 * f2 + s.sh3 * f3;

    // Lunar contribution.
    zm = s.zmol + kZnl * t;
    zf = zm + 2.0 * kZel * std::sin(zm);
    sinzf = std::sin(zf);
    f2 =  0.5 * sinzf * sinzf - 0.25;
    f3 = -0.5 * sinzf * std::cos(zf);
    const double sel  = s.ee2 * f2 + s.e3 * f3;
    const double sil  = s.xi2 * f2 + s.xi3 * f3;
    const double sll  = s.xl2 * f2 + s.xl3 * f3 + s.xl4 * sinzf;
    const double sghl = s.xgh2 * f2 + s.xgh3 * f3 + s.xgh4 * sinzf;
    const double shll = s.xh2 * f2 + s.xh3 * f3;

    double pe   = ses + sel  - s.peo;
    double pinc = sis + sil  - s.pinco;
    double pl   = sls + sll  - s.plo;
    double pgh  = sghs + sghl - s.pgho;
    double ph   = shs + shll - s.pho;

    inclp += pinc;
    ep    += pe;
    const double sinip = std::sin(inclp);
    const double cosip = std::cos(inclp);

    if (inclp >= 0.2) {
        // Well away from the equator the node is well conditioned.
        ph   = ph / sinip;
        pgh  = pgh - cosip * ph;
        argpp += pgh;
        nodep += ph;
        mp    += pl;
    } else {
        // Lyddane modification: near-equatorial orbits are reformulated in
        // terms of the equinoctial elements to dodge the 1/sin(i) singularity.
        const double sinop = std::sin(nodep), cosop = std::cos(nodep);
        double alfdp = sinip * sinop;
        double betdp = sinip * cosop;
        const double dalf =  ph * cosop + pinc * cosip * sinop;
        const double dbet = -ph * sinop + pinc * cosip * cosop;
        alfdp += dalf;
        betdp += dbet;
        nodep = std::fmod(nodep, kTwoPi);
        if (nodep < 0.0) nodep += kTwoPi;
        double xls = mp + argpp + cosip * nodep;
        xls += pl + pgh - pinc * nodep * sinip;
        const double xnoh = nodep;
        nodep = std::atan2(alfdp, betdp);
        if (nodep < 0.0) nodep += kTwoPi;
        if (std::fabs(xnoh - nodep) > kPi)
            nodep += (nodep < xnoh) ? kTwoPi : -kTwoPi;
        mp    += pl;
        argpp  = xls - mp - cosip * nodep;
    }
}

void dspace(const Satellite& s, double t, double tc,
            double& em, double& argpm, double& inclm, double& mm,
            double& nodem, double& nm) {
    constexpr double fasx2 = 0.13130908, fasx4 = 2.8843198, fasx6 = 0.37448087;
    constexpr double g22 = 5.7686396, g32 = 0.95240898, g44 = 1.8014998;
    constexpr double g52 = 1.0508330, g54 = 4.4108898;
    constexpr double stepp = 720.0, stepn = -720.0, step2 = 259200.0;

    const double theta = fmod2p(s.gsto + tc * kRptim);

    // Lunar-solar secular drift, linear in time.
    em    += s.dedt * t;
    inclm += s.didt * t;
    argpm += s.domdt * t;
    nodem += s.dnodt * t;
    mm    += s.dmdt * t;

    if (s.irez == 0) return;

    // Resonance is integrated numerically in fixed 720-minute steps. The
    // integration state is local, so propagation stays pure and thread-safe:
    // every call re-integrates from epoch rather than resuming shared state.
    double atime = 0.0;
    double xni   = s.no_unkozai;
    double xli   = s.xlamo;
    const double delt = (t > 0.0) ? stepp : stepn;

    double xndt = 0.0, xldot = 0.0, xnddt = 0.0, ft = 0.0;
    for (;;) {
        if (s.irez != 2) {
            // Synchronous resonance: three tesseral harmonics.
            xndt  = s.del1 * std::sin(xli - fasx2)
                  + s.del2 * std::sin(2.0 * (xli - fasx4))
                  + s.del3 * std::sin(3.0 * (xli - fasx6));
            xldot = xni + s.xfact;
            xnddt = s.del1 * std::cos(xli - fasx2)
                  + 2.0 * s.del2 * std::cos(2.0 * (xli - fasx4))
                  + 3.0 * s.del3 * std::cos(3.0 * (xli - fasx6));
            xnddt *= xldot;
        } else {
            // Half-day resonance: ten terms, argument-of-perigee dependent.
            const double xomi  = s.argpo + s.argpdot * atime;
            const double x2omi = xomi + xomi;
            const double x2li  = xli + xli;
            xndt = s.d2201 * std::sin(x2omi + xli - g22) + s.d2211 * std::sin(xli - g22)
                 + s.d3210 * std::sin(xomi + xli - g32) + s.d3222 * std::sin(-xomi + xli - g32)
                 + s.d4410 * std::sin(x2omi + x2li - g44) + s.d4422 * std::sin(x2li - g44)
                 + s.d5220 * std::sin(xomi + xli - g52) + s.d5232 * std::sin(-xomi + xli - g52)
                 + s.d5421 * std::sin(xomi + x2li - g54) + s.d5433 * std::sin(-xomi + x2li - g54);
            xldot = xni + s.xfact;
            xnddt = s.d2201 * std::cos(x2omi + xli - g22) + s.d2211 * std::cos(xli - g22)
                  + s.d3210 * std::cos(xomi + xli - g32) + s.d3222 * std::cos(-xomi + xli - g32)
                  + s.d5220 * std::cos(xomi + xli - g52) + s.d5232 * std::cos(-xomi + xli - g52)
                  + 2.0 * (s.d4410 * std::cos(x2omi + x2li - g44)
                         + s.d4422 * std::cos(x2li - g44)
                         + s.d5421 * std::cos(xomi + x2li - g54)
                         + s.d5433 * std::cos(-xomi + x2li - g54));
            xnddt *= xldot;
        }

        if (std::fabs(t - atime) < stepp) { ft = t - atime; break; }
        xli   += xldot * delt + xndt * step2;
        xni   += xndt * delt + xnddt * step2;
        atime += delt;
    }

    nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
    const double xl = xli + xldot * ft + xndt * ft * ft * 0.5;
    if (s.irez != 1) mm = xl - 2.0 * nodem + 2.0 * theta;
    else             mm = xl - nodem - argpm + theta;
}

}  // namespace ds
}  // namespace sgp4
