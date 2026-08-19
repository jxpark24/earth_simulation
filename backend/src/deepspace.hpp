// SDP4 deep-space extensions to SGP4.
//
// Objects with a period of 225 minutes or more (GEO, GPS/Galileo, Molniya)
// are perturbed enough by the Sun and Moon -- and by resonance with the
// Earth's tesseral harmonics -- that the near-Earth model alone drifts badly.
// These routines add those effects.
#pragma once

#include "sgp4.hpp"

namespace sgp4 {
namespace ds {

// Intermediate quantities shared between dscom(), dsinit() and dpper().
struct Common {
    double snodm, cnodm, sinim, cosim, sinomm, cosomm, day, rtemsq, gam;
    double em, emsq, nm;
    double s1, s2, s3, s4, s5, s6, s7;
    double ss1, ss2, ss3, ss4, ss5, ss6, ss7;
    double sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33;
    double z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33;
};

// Solar and lunar terms evaluated once at epoch.
Common dscom(double epoch, double ep, double argpp, double tc,
             double inclp, double nodep, double np, Satellite& s);

// Resonance setup and lunar-solar secular rates.
void dsinit(Satellite& s, const Common& c, double tc, double xpidot,
            double& em, double& argpm, double& inclm, double& mm,
            double& nm, double& nodem);

// Lunar-solar long-period periodics, applied every propagation step.
void dpper(const Satellite& s, double t,
           double& ep, double& inclp, double& nodep, double& argpp, double& mp);

// Deep-space secular effects, including numerical integration of resonance.
void dspace(const Satellite& s, double t, double tc,
            double& em, double& argpm, double& inclm, double& mm,
            double& nodem, double& nm);

}  // namespace ds
}  // namespace sgp4
