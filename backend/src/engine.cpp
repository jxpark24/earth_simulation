// pybind11 bindings: a Propagator that holds the whole satellite catalogue and
// advances every object to a requested epoch in one batched, parallel call.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <string>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sgp4.hpp"

namespace py = pybind11;

class Propagator {
public:
    // Load a catalogue from parallel arrays of TLE names/lines.
    void load(const std::vector<std::string>& names,
              const std::vector<std::string>& l1,
              const std::vector<std::string>& l2) {
        sats_.clear();
        sats_.reserve(names.size());
        for (size_t i = 0; i < names.size(); ++i)
            sats_.push_back(sgp4::tle_to_satellite(names[i], l1[i], l2[i]));
    }

    size_t size()  const { return sats_.size(); }
    size_t valid() const {
        size_t n = 0;
        for (const auto& s : sats_) if (s.valid) ++n;
        return n;
    }

    // Propagate every satellite to Julian date `jd`.
    // Returns an (N,3) float32 array of TEME positions in km. Satellites that
    // error out (decayed, numerically diverged) are emitted as NaN so the
    // frontend can skip them without a separate mask.
    py::array_t<float> propagate_all(double jd) {
        const size_t n = sats_.size();
        py::array_t<float> result({static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(3)});
        auto buf = result.mutable_unchecked<2>();

        {
            py::gil_scoped_release release;   // let the C++ loop run truly parallel
            const Satellites& sats = sats_;
            #pragma omp parallel for schedule(static)
            for (long long i = 0; i < static_cast<long long>(n); ++i) {
                const sgp4::Satellite& s = sats[i];
                const double tsince = (jd - s.jdsatepoch) * 1440.0;   // days -> minutes
                sgp4::StateVector sv = sgp4::propagate(s, tsince);
                if (sv.error != 0 || !std::isfinite(sv.r[0])) {
                    buf(i, 0) = buf(i, 1) = buf(i, 2) = std::nanf("");
                } else {
                    buf(i, 0) = static_cast<float>(sv.r[0]);
                    buf(i, 1) = static_cast<float>(sv.r[1]);
                    buf(i, 2) = static_cast<float>(sv.r[2]);
                }
            }
        }
        return result;
    }

    // Position + velocity for a single satellite — used for the info panel.
    py::dict state_of(size_t index, double jd) {
        py::dict d;
        if (index >= sats_.size()) { d["error"] = 99; return d; }
        const sgp4::Satellite& s = sats_[index];
        sgp4::StateVector sv = sgp4::propagate(s, (jd - s.jdsatepoch) * 1440.0);
        const double R = std::sqrt(sv.r[0]*sv.r[0] + sv.r[1]*sv.r[1] + sv.r[2]*sv.r[2]);
        const double V = std::sqrt(sv.v[0]*sv.v[0] + sv.v[1]*sv.v[1] + sv.v[2]*sv.v[2]);
        d["error"]        = sv.error;
        d["altitude_km"]  = R - sgp4::kRadiusEarthKm;
        d["speed_kms"]    = V;
        d["period_min"]   = s.no_unkozai > 0 ? sgp4::kTwoPi / s.no_unkozai : 0.0;
        d["inclination"]  = s.inclo * 180.0 / sgp4::kPi;
        d["eccentricity"] = s.ecco;
        return d;
    }

    // Sample one full orbital period into an (steps,3) array, for drawing the
    // selected satellite's ground-independent orbit path.
    py::array_t<float> orbit_track(size_t index, double jd, int steps) {
        if (steps < 2) steps = 2;
        py::array_t<float> result({static_cast<py::ssize_t>(steps), static_cast<py::ssize_t>(3)});
        auto buf = result.mutable_unchecked<2>();
        if (index >= sats_.size()) return result;
        const sgp4::Satellite& s = sats_[index];
        const double period = s.no_unkozai > 0 ? sgp4::kTwoPi / s.no_unkozai : 90.0;
        const double t0 = (jd - s.jdsatepoch) * 1440.0;
        for (int i = 0; i < steps; ++i) {
            sgp4::StateVector sv = sgp4::propagate(s, t0 + period * i / (steps - 1));
            const bool bad = sv.error != 0 || !std::isfinite(sv.r[0]);
            for (int k = 0; k < 3; ++k)
                buf(i, k) = bad ? std::nanf("") : static_cast<float>(sv.r[k]);
        }
        return result;
    }

    // Per-satellite metadata for the catalogue the frontend receives once.
    py::list catalog() {
        py::list out;
        for (size_t i = 0; i < sats_.size(); ++i) {
            const auto& s = sats_[i];
            py::dict d;
            d["i"]      = i;
            d["name"]   = s.name;
            d["id"]     = s.norad_id;
            d["valid"]  = s.valid;
            d["deep"]   = s.deep_space;
            d["incl"]   = s.inclo * 180.0 / sgp4::kPi;
            d["period"] = s.no_unkozai > 0 ? sgp4::kTwoPi / s.no_unkozai : 0.0;
            out.append(d);
        }
        return out;
    }

private:
    using Satellites = std::vector<sgp4::Satellite>;
    Satellites sats_;
};

PYBIND11_MODULE(sat_engine, m) {
    m.doc() = "SGP4 satellite propagation engine (C++)";

    py::class_<Propagator>(m, "Propagator")
        .def(py::init<>())
        .def("load",          &Propagator::load)
        .def("size",          &Propagator::size)
        .def("valid",         &Propagator::valid)
        .def("propagate_all", &Propagator::propagate_all, py::arg("jd"))
        .def("state_of",      &Propagator::state_of, py::arg("index"), py::arg("jd"))
        .def("orbit_track",   &Propagator::orbit_track,
             py::arg("index"), py::arg("jd"), py::arg("steps") = 180)
        .def("catalog",       &Propagator::catalog);

    m.def("julian_date", &sgp4::julian_date,
          py::arg("year"), py::arg("month"), py::arg("day"),
          py::arg("hour"), py::arg("minute"), py::arg("second"));
    m.def("gmst", &sgp4::gmst, py::arg("jd"),
          "Greenwich Mean Sidereal Time in radians.");

    m.attr("EARTH_RADIUS_KM") = sgp4::kRadiusEarthKm;
#ifdef _OPENMP
    m.attr("openmp") = true;
#else
    m.attr("openmp") = false;
#endif
}
