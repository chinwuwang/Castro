#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();
    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        if (coord_type == 1 || coord_type == 2) {
            problem::center[n] = 0.0_rt;
        }
        else {
            problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
        }
    }

}

#endif
