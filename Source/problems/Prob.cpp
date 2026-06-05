/* Implementations of functions in Problem.H go here */

#include <Castro.H>

using namespace amrex;

/* Implementations of functions in Problem.H go here */

#include <iomanip>

#ifdef DO_PROBLEM_POST_TIMESTEP
void
Castro::problem_post_timestep ()
{
    BL_PROFILE("Castro::problem_post_timestep()");

    if (level != 0) return;   // run once per coarse step, after the sync

#ifdef MHD
    Vector<Real> divB_lev(parent->finestLevel()+1, 0.0_rt);

    for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
        Castro& clev = getLevel(lev);

        const int ng = 1;
        const DistributionMapping& dm = clev.get_new_data(State_Type).DistributionMap();
        const BoxArray& ba           = clev.get_new_data(State_Type).boxArray();

        MultiFab Bx(clev.get_new_data(Mag_Type_x).boxArray(), dm, 1, ng);
        MultiFab By(clev.get_new_data(Mag_Type_y).boxArray(), dm, 1, ng);
        MultiFab Bz(clev.get_new_data(Mag_Type_z).boxArray(), dm, 1, ng);

        Real t = clev.state[State_Type].curTime();
        clev.FillPatchMHD(t, Bx, By, Bz, ng);

        iMultiFab mask(ba, dm, 1, 0);
        if (lev < parent->finestLevel()) {
            mask = amrex::makeFineMask(ba, dm,
                       getLevel(lev+1).get_new_data(State_Type).boxArray(),
                       parent->refRatio(lev), 1, 0);
        } else {
            mask.setVal(1);
        }

        const auto dx = clev.geom.CellSizeArray();

        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(clev.get_new_data(State_Type), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& box = mfi.tilebox();
            auto bxa = Bx.array(mfi);
            auto bya = By.array(mfi);
            auto bza = Bz.array(mfi);
            auto m   = mask.const_array(mfi);
            reduce_op.eval(box, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                if (m(i,j,k) == 0) { return { 0.0_rt }; }   // skip covered cells
                Real d = (bxa(i+1,j,k) - bxa(i,j,k))/dx[0]
                       + (bya(i,j+1,k) - bya(i,j,k))/dx[1]
                       + (bza(i,j,k+1) - bza(i,j,k))/dx[2];
                return { std::abs(d) };
            });
        }   // <-- close MFIter loop

        divB_lev[lev] = amrex::get<0>(reduce_data.value());   // after MFIter, inside lev loop
    }   // <-- close lev loop

    for (auto& v : divB_lev) { ParallelDescriptor::ReduceRealMax(v); }

    Real max_divB = 0.0_rt;
    for (auto v : divB_lev) { max_divB = amrex::max(max_divB, v); }

    if (ParallelDescriptor::IOProcessor() && parent->NumDataLogs() > 0) {
        std::ostream& log = parent->DataLog(0);
        if (log.good()) {
            static bool wrote_header = false;
            if (!wrote_header) {
                log << std::setw(10) << "# nstep"
                    << std::setw(24) << "time"
                    << std::setw(24) << "max|divB|";

                for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                std::string lev_str = "divB_lev_" + std::to_string(lev);
                log << std::setw(24) << lev_str;
                }
                log << "\n";
                wrote_header = true;
            }

            log << std::setw(10) << parent->levelSteps(0)
                << std::setw(24) << std::setprecision(14) << std::scientific
                                 << state[State_Type].curTime()
                << std::setw(24) << std::setprecision(14) << std::scientific
                                 << max_divB;

            for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                log << std::setw(24) << std::setprecision(14) << std::scientific
                    << divB_lev[lev];
            }
            log << std::endl;
        }
    }
#endif
}
#endif


#ifdef DO_PROBLEM_POST_SIMULATION
void Castro::problem_post_simulation(Vector<std::unique_ptr<AmrLevel> >& amr_level) {

  // this is a stub post_simulation() routine

  // you can put this in your problem directory to create a custom
  // diagnostic to run at the end of a simulation

  // all of the needed data comes in though the amr_level array,
  // which is a PArray of AmrLevel objects.  The number of levels
  // is simply the size of this array:

  // int nlevels = amr_level.size();

  // To access data, cast a level to a Castro object, e.g. for
  // level 0:

  // Castro* castro = dynamic_cast<Castro*>(&amr_level[0]);

  // then you can get the data, e.g. for state data as:

  // MultiFab& S = castro->get_new_data(State_Type);

  // and if needed, the state descriptor:

  // const StateDescriptor* desc = &castro->desc_lst[State_Type];

  // and then you can get the names of the state data as desc->name(comp)

}
#endif

