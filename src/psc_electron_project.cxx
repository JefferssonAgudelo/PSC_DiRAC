
#include <psc_config.h>

#include <psc.h>
#include <psc.hxx>
#include <output_particles.hxx>

#include <balance.hxx>
#include <particles.hxx>
#include <fields3d.hxx>
#include <push_particles.hxx>
#include <push_fields.hxx>
#include <sort.hxx>
#include <collision.hxx>
#include <bnd_particles.hxx>
#include <bnd.hxx>
#include <bnd_fields.hxx>
#include <marder.hxx>
#include <inject.hxx>
#include <heating.hxx>
#include <setup_particles.hxx>
#include <setup_fields.hxx>

#include "inject_impl.hxx"
#include "../libpsc/psc_heating/psc_heating_impl.hxx"
#include "../libpsc/psc_checks/checks_impl.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#endif

#include "psc_config.hxx"

#include "heating_spot_foil.hxx"



// =================================================================================
// This one is the file that is supossed to be the base file
// In this file the information is given in terms of the inertial electron lenght
// and electron times

// =================================================================================

// =================================================================================

// EDIT to change order / floating point type / cuda / 2d/3d
using dim_t = dim_yz;
//using dim_t = dim_xyz;
using PscConfig = PscConfig1vbecSingle<dim_t>;
//using PscConfig = PscConfig1vbecCuda<dim_t>;

// ======================================================================
// ElectronProject

struct ElectronProject : Psc<PscConfig>
{
  using DIM = PscConfig::dim_t;
  //using Heating_t = typename HeatingSelector<Mparticles_t>::Heating;
  //using Inject_t = typename InjectSelector<Mparticles_t, InjectFoil, DIM>::Inject;

	enum {
  		MY_ION, //0
  		MY_ELECTRON_CORE,
  		MY_ELECTRON_STRAHL,
	};

  // ----------------------------------------------------------------------
  // ctor
  
  ElectronProject()
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");


  // ----------------------------------------------------------------------
    // Parameters for background plasma
    kb = 1.;
    mu0 =1.;
    mi_over_me_ = 100.;
    vA_over_c_ = .1;
    
    B0_ = vA_over_c_;
    double mi = 1.;
    double me = 1. / mi_over_me_;

    ni_ = 1.;
    ne_core_ = .97;
    ne_strahl_ = .03;    

    vA_i_over_c_ = vA_over_c_;
    vA_e_core_over_c_ = vA_over_c_ * ((mi * ni_) / (me * ne_core_));    
    vA_e_strahl_over_c_ = vA_over_c_ * ((mi * ni_) / (me * ne_strahl_)); 
    
    
    // These values are in units of c from the fitting    
    Vi_th_per_ = .00016;
    Ve_core_th_per_ = .00826;
    Ve_strahl_th_per_ = .0112;
    Vi_th_par_ = .00016;
    Ve_core_th_par_ = .00821;
    Ve_strahl_th_par_ = .0172;
 
    //beta_i_ = 1.;
    //beta_e_core_ = .97;
    //beta_e_strahl_ = .03; 
    
    beta_i_per_ = .768;
    beta_i_par_ = .768;
    beta_e_core_per_ = 1.07;
    beta_e_core_par_ = 1.06;
    beta_e_strahl_per_ = .0518;
    beta_e_strahl_par_ = .122;
        

	Ti_per_ = beta_i_per_ * sqr(B0_) / ( 2. * kb * mu0 *ni_ );
	Te_core_per_ = beta_e_core_per_ * sqr(B0_) / ( 2. * kb * mu0 *ne_core_ ); 
	Te_strahl_per_ = beta_e_strahl_per_ * sqr(B0_) / (2. * kb * mu0 *ne_strahl_);     
	
	Ti_par_ = beta_i_par_ * sqr(B0_) / ( 2. * kb * mu0 *ni_ );
	Te_core_par_ = beta_e_core_par_ * sqr(B0_) / ( 2. * kb * mu0 *ne_core_ ); 
	Te_strahl_par_ = beta_e_strahl_par_ * sqr(B0_) / (2. * kb * mu0 *ne_strahl_);  
	
	// Bulk velocities. These values are in units of c from the fitting 
	Vi_per_ = 0.;
	Vi_par_ = 0.;	
	Ve_core_per_ = 0.;
	Ve_core_par_ = 0.;
	Ve_strahl_per_ = 0.;
	Ve_strahl_par_ = -0.0124;	   
		
  // ----------------------------------------------------------------------	
	
            
    //double d_i = sqrt(kinds[MY_ION].m / kinds[MY_ION].q);
    double output_particle_interval;
    double d_i =1.;
    
    mpi_printf(comm, "d_e = %g, d_i = %g\n", 1./sqrt(mi_over_me_) , d_i); //assuming ne = ni
    mpi_printf(comm, "lambda_Di (background) = %g\n", sqrt(beta_i_per_ / 2.) * vA_i_over_c_ * d_i);
    p_.nmax = 101;
    p_.cfl = 0.75;

    // -- setup particle kinds
    // last population ("e") is neutralizing
    Grid_t::Kinds kinds = {{1., mi, "i"}, { -1., me, "e"}}; 
    
    // --- setup domain
#if 0    
    Grid_t::Real3 LL = { 400., 800., 400.*6 }; // domain size (in d_e)
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 40, 80, 4 }; // division into patches
#else
    Grid_t::Real3 LL = { 1., 4., 4. }; // domain size (in d_e). is it really d_e question Jeff?
    Int3 gdims = { 1, 80, 80}; // global number of grid points
    Int3 np = { 1, 4, 4 }; // division into patches
#endif

    if (dim::InvarX::value) { ibn[0] = 0; } // FIXME, wrong place, not for VPIC...
    
    //auto grid_domain = Grid_t::Domain{gdims, LL, -.5 * LL, np}; // What is this LL, -.5? question Jeff
    auto grid_domain = Grid_t::Domain{gdims, LL, {}, np}; // What is this LL, -.5? question Jeff
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 10;

    double dt = p_.cfl * courant_length(grid_domain);
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles_t{grid()});

    // -- Balance
    balance_interval = 50;
    balance_.reset(new Balance_t{balance_interval, .1, false});

    // -- Sort
    // FIXME, needs a way to make sure it gets set?
    sort_interval = 10;

    // -- Collision
    int collision_interval = -1; //10; with no collision is collisionless ? 
    double collision_nu = .1;
    collision_.reset(new Collision_t{grid(), collision_interval, collision_nu});

    // -- Checks
    ChecksParams checks_params{};
    checks_params.continuity_every_step = 50;
    checks_params.continuity_threshold = 1e-4;
    checks_params.continuity_verbose = false;
    checks_.reset(new Checks_t{grid(), comm, checks_params});

    // -- Marder correction
    double marder_diffusion = 0.9;
    int marder_loop = 3;
    bool marder_dump = false;
    marder_interval = 0*5; // What about this. Is this doing something question Jeff
    marder_.reset(new Marder_t(grid(), marder_diffusion, marder_loop, marder_dump));

    
    // -- output fields // In the second simulation do not ask for all the fields information as it is not necessary.
    OutputFieldsCParams outf_params;
    outf_params.output_fields = "e,h,j,n_1st_single,v_1st_single,T_1st_single";
    outf_params.pfield_step = 50;
    outf_.reset(new OutputFieldsC{grid(), outf_params});

    
    // -- output particles 
    // Question: Currently it returns x,y,z,px,py,pz,q,m and w
    // How to ask just for a x,y,z,px,py,pz,q ?
    // How to ask for information of a few particles of each cell ?
    
    //psc_output_particles_descr
	output_particle_interval = 50; // Interval to request particle output   
    if (output_particle_interval > 0) {
	//int lo[3] = { 0, 1, 2 };
	//int hi[3] = { 1, 40, 40 }; 
        psc_output_particles_set_param_int(outp_,"every_step", (int) (output_particle_interval));
		//psc_output_particles_set_param_int3(outp_,"lo", lo);
		//psc_output_particles_set_param_int3(outp_,"hi", hi);
		//psc_output_particles_set_param_int3(outp_,"lo", // something in here is not working
		//(int[3]) { 1, 1, 1}); // when setting (int[3]) { 1, 1, 1}); it runs but the outputs are empty
		//when setting (int[3]) { 0, 0, 0}); there are errors reported in the file output_error.txt	
		//psc_output_particles_set_param_int3(outp_,"hi",
		//int[3]) { 1, 4, 4}); //(int[3]) { 1, 4, 4});
		}
    
    
    // --- partition particles and initial balancing
    mpi_printf(comm, "**** Partitioning...\n");
    auto n_prts_by_patch_old = setup_initial_partition();
    auto n_prts_by_patch_new = balance_->initial(n_prts_by_patch_old);
    // balance::initial does not rebalance particles, because the old way of doing this
    // does't even have the particle data structure created yet -- FIXME?
    mprts_->reset(grid());
    
    mpi_printf(comm, "**** Setting up particles...\n");
    setup_initial_particles(*mprts_, n_prts_by_patch_new);
    
    mpi_printf(comm, "**** Setting up fields...\n");
    setup_initial_fields(*mflds_);

    init();
  }

  void init_npt(int kind, double crd[3], psc_particle_npt& npt)
  {
    switch (kind) {
    case MY_ION:
      npt.n    = ni_;
      npt.T[0] = Ti_per_;
      npt.T[1] = Ti_per_;
      npt.T[2] = Ti_par_;
      npt.p[0] = 0.; 
      npt.p[1] = Vi_per_;
      npt.p[2] = Vi_par_;
      
      break;
    case MY_ELECTRON_CORE:
      npt.n    = ne_core_;
      npt.T[0] = Te_core_per_; // 
      npt.T[1] = Te_core_per_;
      npt.T[2] = Te_core_par_; //parallel
      npt.p[0] = 0.; 
      npt.p[1] = Ve_core_per_;
      npt.p[2] = Ve_core_par_; //parallel
      
      break;
    case MY_ELECTRON_STRAHL:
      npt.n    = ne_strahl_;
      npt.T[0] = Te_strahl_per_;
      npt.T[1] = Te_strahl_per_;
      npt.T[2] = Te_strahl_par_; //parallel
      npt.p[0] = 0.; 
      npt.p[1] = Ve_strahl_per_;
      npt.p[2] = Ve_strahl_par_; //parallel

      break;
    default:
      assert(0);
    }
  }
  
  
  // ----------------------------------------------------------------------
  // setup_initial_partition
  
  std::vector<uint> setup_initial_partition()
  {
    SetupParticles<Mparticles_t> setup_particles;
    return setup_particles.setup_partition(grid(), [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_particles
  
  void setup_initial_particles(Mparticles_t& mprts, std::vector<uint>& n_prts_by_patch)
  {
    SetupParticles<Mparticles_t> setup_particles;
    setup_particles.setup_particles(mprts, n_prts_by_patch, [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }

  // ----------------------------------------------------------------------
  // setup_initial_fields
  
  void setup_initial_fields(MfieldsState& mflds)
  {
    SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
	switch (m) {
	case HZ: return B0_; // Magnetic field in z-direction
	default: return 0.;
	}
      });
  }

  // ----------------------------------------------------------------------


private:
  double B0_;
  double kb;
  double mu0;
  double mi_over_me_;
  double vA_over_c_;  
  
  double vA_i_over_c_ ;
  double vA_e_core_over_c_ ;    
  double vA_e_strahl_over_c_ ; 
  
  double ni_;
  double ne_core_;
  double ne_strahl_;
  
  double beta_i_per_;
  double beta_e_core_per_;
  double beta_e_strahl_per_;
  double beta_i_par_;
  double beta_e_core_par_;
  double beta_e_strahl_par_;
  
  double Vi_th_per_;
  double Ve_core_th_per_;
  double Ve_strahl_th_per_;
  double Vi_th_par_;
  double Ve_core_th_par_;
  double Ve_strahl_th_par_;
    
  double Ti_per_;
  double Te_core_per_;
  double Te_strahl_per_; 
  double Ti_par_;
  double Te_core_par_;
  double Te_strahl_par_;
  
  double Vi_per_;
  double Vi_par_;	
  double Ve_core_per_;
  double Ve_core_par_;
  double Ve_strahl_per_;
  double Ve_strahl_par_;	
   
  double output_particle_interval; 
  //
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);
  
  auto psc = new ElectronProject;

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
