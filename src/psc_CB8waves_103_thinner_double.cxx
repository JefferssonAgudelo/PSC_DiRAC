
#include <psc.hxx>
#include <setup_particles.hxx>
#include <setup_fields.hxx>
#include "psc_config.hxx"

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include <psc_config.h>
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

// EDIT to change particle shape order / floating point type / 2d/3d / ...
using dim_t = dim_xyz;
using PscConfig = PscConfig1vbecDouble<dim_t>;

// ======================================================================
// PscTurAW

struct PscTurAW : Psc<PscConfig>
{
  using Dim = PscConfig::dim_t;

  enum {
    TurAW_ELECTRON,
    TurAW_ION,
  };
         
  // ----------------------------------------------------------------------
  // ctor
  
  PscTurAW()
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    // Parameters
    mi_over_me_ = 100.; //Ques_jeff from 10 to 50? Do I have to change something else when I change mi/me?
    vA_over_c_ = .1; //Why 0.1?? 
    amplitude_ = .5;     
    
    beta_e_par_ = 1.; //Ques_jeff what beta 0.1, 1??
    beta_i_par_ = 1.;
    Ti_perp_over_Ti_par_ = 1.;
    Te_perp_over_Te_par_ = 1.;
    
    double B0 = vA_over_c_; //How to define this??
    double Te_par = beta_e_par_ * sqr(B0) / 2.;
    double Te_perp = Te_perp_over_Te_par_ * Te_par;
    double Ti_par = beta_i_par_ * sqr(B0) / 2.;
    double Ti_perp = Ti_perp_over_Ti_par_ * Ti_par;
    double mi = 1.;
    double me = 1. / mi_over_me_;
    //double output_particle_interval;  
    double debye_length_ = 1.* vA_over_c_ *sqrt(beta_i_par_ / 2. );

    mpi_printf(comm, "d_i = 1., d_e = %g\n", sqrt(me));
    mpi_printf(comm, "om_ci = %g, om_ce = %g\n", B0, B0 / me);
    mpi_printf(comm, "\n");
    mpi_printf(comm, "v_i,perp = %g [c] T_i,perp = %g\n", sqrt(2.*Ti_perp), Ti_perp);
    mpi_printf(comm, "v_i,par  = %g [c] T_i,par = %g\n", sqrt(2.*Ti_par), Ti_par);
    mpi_printf(comm, "v_e,perp = %g [c] T_e,perp = %g\n", sqrt(2*Te_perp / me), Te_perp);
    mpi_printf(comm, "v_e,par  = %g [c] T_e,par = %g\n", sqrt(2.*Te_par / me), Te_par);
    mpi_printf(comm, "\n");
    mpi_printf(comm, "beta_e_par_  = %g [c] mi_over_me_ = %g\n", beta_e_par_, mi_over_me_);   
    mpi_printf(comm, " debye_length_ = %g [c] mi_over_me_ = %g\n", debye_length_, mi_over_me_);   
    mpi_printf(comm, "\n");
    p_.nmax = 3201;  // Ques_jeff How long should I run this?
    p_.cfl = 0.98;

    // -- setup particle kinds
    Grid_t::Kinds kinds = {{ -1., me, "e"}, {1., mi, "i"}, };

    // --- setup domain
    Grid_t::Real3 LL = {17.77, 17.77, 99.74}; // domain size (normalized units, ie, in d_i)
    Int3 gdims = {294, 297, 1648}; // global number of grid points
    Int3 np = {6, 9, 16}; // division into patches
    
    
    auto grid_domain = Grid_t::Domain{gdims, LL, {}, np};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 100; // Ques_jeff 100 particles ??

    double dt = p_.cfl * courant_length(grid_domain);
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles_t{grid()});

    // -- Balance
    balance_interval = -1; //Ques_jeff critical balance on (false -> true)?
    balance_.reset(new Balance_t{balance_interval, .1, false});

    // -- Sort
    // FIXME, needs a way to make sure it gets set?
    sort_interval = 100;

    // -- Collision
    int collision_interval = -1;
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
    marder_interval = -1;
    marder_.reset(new Marder_t(grid(), marder_diffusion, marder_loop, marder_dump));

    // -- output fields
    OutputFieldsCParams outf_params;
    //outf_params.output_fields = "e,h,j,n_1st_single,v_1st_single,T_1st_single";
    outf_params.output_fields = "h";
    outf_params.pfield_step = 400 ;// 1250; // it was every 100
    outf_.reset(new OutputFieldsC{grid(), outf_params});

    // -- output particles
     //output_particle_interval = 2000;
	//	if (output_particle_interval > 0) {
        //	psc_output_particles_set_param_int(outp_,"every_step",
        //                                           (int) (output_particle_interval);
        //                                                }

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

    // do remainder of generic initialization
    init();
  }
  
  // ctor
  // ----------------------------------------------------------------------
  
 
private:
 
  void init_npt(int kind, double crd[3], psc_particle_npt& npt)
  {	

    double B0 = vA_over_c_; //How to define this??
    double Te_par = beta_e_par_ * sqr(B0) / 2.;
    double Te_perp = Te_perp_over_Te_par_ * Te_par;
    double Ti_par = beta_i_par_ * sqr(B0) / 2.;
    double Ti_perp = Ti_perp_over_Ti_par_ * Ti_par;    
    double x = crd[0], y = crd[1], z = crd[2];    
    double Lx = grid().domain.length[0];  
    double Ly = grid().domain.length[1];       
    double Lz = grid().domain.length[2];
    double k_x1, k_x2, k_x3, k_x4, k_x5, k_x6, k_x7, k_x8, k_x9, k_x10, k_x11, k_x12, k_x13, k_x14, k_x15, k_x16;
    double k_y1, k_y2, k_y3, k_y4, k_y5, k_y6, k_y7, k_y8, k_y9, k_y10, k_y11, k_y12, k_y13, k_y14, k_y15, k_y16;
    double k_z1, k_z2, k_z3, k_z4, k_z5, k_z6, k_z7, k_z8, k_z9, k_z10, k_z11, k_z12, k_z13, k_z14, k_z15, k_z16;;
    
/////////////////////////////////////////////////////////////////////////////////////////////    
    double rho_i = sqrt(beta_i_par_)*1.; //check how to set this as di !!
    double sp = 1./3.;  //Spectral index for the alfven wave (1/3 for AW and 2/3 for KAW)
    int p = 1;// fraction of B0 so dB=(1/p)B0
    double crit_fact = 0.1;   //critical balance /normalization coefficient       
    double Amp1, Amp2, Amp3, Amp4, Amp5, Amp6, Amp7, Amp8, Amp9, Amp10, Amp11, Amp12, Amp13, Amp14, Amp15, Amp16;
    double C1; // normalization factor
    int m_per[2]={1,2}; //modes in the perpendicular directions
    int m_par[2]={1,2}; //modes in the parallel direction acording to critical balance    
    double Thephase1[8]={0.401*2.*M_PI, 0.566*2.*M_PI, 0.714*2.*M_PI, 0.729*2.*M_PI, 0.101*2.*M_PI, 0.149*2.*M_PI, 0.584*2.*M_PI, 0.293*2.*M_PI};
    double Thephase2[8]={0.861*2.*M_PI, 0.010*2.*M_PI, 0.519*2.*M_PI, 0.736*2.*M_PI, 0.204*2.*M_PI, 0.848*2.*M_PI, 0.282*2.*M_PI, 0.954*2.*M_PI};
    double phase1, phase2, phase3, phase4, phase5, phase6, phase7, phase8, phase9, phase10, phase11, phase12, phase13, phase14, phase15, phase16 ;
    double mn_per1, mn_per2, mn_per3, mn_per4, mn_per5, mn_per6, mn_per7, mn_per8, mn_per9, mn_per10, mn_per11, mn_per12, mn_per13, mn_per14, mn_per15, mn_per16;  
    double mn_par1, mn_par2, mn_par3, mn_par4, mn_par5, mn_par6, mn_par7, mn_par8, mn_par9, mn_par10, mn_par11, mn_par12, mn_par13, mn_par14, mn_par15, mn_par16;
    mn_per1 = m_per[0] ; mn_per2 = m_per[0] ; mn_per3 = m_per[0] ; mn_per4 = m_per[0] ; mn_per5 = m_per[0] ; mn_per6 = m_per[0] ; mn_per7 = m_per[0] ; mn_per8 = m_per[0] ; 
    mn_per9 = m_per[1] ; mn_per10 = m_per[1] ; mn_per11 = m_per[1] ; mn_per12 = m_per[1] ; mn_per13 = m_per[1] ; mn_per14 = m_per[1] ; mn_per15 = m_per[1] ; mn_per16 = m_per[1] ; 
    mn_par1 = m_par[0] ; mn_par2 = m_par[0] ; mn_par3 = m_par[0] ; mn_par4 = m_par[0] ; mn_par5 = m_par[0] ; mn_par6 = m_par[0] ; mn_par7 = m_par[0] ; mn_par8 = m_par[0] ;
    mn_par9 = m_par[1] ; mn_par10 = m_par[1] ; mn_par11 = m_par[1] ; mn_par12 = m_par[1] ; mn_par13 = m_par[1] ; mn_par14 = m_par[1] ; mn_par15 = m_par[1] ; mn_par16 = m_par[1] ;
    phase1 = Thephase1[0]; phase2 = Thephase1[1]; phase3 = Thephase1[2]; phase4 = Thephase1[3]; phase5 = Thephase1[4]; phase6 = Thephase1[5]; phase7 = Thephase1[6]; phase8 = Thephase1[7];    
    phase9 = Thephase2[0]; phase10 = Thephase2[1]; phase11 = Thephase2[2]; phase12 = Thephase2[3]; phase13 = Thephase2[4]; phase14 = Thephase2[5]; phase15 = Thephase2[6]; phase16 = Thephase2[7];

    double phi1=0.;
    double phi2=M_PI/2.;
    double phi3=M_PI;    
    double phi4=3.*M_PI/2.;
    double phi5=M_PI/4.;
    double phi6=3.*M_PI/4.;
    double phi7=5.*M_PI/4.;    
    double phi8=7.*M_PI/4.;
    double dB_ax1, dB_ax2, dB_ax3, dB_ax4, dB_ax5, dB_ax6, dB_ax7, dB_ax8, dB_ax9, dB_ax10, dB_ax11, dB_ax12, dB_ax13, dB_ax14, dB_ax15, dB_ax16 ;
    double dB_ay1, dB_ay2, dB_ay3, dB_ay4, dB_ay5, dB_ay6, dB_ay7, dB_ay8, dB_ay9, dB_ay10, dB_ay11, dB_ay12, dB_ay13, dB_ay14, dB_ay15, dB_ay16 ;
    double dB_axT, dB_ayT, dv_axT, dv_ayT; 
    double k_a_per1, k_a_per2, k_a_per3, k_a_per4, k_a_per5, k_a_per6, k_a_per7, k_a_per8, k_a_per9, k_a_per10, k_a_per11, k_a_per12, k_a_per13, k_a_per14, k_a_per15, k_a_per16;
    double k_a_par1, k_a_par2, k_a_par3, k_a_par4, k_a_par5, k_a_par6, k_a_par7, k_a_par8, k_a_par9, k_a_par10, k_a_par11, k_a_par12, k_a_par13, k_a_par14, k_a_par15, k_a_par16;
    
    k_x1= mn_per1*(2. * M_PI /Lx) * (1.) ; k_y1 = mn_per1*(2. * M_PI /Ly)*(0.) ; k_z1 = mn_par1*(2. * M_PI /Lz);  
    k_x2= mn_per2*(2. * M_PI /Lx) * (0.) ; k_y2 = mn_per2*(2. * M_PI /Ly)*(1.) ; k_z2 = mn_par2*(2. * M_PI /Lz);  
    k_x3= mn_per3*(2. * M_PI /Lx) * (-1.) ; k_y3 = mn_per3*(2. * M_PI /Ly)*(0.) ; k_z3 = mn_par3*(2. * M_PI /Lz);  
    k_x4= mn_per4*(2. * M_PI /Lx) * (0.) ; k_y4 = mn_per4*(2. * M_PI /Ly)*(-1.) ; k_z4 = mn_par4*(2. * M_PI /Lz);      
    k_x5= mn_per5*(2. * M_PI /Lx) * (1.) ; k_y5 = mn_per5*(2. * M_PI /Ly)*(1.) ; k_z5 = mn_par5*(2. * M_PI /Lz);  
    k_x6= mn_per6*(2. * M_PI /Lx) * (-1.) ; k_y6 = mn_per6*(2. * M_PI /Ly)*(1.) ; k_z6 = mn_par6*(2. * M_PI /Lz);  
    k_x7= mn_per7*(2. * M_PI /Lx) * (-1.) ; k_y7 = mn_per7*(2. * M_PI /Ly)*(-1.) ; k_z7 = mn_par7*(2. * M_PI /Lz);  
    k_x8= mn_per8*(2. * M_PI /Lx) * (1.) ; k_y8 = mn_per8*(2. * M_PI /Ly)*(-1.) ; k_z8 = mn_par8*(2. * M_PI /Lz);
    k_x9= mn_per9*(2. * M_PI /Lx) * (1.) ; k_y9 = mn_per9*(2. * M_PI /Ly)*(0.) ; k_z9 = mn_par9*(2. * M_PI /Lz);
    k_x10= mn_per10*(2. * M_PI /Lx) * (0.) ; k_y10 = mn_per10*(2. * M_PI /Ly)*(1.) ; k_z10 = mn_par10*(2. * M_PI /Lz);
    k_x11= mn_per11*(2. * M_PI /Lx) * (-1.) ; k_y11 = mn_per11*(2. * M_PI /Ly)*(0.) ; k_z11 = mn_par11*(2. * M_PI /Lz);
    k_x12= mn_per12*(2. * M_PI /Lx) * (0.) ; k_y12 = mn_per12*(2. * M_PI /Ly)*(-1.) ; k_z12 = mn_par12*(2. * M_PI /Lz);
    k_x13= mn_per13*(2. * M_PI /Lx) * (1.) ; k_y13 = mn_per13*(2. * M_PI /Ly)*(1.) ; k_z13 = mn_par13*(2. * M_PI /Lz);
    k_x14= mn_per14*(2. * M_PI /Lx) * (-1.) ; k_y14 = mn_per14*(2. * M_PI /Ly)*(1.) ; k_z14 = mn_par14*(2. * M_PI /Lz);
    k_x15= mn_per15*(2. * M_PI /Lx) * (-1.) ; k_y15 = mn_per15*(2. * M_PI /Ly)*(-1.) ; k_z15 = mn_par15*(2. * M_PI /Lz);
    k_x16= mn_per16*(2. * M_PI /Lx) * (1.) ; k_y16 = mn_per16*(2. * M_PI /Ly)*(-1.) ; k_z16 = mn_par16*(2. * M_PI /Lz);

    k_a_per1 = sqrt(k_x1*k_x1 + k_y1*k_y1); k_a_per2 = sqrt(k_x2*k_x2 + k_y2*k_y2); k_a_per3 = sqrt(k_x3*k_x3 + k_y3*k_y3); k_a_per4 = sqrt(k_x4*k_x4 + k_y4*k_y4);
    k_a_per5 = sqrt(k_x5*k_x5 + k_y5*k_y5); k_a_per6 = sqrt(k_x6*k_x6 + k_y6*k_y6); k_a_per7 = sqrt(k_x7*k_x7 + k_y7*k_y7); k_a_per8 = sqrt(k_x8*k_x8 + k_y8*k_y8);
    k_a_per9 = sqrt(k_x9*k_x9 + k_y9*k_y9); k_a_per10 = sqrt(k_x10*k_x10 + k_y10*k_y10); k_a_per11 = sqrt(k_x11*k_x11 + k_y11*k_y11); k_a_per12 = sqrt(k_x12*k_x12 + k_y12*k_y12);
    k_a_per13 = sqrt(k_x13*k_x13 + k_y13*k_y13); k_a_per14 = sqrt(k_x14*k_x14 + k_y14*k_y14); k_a_per15 = sqrt(k_x15*k_x15 + k_y15*k_y15); k_a_per16 = sqrt(k_x16*k_x16 + k_y16*k_y16); 
    k_a_par1 = k_z1; k_a_par2 = k_z2; k_a_par3 = k_z3; k_a_par4 = k_z4; k_a_par5 = k_z5; k_a_par6 = k_z6; k_a_par7 = k_z7; k_a_par8 = k_z8;  
    k_a_par9 = k_z9; k_a_par10 = k_z10; k_a_par11 = k_z11; k_a_par12 = k_z12; k_a_par13 = k_z13; k_a_par14 = k_z14; k_a_par15 = k_z15; k_a_par16 = k_z16; 

    Amp1 = B0 * crit_fact *pow((k_a_per1), -sp) ; Amp2 = B0 * crit_fact *pow((k_a_per2), -sp) ; Amp3 = B0 * crit_fact *pow((k_a_per3), -sp) ; Amp4 = B0 * crit_fact *pow((k_a_per4), -sp) ; 
    Amp5 = B0 * crit_fact *pow((k_a_per5), -sp) ; Amp6 = B0 * crit_fact *pow((k_a_per6), -sp) ; Amp7 = B0 * crit_fact *pow((k_a_per7), -sp) ; Amp8 = B0 * crit_fact *pow((k_a_per8), -sp) ;
    Amp9 = B0 * crit_fact *pow((k_a_per9), -sp) ; Amp10 = B0 * crit_fact *pow((k_a_per10), -sp) ; Amp11 = B0 * crit_fact *pow((k_a_per11), -sp) ; Amp12 = B0 * crit_fact *pow((k_a_per12), -sp) ; 
    Amp13 = B0 * crit_fact *pow((k_a_per13), -sp) ; Amp14 = B0 * crit_fact *pow((k_a_per14), -sp) ; Amp15 = B0 * crit_fact *pow((k_a_per15), -sp) ; Amp16 = B0 * crit_fact *pow((k_a_per16), -sp) ;

    dB_ax1 = -Amp1 * cos ( k_x1* x + k_y1 * y + k_z1* z + phase1 )*sin(phi1) ; dB_ax2 = -Amp2 * cos ( k_x2* x + k_y2 * y - k_z2* z + phase2 )*sin(phi2) ;
    dB_ax3 = -Amp3 * cos ( k_x3* x + k_y3 * y + k_z3* z + phase3 )*sin(phi3) ; dB_ax4 = -Amp4 * cos ( k_x4* x + k_y4 * y - k_z4* z + phase4 )*sin(phi4) ;
    dB_ax5 = -Amp5 * cos ( k_x5* x + k_y5 * y + k_z5* z + phase5 )*sin(phi5) ; dB_ax6 = -Amp6 * cos ( k_x6* x + k_y6 * y - k_z6* z + phase6 )*sin(phi6) ;
    dB_ax7 = -Amp7 * cos ( k_x7* x + k_y7 * y + k_z7* z + phase7 )*sin(phi7) ; dB_ax8 = -Amp8 * cos ( k_x8* x + k_y8 * y - k_z8* z + phase8 )*sin(phi8) ;
    dB_ax9 = -Amp9 * cos ( k_x9* x + k_y9 * y + k_z9* z + phase9 )*sin(phi1) ; dB_ax10 = -Amp10 * cos ( k_x10* x + k_y10 * y - k_z10* z + phase10 )*sin(phi2) ;
    dB_ax11 = -Amp11 * cos ( k_x11* x + k_y11 * y + k_z11* z + phase11 )*sin(phi3) ; dB_ax12 = -Amp12 * cos ( k_x12* x + k_y12 * y - k_z12* z + phase12 )*sin(phi4) ;
    dB_ax13 = -Amp13 * cos ( k_x13* x + k_y13 * y + k_z13* z + phase13 )*sin(phi5) ; dB_ax14 = -Amp14 * cos ( k_x14* x + k_y14 * y - k_z14* z + phase14 )*sin(phi6) ;
    dB_ax15 = -Amp15 * cos ( k_x15* x + k_y15 * y + k_z15* z + phase15 )*sin(phi7) ; dB_ax16 = -Amp16 * cos ( k_x16* x + k_y16 * y - k_z16* z + phase16 )*sin(phi8) ;
    
    dB_ay1 =  Amp1 * cos ( k_x1* x + k_y1 * y + k_z1* z + phase1 )*cos(phi1) ; dB_ay2 =  Amp2 * cos ( k_x2* x + k_y2 * y - k_z2* z + phase2 )*cos(phi2) ;    
    dB_ay3 =  Amp3 * cos ( k_x3* x + k_y3 * y + k_z3* z + phase3 )*cos(phi3) ; dB_ay4 =  Amp4 * cos ( k_x4* x + k_y4 * y - k_z4* z + phase4 )*cos(phi4) ;    
    dB_ay5 =  Amp5 * cos ( k_x5* x + k_y5 * y + k_z5* z + phase5 )*cos(phi5) ; dB_ay6 =  Amp6 * cos ( k_x6* x + k_y6 * y - k_z6* z + phase6 )*cos(phi6) ;    
    dB_ay7 =  Amp7 * cos ( k_x7* x + k_y7 * y + k_z7* z + phase7 )*cos(phi7) ; dB_ay8 =  Amp8 * cos ( k_x8* x + k_y8 * y - k_z8* z + phase8 )*cos(phi8) ;    
    dB_ax9 =  Amp9 * cos ( k_x9* x + k_y9 * y + k_z9* z + phase9 )*cos(phi1) ; dB_ax10 = Amp10 * cos ( k_x10* x + k_y10 * y - k_z10* z + phase10 )*cos(phi2) ;
    dB_ax11 = Amp11 * cos ( k_x11* x + k_y11 * y + k_z11* z + phase11 )*cos(phi3) ; dB_ax12 = Amp12 * cos ( k_x12* x + k_y12 * y - k_z12* z + phase12 )*cos(phi4) ;
    dB_ax13 = Amp13 * cos ( k_x13* x + k_y13 * y + k_z13* z + phase13 )*cos(phi5) ; dB_ax14 = Amp14 * cos ( k_x14* x + k_y14 * y - k_z14* z + phase14 )*cos(phi6) ;
    dB_ax15 = Amp15 * cos ( k_x15* x + k_y15 * y + k_z15* z + phase15 )*cos(phi7) ; dB_ax16 = Amp16 * cos ( k_x16* x + k_y16 * y - k_z16* z + phase16 )*cos(phi8) ;    

//magnetic field
    dB_axT = dB_ax1 + dB_ax2 + dB_ax3 + dB_ax4 + dB_ax5 + dB_ax6 + dB_ax7 + dB_ax8 + dB_ax9 + dB_ax10 + dB_ax11 + dB_ax12 + dB_ax13 + dB_ax14 + dB_ax15 + dB_ax16;
    dB_ayT = dB_ay1 + dB_ay2 + dB_ay3 + dB_ay4 + dB_ay5 + dB_ay6 + dB_ay7 + dB_ay8 + dB_ay9 + dB_ay10 + dB_ay11 + dB_ay12 + dB_ay13 + dB_ay14 + dB_ay15 + dB_ay16;
//velocities
    dv_axT = -dB_ax1 + dB_ax2 -dB_ax3 + dB_ax4 -dB_ax5 + dB_ax6 -dB_ax7 + dB_ax8 -dB_ax9 + dB_ax10 -dB_ax11 + dB_ax12 -dB_ax13 + dB_ax14 -dB_ax15 + dB_ax16;
    dv_ayT = -dB_ay1 + dB_ay2 -dB_ay3 + dB_ay4 -dB_ay5 + dB_ay6 -dB_ay7 + dB_ay8 -dB_ay9 + dB_ay10 -dB_ay11 + dB_ay12 -dB_ay13 + dB_ay14 -dB_ay15 + dB_ay16;

    C1 = B0/sqrt(sqr(Amp1)+sqr(Amp2)+sqr(Amp3)+sqr(Amp4)+sqr(Amp5)+sqr(Amp6)+sqr(Amp7)+sqr(Amp8)+sqr(Amp9)+sqr(Amp10)+sqr(Amp11)+sqr(Amp12)+sqr(Amp13)+sqr(Amp14)+sqr(Amp15)+sqr(Amp16)); //    
/////////////////////////////////////////////////////////////////////////////////////////////
    
    switch (kind) {
    case TurAW_ELECTRON:
      npt.n = 1.;
      
      npt.T[0] = Te_perp; //Should I change the temperature?
      npt.T[1] = Te_perp;
      npt.T[2] = Te_par;
      
      // Set velocities for first wave:
      npt.p[0] = C1* dv_axT; //this - is the one of the direction of propagation
      npt.p[1] = C1* dv_ayT;
      npt.p[2] = 0.;
      
      break;
    case TurAW_ION:
      npt.n = 1.;
      
      npt.T[0] = Ti_perp;
      npt.T[1] = Ti_perp;
      npt.T[2] = Ti_par;
      
      // Set velocities for first wave:
      npt.p[0] = C1* dv_axT; //this - is the one of the direction of propagation
      npt.p[1] = C1* dv_ayT;
      npt.p[2] = 0.;
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
    double B0 = vA_over_c_; //How to define this?? 

    SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
    double x = crd[0], y = crd[1], z = crd[2];    

    double Lx = grid().domain.length[0];  
    double Ly = grid().domain.length[1];       
    double Lz = grid().domain.length[2];
    double k_x1, k_x2, k_x3, k_x4, k_x5, k_x6, k_x7, k_x8, k_x9, k_x10, k_x11, k_x12, k_x13, k_x14, k_x15, k_x16;
    double k_y1, k_y2, k_y3, k_y4, k_y5, k_y6, k_y7, k_y8, k_y9, k_y10, k_y11, k_y12, k_y13, k_y14, k_y15, k_y16;
    double k_z1, k_z2, k_z3, k_z4, k_z5, k_z6, k_z7, k_z8, k_z9, k_z10, k_z11, k_z12, k_z13, k_z14, k_z15, k_z16;;
    
/////////////////////////////////////////////////////////////////////////////////////////////    
    double rho_i = sqrt(beta_i_par_)*1.; //check how to set this as di !!
    double sp = 1./3.;  //Spectral index for the alfven wave (1/3 for AW and 2/3 for KAW)
    int p = 1;// fraction of B0 so dB=(1/p)B0
    double crit_fact = 0.1;   //critical balance /normalization coefficient       
    double Amp1, Amp2, Amp3, Amp4, Amp5, Amp6, Amp7, Amp8, Amp9, Amp10, Amp11, Amp12, Amp13, Amp14, Amp15, Amp16;
    double C1; // normalization factor
    int m_per[2]={1,2}; //modes in the perpendicular directions
    int m_par[2]={1,2}; //modes in the parallel direction acording to critical balance    
    double Thephase1[8]={0.401*2.*M_PI, 0.566*2.*M_PI, 0.714*2.*M_PI, 0.729*2.*M_PI, 0.101*2.*M_PI, 0.149*2.*M_PI, 0.584*2.*M_PI, 0.293*2.*M_PI};
    double Thephase2[8]={0.861*2.*M_PI, 0.010*2.*M_PI, 0.519*2.*M_PI, 0.736*2.*M_PI, 0.204*2.*M_PI, 0.848*2.*M_PI, 0.282*2.*M_PI, 0.954*2.*M_PI};
    double phase1, phase2, phase3, phase4, phase5, phase6, phase7, phase8, phase9, phase10, phase11, phase12, phase13, phase14, phase15, phase16 ;
    double mn_per1, mn_per2, mn_per3, mn_per4, mn_per5, mn_per6, mn_per7, mn_per8, mn_per9, mn_per10, mn_per11, mn_per12, mn_per13, mn_per14, mn_per15, mn_per16;  
    double mn_par1, mn_par2, mn_par3, mn_par4, mn_par5, mn_par6, mn_par7, mn_par8, mn_par9, mn_par10, mn_par11, mn_par12, mn_par13, mn_par14, mn_par15, mn_par16;
    mn_per1 = m_per[0] ; mn_per2 = m_per[0] ; mn_per3 = m_per[0] ; mn_per4 = m_per[0] ; mn_per5 = m_per[0] ; mn_per6 = m_per[0] ; mn_per7 = m_per[0] ; mn_per8 = m_per[0] ; 
    mn_per9 = m_per[1] ; mn_per10 = m_per[1] ; mn_per11 = m_per[1] ; mn_per12 = m_per[1] ; mn_per13 = m_per[1] ; mn_per14 = m_per[1] ; mn_per15 = m_per[1] ; mn_per16 = m_per[1] ; 
    mn_par1 = m_par[0] ; mn_par2 = m_par[0] ; mn_par3 = m_par[0] ; mn_par4 = m_par[0] ; mn_par5 = m_par[0] ; mn_par6 = m_par[0] ; mn_par7 = m_par[0] ; mn_par8 = m_par[0] ;
    mn_par9 = m_par[1] ; mn_par10 = m_par[1] ; mn_par11 = m_par[1] ; mn_par12 = m_par[1] ; mn_par13 = m_par[1] ; mn_par14 = m_par[1] ; mn_par15 = m_par[1] ; mn_par16 = m_par[1] ;
    phase1 = Thephase1[0]; phase2 = Thephase1[1]; phase3 = Thephase1[2]; phase4 = Thephase1[3]; phase5 = Thephase1[4]; phase6 = Thephase1[5]; phase7 = Thephase1[6]; phase8 = Thephase1[7];    
    phase9 = Thephase2[0]; phase10 = Thephase2[1]; phase11 = Thephase2[2]; phase12 = Thephase2[3]; phase13 = Thephase2[4]; phase14 = Thephase2[5]; phase15 = Thephase2[6]; phase16 = Thephase2[7];

    double phi1=0.;
    double phi2=M_PI/2.;
    double phi3=M_PI;    
    double phi4=3.*M_PI/2.;
    double phi5=M_PI/4.;
    double phi6=3.*M_PI/4.;
    double phi7=5.*M_PI/4.;    
    double phi8=7.*M_PI/4.;
    double dB_ax1, dB_ax2, dB_ax3, dB_ax4, dB_ax5, dB_ax6, dB_ax7, dB_ax8, dB_ax9, dB_ax10, dB_ax11, dB_ax12, dB_ax13, dB_ax14, dB_ax15, dB_ax16 ;
    double dB_ay1, dB_ay2, dB_ay3, dB_ay4, dB_ay5, dB_ay6, dB_ay7, dB_ay8, dB_ay9, dB_ay10, dB_ay11, dB_ay12, dB_ay13, dB_ay14, dB_ay15, dB_ay16 ;
    double dB_axT, dB_ayT, dv_axT, dv_ayT; 
    double k_a_per1, k_a_per2, k_a_per3, k_a_per4, k_a_per5, k_a_per6, k_a_per7, k_a_per8, k_a_per9, k_a_per10, k_a_per11, k_a_per12, k_a_per13, k_a_per14, k_a_per15, k_a_per16;
    double k_a_par1, k_a_par2, k_a_par3, k_a_par4, k_a_par5, k_a_par6, k_a_par7, k_a_par8, k_a_par9, k_a_par10, k_a_par11, k_a_par12, k_a_par13, k_a_par14, k_a_par15, k_a_par16;
    
    k_x1= mn_per1*(2. * M_PI /Lx) * (1.) ; k_y1 = mn_per1*(2. * M_PI /Ly)*(0.) ; k_z1 = mn_par1*(2. * M_PI /Lz);  
    k_x2= mn_per2*(2. * M_PI /Lx) * (0.) ; k_y2 = mn_per2*(2. * M_PI /Ly)*(1.) ; k_z2 = mn_par2*(2. * M_PI /Lz);  
    k_x3= mn_per3*(2. * M_PI /Lx) * (-1.) ; k_y3 = mn_per3*(2. * M_PI /Ly)*(0.) ; k_z3 = mn_par3*(2. * M_PI /Lz);  
    k_x4= mn_per4*(2. * M_PI /Lx) * (0.) ; k_y4 = mn_per4*(2. * M_PI /Ly)*(-1.) ; k_z4 = mn_par4*(2. * M_PI /Lz);      
    k_x5= mn_per5*(2. * M_PI /Lx) * (1.) ; k_y5 = mn_per5*(2. * M_PI /Ly)*(1.) ; k_z5 = mn_par5*(2. * M_PI /Lz);  
    k_x6= mn_per6*(2. * M_PI /Lx) * (-1.) ; k_y6 = mn_per6*(2. * M_PI /Ly)*(1.) ; k_z6 = mn_par6*(2. * M_PI /Lz);  
    k_x7= mn_per7*(2. * M_PI /Lx) * (-1.) ; k_y7 = mn_per7*(2. * M_PI /Ly)*(-1.) ; k_z7 = mn_par7*(2. * M_PI /Lz);  
    k_x8= mn_per8*(2. * M_PI /Lx) * (1.) ; k_y8 = mn_per8*(2. * M_PI /Ly)*(-1.) ; k_z8 = mn_par8*(2. * M_PI /Lz);
    k_x9= mn_per9*(2. * M_PI /Lx) * (1.) ; k_y9 = mn_per9*(2. * M_PI /Ly)*(0.) ; k_z9 = mn_par9*(2. * M_PI /Lz);
    k_x10= mn_per10*(2. * M_PI /Lx) * (0.) ; k_y10 = mn_per10*(2. * M_PI /Ly)*(1.) ; k_z10 = mn_par10*(2. * M_PI /Lz);
    k_x11= mn_per11*(2. * M_PI /Lx) * (-1.) ; k_y11 = mn_per11*(2. * M_PI /Ly)*(0.) ; k_z11 = mn_par11*(2. * M_PI /Lz);
    k_x12= mn_per12*(2. * M_PI /Lx) * (0.) ; k_y12 = mn_per12*(2. * M_PI /Ly)*(-1.) ; k_z12 = mn_par12*(2. * M_PI /Lz);
    k_x13= mn_per13*(2. * M_PI /Lx) * (1.) ; k_y13 = mn_per13*(2. * M_PI /Ly)*(1.) ; k_z13 = mn_par13*(2. * M_PI /Lz);
    k_x14= mn_per14*(2. * M_PI /Lx) * (-1.) ; k_y14 = mn_per14*(2. * M_PI /Ly)*(1.) ; k_z14 = mn_par14*(2. * M_PI /Lz);
    k_x15= mn_per15*(2. * M_PI /Lx) * (-1.) ; k_y15 = mn_per15*(2. * M_PI /Ly)*(-1.) ; k_z15 = mn_par15*(2. * M_PI /Lz);
    k_x16= mn_per16*(2. * M_PI /Lx) * (1.) ; k_y16 = mn_per16*(2. * M_PI /Ly)*(-1.) ; k_z16 = mn_par16*(2. * M_PI /Lz);

    k_a_per1 = sqrt(k_x1*k_x1 + k_y1*k_y1); k_a_per2 = sqrt(k_x2*k_x2 + k_y2*k_y2); k_a_per3 = sqrt(k_x3*k_x3 + k_y3*k_y3); k_a_per4 = sqrt(k_x4*k_x4 + k_y4*k_y4);
    k_a_per5 = sqrt(k_x5*k_x5 + k_y5*k_y5); k_a_per6 = sqrt(k_x6*k_x6 + k_y6*k_y6); k_a_per7 = sqrt(k_x7*k_x7 + k_y7*k_y7); k_a_per8 = sqrt(k_x8*k_x8 + k_y8*k_y8);
    k_a_per9 = sqrt(k_x9*k_x9 + k_y9*k_y9); k_a_per10 = sqrt(k_x10*k_x10 + k_y10*k_y10); k_a_per11 = sqrt(k_x11*k_x11 + k_y11*k_y11); k_a_per12 = sqrt(k_x12*k_x12 + k_y12*k_y12);
    k_a_per13 = sqrt(k_x13*k_x13 + k_y13*k_y13); k_a_per14 = sqrt(k_x14*k_x14 + k_y14*k_y14); k_a_per15 = sqrt(k_x15*k_x15 + k_y15*k_y15); k_a_per16 = sqrt(k_x16*k_x16 + k_y16*k_y16); 
    k_a_par1 = k_z1; k_a_par2 = k_z2; k_a_par3 = k_z3; k_a_par4 = k_z4; k_a_par5 = k_z5; k_a_par6 = k_z6; k_a_par7 = k_z7; k_a_par8 = k_z8;  
    k_a_par9 = k_z9; k_a_par10 = k_z10; k_a_par11 = k_z11; k_a_par12 = k_z12; k_a_par13 = k_z13; k_a_par14 = k_z14; k_a_par15 = k_z15; k_a_par16 = k_z16; 

    Amp1 = B0 * crit_fact *pow((k_a_per1), -sp) ; Amp2 = B0 * crit_fact *pow((k_a_per2), -sp) ; Amp3 = B0 * crit_fact *pow((k_a_per3), -sp) ; Amp4 = B0 * crit_fact *pow((k_a_per4), -sp) ; 
    Amp5 = B0 * crit_fact *pow((k_a_per5), -sp) ; Amp6 = B0 * crit_fact *pow((k_a_per6), -sp) ; Amp7 = B0 * crit_fact *pow((k_a_per7), -sp) ; Amp8 = B0 * crit_fact *pow((k_a_per8), -sp) ;
    Amp9 = B0 * crit_fact *pow((k_a_per9), -sp) ; Amp10 = B0 * crit_fact *pow((k_a_per10), -sp) ; Amp11 = B0 * crit_fact *pow((k_a_per11), -sp) ; Amp12 = B0 * crit_fact *pow((k_a_per12), -sp) ; 
    Amp13 = B0 * crit_fact *pow((k_a_per13), -sp) ; Amp14 = B0 * crit_fact *pow((k_a_per14), -sp) ; Amp15 = B0 * crit_fact *pow((k_a_per15), -sp) ; Amp16 = B0 * crit_fact *pow((k_a_per16), -sp) ;

    dB_ax1 = -Amp1 * cos ( k_x1* x + k_y1 * y + k_z1* z + phase1 )*sin(phi1) ; dB_ax2 = -Amp2 * cos ( k_x2* x + k_y2 * y - k_z2* z + phase2 )*sin(phi2) ;
    dB_ax3 = -Amp3 * cos ( k_x3* x + k_y3 * y + k_z3* z + phase3 )*sin(phi3) ; dB_ax4 = -Amp4 * cos ( k_x4* x + k_y4 * y - k_z4* z + phase4 )*sin(phi4) ;
    dB_ax5 = -Amp5 * cos ( k_x5* x + k_y5 * y + k_z5* z + phase5 )*sin(phi5) ; dB_ax6 = -Amp6 * cos ( k_x6* x + k_y6 * y - k_z6* z + phase6 )*sin(phi6) ;
    dB_ax7 = -Amp7 * cos ( k_x7* x + k_y7 * y + k_z7* z + phase7 )*sin(phi7) ; dB_ax8 = -Amp8 * cos ( k_x8* x + k_y8 * y - k_z8* z + phase8 )*sin(phi8) ;
    dB_ax9 = -Amp9 * cos ( k_x9* x + k_y9 * y + k_z9* z + phase9 )*sin(phi1) ; dB_ax10 = -Amp10 * cos ( k_x10* x + k_y10 * y - k_z10* z + phase10 )*sin(phi2) ;
    dB_ax11 = -Amp11 * cos ( k_x11* x + k_y11 * y + k_z11* z + phase11 )*sin(phi3) ; dB_ax12 = -Amp12 * cos ( k_x12* x + k_y12 * y - k_z12* z + phase12 )*sin(phi4) ;
    dB_ax13 = -Amp13 * cos ( k_x13* x + k_y13 * y + k_z13* z + phase13 )*sin(phi5) ; dB_ax14 = -Amp14 * cos ( k_x14* x + k_y14 * y - k_z14* z + phase14 )*sin(phi6) ;
    dB_ax15 = -Amp15 * cos ( k_x15* x + k_y15 * y + k_z15* z + phase15 )*sin(phi7) ; dB_ax16 = -Amp16 * cos ( k_x16* x + k_y16 * y - k_z16* z + phase16 )*sin(phi8) ;
    
    dB_ay1 =  Amp1 * cos ( k_x1* x + k_y1 * y + k_z1* z + phase1 )*cos(phi1) ; dB_ay2 =  Amp2 * cos ( k_x2* x + k_y2 * y - k_z2* z + phase2 )*cos(phi2) ;    
    dB_ay3 =  Amp3 * cos ( k_x3* x + k_y3 * y + k_z3* z + phase3 )*cos(phi3) ; dB_ay4 =  Amp4 * cos ( k_x4* x + k_y4 * y - k_z4* z + phase4 )*cos(phi4) ;    
    dB_ay5 =  Amp5 * cos ( k_x5* x + k_y5 * y + k_z5* z + phase5 )*cos(phi5) ; dB_ay6 =  Amp6 * cos ( k_x6* x + k_y6 * y - k_z6* z + phase6 )*cos(phi6) ;    
    dB_ay7 =  Amp7 * cos ( k_x7* x + k_y7 * y + k_z7* z + phase7 )*cos(phi7) ; dB_ay8 =  Amp8 * cos ( k_x8* x + k_y8 * y - k_z8* z + phase8 )*cos(phi8) ;    
    dB_ax9 =  Amp9 * cos ( k_x9* x + k_y9 * y + k_z9* z + phase9 )*cos(phi1) ; dB_ax10 = Amp10 * cos ( k_x10* x + k_y10 * y - k_z10* z + phase10 )*cos(phi2) ;
    dB_ax11 = Amp11 * cos ( k_x11* x + k_y11 * y + k_z11* z + phase11 )*cos(phi3) ; dB_ax12 = Amp12 * cos ( k_x12* x + k_y12 * y - k_z12* z + phase12 )*cos(phi4) ;
    dB_ax13 = Amp13 * cos ( k_x13* x + k_y13 * y + k_z13* z + phase13 )*cos(phi5) ; dB_ax14 = Amp14 * cos ( k_x14* x + k_y14 * y - k_z14* z + phase14 )*cos(phi6) ;
    dB_ax15 = Amp15 * cos ( k_x15* x + k_y15 * y + k_z15* z + phase15 )*cos(phi7) ; dB_ax16 = Amp16 * cos ( k_x16* x + k_y16 * y - k_z16* z + phase16 )*cos(phi8) ;    

//magnetic field
    dB_axT = dB_ax1 + dB_ax2 + dB_ax3 + dB_ax4 + dB_ax5 + dB_ax6 + dB_ax7 + dB_ax8 + dB_ax9 + dB_ax10 + dB_ax11 + dB_ax12 + dB_ax13 + dB_ax14 + dB_ax15 + dB_ax16;
    dB_ayT = dB_ay1 + dB_ay2 + dB_ay3 + dB_ay4 + dB_ay5 + dB_ay6 + dB_ay7 + dB_ay8 + dB_ay9 + dB_ay10 + dB_ay11 + dB_ay12 + dB_ay13 + dB_ay14 + dB_ay15 + dB_ay16;
//velocities
    dv_axT = -dB_ax1 + dB_ax2 -dB_ax3 + dB_ax4 -dB_ax5 + dB_ax6 -dB_ax7 + dB_ax8 -dB_ax9 + dB_ax10 -dB_ax11 + dB_ax12 -dB_ax13 + dB_ax14 -dB_ax15 + dB_ax16;
    dv_ayT = -dB_ay1 + dB_ay2 -dB_ay3 + dB_ay4 -dB_ay5 + dB_ay6 -dB_ay7 + dB_ay8 -dB_ay9 + dB_ay10 -dB_ay11 + dB_ay12 -dB_ay13 + dB_ay14 -dB_ay15 + dB_ay16;

    C1 = B0/sqrt(sqr(Amp1)+sqr(Amp2)+sqr(Amp3)+sqr(Amp4)+sqr(Amp5)+sqr(Amp6)+sqr(Amp7)+sqr(Amp8)+sqr(Amp9)+sqr(Amp10)+sqr(Amp11)+sqr(Amp12)+sqr(Amp13)+sqr(Amp14)+sqr(Amp15)+sqr(Amp16)); //    
/////////////////////////////////////////////////////////////////////////////////////////////

	  switch (m) {
	  case HX: return  C1*dB_axT;
	  case HY: return  C1*dB_ayT;
	  case HZ: return B0;
	  
	  default: return 0.;
	  }
      });
  }
  
  // Definition of variables
  /////////////////////////////////////////////////////////////////////////
  
private:
    
    double mi_over_me_;
    double vA_over_c_;
    double amplitude_;
    double beta_e_par_;
    double beta_i_par_;
    double Ti_perp_over_Ti_par_;
    double Te_perp_over_Te_par_;
    double  output_particle_interval;
  
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);

  auto psc = new PscTurAW;

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}

