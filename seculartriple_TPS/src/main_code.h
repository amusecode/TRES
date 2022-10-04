#include <math.h>
#include <cstdlib>

/*********
 * CVODE *
 *********/
 
/* general */
#include "cvode/cvode.h"					/* prototypes for CVODE fcts., consts. */
#include "cvode/nvector_serial.h"			/* serial N_Vector types, fcts., macros */
#include "cvode/sundials_dense.h"			/* definitions DlsMat DENSE_ELEM */
#include "cvode/sundials_types.h"			/* definition of type realtype */
#include "cvode/sundials_math.h"

/* linear solvers */
#include "cvode/cvode_dense.h"				
#include "cvode/cvode_diag.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"
#include "cvode/cvode_sptfqmr.h"

/* macros */
#define Ith(v,i)    NV_Ith_S(v,i-1)       		/* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 		/* IJth numbers rows,cols 1..NEQ */
#ifndef max
	#define max(a,b) ( ((a) > (b)) ? (a):(b) )
#endif
#ifndef min
    #define min(a,b) ( ((a) < (b)) ? (a):(b) )
#endif

/*********
 * other *
 * *******/
 
#include "helper_routines.h"
#include "tidal_friction_parameters.h"
#include "ODE_system.h"

#define CONST_G			    (double)	3.94852492465e+13
#define CONST_G_P2          (double)    CONST_G*CONST_G
#define CONST_G_P3          (double)    CONST_G*CONST_G_P2
#define CONST_C_LIGHT		(double)	63239726386.8
#define CONST_C_LIGHT_P2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_P4	(double)	CONST_C_LIGHT_P2*CONST_C_LIGHT_P2
#define CONST_C_LIGHT_P5	(double)	CONST_C_LIGHT_P4*CONST_C_LIGHT
#define CONST_MSUN          (double)    1.0
#define CONST_R_SUN         (double)    0.00464913034382
#define CONST_L_SUN         (double)    2.71040410975e+14

#ifdef IGNORE
#define CONST_G			    (double)	6.674279999999998e-08
#define CONST_G_P2          (double)    CONST_G*CONST_G
#define CONST_G_P3          (double)    CONST_G*CONST_G_P2
#define CONST_C_LIGHT		(double)	2.99792458e10
#define CONST_C_LIGHT_P2	(double)	CONST_C_LIGHT*CONST_C_LIGHT
#define CONST_C_LIGHT_P4	(double)	CONST_C_LIGHT_P2*CONST_C_LIGHT_P2
#define CONST_C_LIGHT_P5	(double)	CONST_C_LIGHT_P4*CONST_C_LIGHT
#define CONST_MSUN          (double)    1.98892e33
#define CONST_R_SUN         (double)    6.955e10
#define CONST_L_SUN         (double)    3.838999999999999e33
#endif

#define c_1div2             (double)    1.0/2.0
#define c_1div3             (double)    1.0/3.0
#define c_1div4             (double)    1.0/4.0
#define c_1div5             (double)    1.0/5.0
#define c_1div6             (double)    1.0/6.0
#define c_1div7             (double)    1.0/7.0
#define c_1div8             (double)    1.0/8.0
#define c_1div16            (double)    1.0/16.0
#define c_2div3             (double)    2.0/3.0
#define c_2div5             (double)    2.0/5.0
#define c_3div2             (double)    3.0/2.0
#define c_3div5             (double)    3.0/5.0
#define c_3div8             (double)    3.0/8.0
#define c_5div2             (double)    5.0/2.0
#define c_5div16            (double)    5.0/16.0
#define c_5div64            (double)    5.0/64.0
#define c_8div5             (double)    8.0/5.0
#define c_9div16            (double)    9.0/16.0
#define c_11div18           (double)    11.0/18.0
#define c_15div2            (double)    15.0/2.0
#define c_15div4            (double)    15.0/4.0
#define c_15div8            (double)    15.0/8.0
#define c_15div16           (double)    15.0/16.0
#define c_25div64           (double)    25.0/64.0
#define c_31div2            (double)    31.0/2.0
#define c_37div96           (double)    37.0/96.0
#define c_45div8            (double)    45.0/8.0
#define c_64div5            (double)    64.0/5.0
#define c_73div24           (double)    73.0/24.0
#define c_121div304         (double)    121.0/304.0
#define c_185div16          (double)    185.0/16.0
#define c_255div8           (double)    255.0/8.0
#define c_304div15          (double)    304.0/15.0

typedef struct {
    double global_time_step; // the global time-step
    
    int stellar_type1,stellar_type2,stellar_type3;
	double m1,m2,m3;
    double m1_convective_envelope,m2_convective_envelope,m3_convective_envelope;
    double R1,R2,R3;
    double R1_convective_envelope,R2_convective_envelope,R3_convective_envelope;
    
    bool include_quadrupole_terms,include_octupole_terms;
    bool include_1PN_inner_terms,include_1PN_outer_terms,include_1PN_inner_outer_terms,include_25PN_inner_terms,include_25PN_outer_terms;
    bool include_inner_tidal_terms,include_outer_tidal_terms;
    bool ignore_tertiary;
    bool include_inner_wind_terms,include_outer_wind_terms;
    bool include_magnetic_braking_terms;
    bool include_spin_radius_mass_coupling_terms_star1;
    bool include_spin_radius_mass_coupling_terms_star2;
    bool include_spin_radius_mass_coupling_terms_star3;
    bool include_inner_RLOF_terms,include_outer_RLOF_terms;
    bool check_for_dynamical_stability;
    bool check_for_inner_collision,check_for_outer_collision;
    bool check_for_inner_RLOF,check_for_outer_RLOF;
    bool include_linear_mass_change,include_linear_radius_change;
    int roche_radius_specification;
    int stability_limit_specification;
    
    
    /* addition Jan 2017 */
    bool check_for_semisecular_regime;
    double check_for_semisecular_regime_parameter;
    
    bool stop_after_error_bool;

    double AMC_star1,AMC_star2,AMC_star3;
    double luminosity_star1,luminosity_star2,luminosity_star3;
    double gyration_radius_star1,gyration_radius_star2,gyration_radius_star3;
    double k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3;
    double moment_of_inertia_star1,moment_of_inertia_star2,moment_of_inertia_star3;
    double moment_of_inertia_dot_star1,moment_of_inertia_dot_star2,moment_of_inertia_dot_star3;
    
    bool star1_is_donor,star2_is_donor,star3_is_donor;
    double threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero;
    double threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero;
    double wind_mass_loss_rate_star1,wind_mass_loss_rate_star2,wind_mass_loss_rate_star3;
    double time_derivative_of_radius_star1,time_derivative_of_radius_star2,time_derivative_of_radius_star3;
    double inner_mass_transfer_rate,outer_mass_transfer_rate;
    double inner_accretion_efficiency_wind_child1_to_child2,inner_accretion_efficiency_wind_child2_to_child1;
    double outer_accretion_efficiency_wind_child1_to_child2,outer_accretion_efficiency_wind_child2_to_child1;
    double inner_accretion_efficiency_mass_transfer,outer_accretion_efficiency_mass_transfer;
    double inner_specific_AM_loss_mass_transfer,outer_specific_AM_loss_mass_transfer;
    double inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2,inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1;
    
} *UserData;

int evolve(
    int stellar_type1, int stellar_type2, int stellar_type3,
    double m1, double m2, double m3,
    double m1_convective_envelope, double m2_convective_envelope, double m3_convective_envelope,
    double R1, double R2, double R3,
    double R1_convective_envelope, double R2_convective_envelope, double R3_convective_envelope,    
    double luminosity_star1, double luminosity_star2, double luminosity_star3,
    double spin_angular_frequency1, double spin_angular_frequency2, double spin_angular_frequency3,
    double AMC_star1, double AMC_star2, double AMC_star3,
    double gyration_radius_star1, double gyration_radius_star2, double gyration_radius_star3,
    double moment_of_inertia_star1, double moment_of_inertia_star2, double moment_of_inertia_star3,
    double moment_of_inertia_dot_star1, double moment_of_inertia_dot_star2, double moment_of_inertia_dot_star3,
//    double k_div_T_tides_star1, double k_div_T_tides_star2, double k_div_T_tides_star3,
    double a_in, double a_out,
    double e_in, double e_out,
    double INCL_in, double INCL_out, double AP_in, double AP_out, double LAN_in, double LAN_out,
    bool star1_is_donor, bool star2_is_donor, bool star3_is_donor,
    double wind_mass_loss_rate_star1, double wind_mass_loss_rate_star2, double wind_mass_loss_rate_star3,
    double time_derivative_of_radius_star1, double time_derivative_of_radius_star2, double time_derivative_of_radius_star3,
    double inner_mass_transfer_rate, double outer_mass_transfer_rate,
    double inner_accretion_efficiency_wind_child1_to_child2, double inner_accretion_efficiency_wind_child2_to_child1,
    double outer_accretion_efficiency_wind_child1_to_child2,double outer_accretion_efficiency_wind_child2_to_child1,
    double inner_accretion_efficiency_mass_transfer, double outer_accretion_efficiency_mass_transfer,
    double inner_specific_AM_loss_mass_transfer, double outer_specific_AM_loss_mass_transfer,    
    double inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2, double inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1,
    double t, double global_time_step,
    double * m1_output, double * m2_output, double * m3_output,
    double * R1_output, double * R2_output, double * R3_output,
    double * spin_angular_frequency1_output, double * spin_angular_frequency2_output, double * spin_angular_frequency3_output,
    double * a_in_output, double * a_out_output,
    double * e_in_output, double * e_out_output, double * max_delta_e_output,
    double *INCL_in_output, double *INCL_out_output, double *INCL_in_out_output, double * AP_in_output, double * AP_out_output, double *LAN_in_output, double *LAN_out_output,
    double * t_output,
    int * CVODE_flag, int * root_finding_flag
);
static int check_flag(void *flagvalue, char *funcname, int opt);
void error_handling_function(int error_code, const char *module, const char *function, char *message, void *data_f);
int get_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero(double *value);
int set_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero(double value);
int get_threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero(double *value);
int set_threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero(double value);
int get_equations_of_motion_specification(int *value);
int set_equations_of_motion_specification(int value);
int get_roche_radius_specification(int *value);
int set_roche_radius_specification(int value);
int get_stability_limit_specification(int *value);
int set_stability_limit_specification(int value);
int get_input_precision(double *input_precision_t);
int set_input_precision(double input_precision_t);
int get_linear_solver(int *linear_solver_t);
int set_linear_solver(int linear_solver_t);
int get_check_for_dynamical_stability(bool *value);
int set_check_for_dynamical_stability(int value);
int get_check_for_inner_collision(bool *value);
int set_check_for_inner_collision(int value);
int get_check_for_outer_collision(bool *value);
int set_check_for_outer_collision(int value);
int get_check_for_inner_RLOF(bool *value);
int set_check_for_inner_RLOF(int value);
int get_check_for_outer_RLOF(bool *value);
int set_check_for_outer_RLOF(int value);
int get_include_quadrupole_terms(bool *value);
int set_include_quadrupole_terms(int value);
int get_include_octupole_terms(bool *value);
int set_include_octupole_terms(int value);
int get_include_1PN_inner_terms(bool *value);
int set_include_1PN_inner_terms(int value);
int get_include_1PN_outer_terms(bool *value);
int set_include_1PN_outer_terms(int value);
int get_include_1PN_inner_outer_terms(bool *value);
int set_include_1PN_inner_outer_terms(int value);
int get_include_25PN_inner_terms(bool *value);
int set_include_25PN_inner_terms(int value);
int get_include_25PN_outer_terms(bool *value);
int set_include_25PN_outer_terms(int value);
int get_include_inner_tidal_terms(bool *value);
int set_include_inner_tidal_terms(int value);
int get_include_outer_tidal_terms(bool *value);
int set_include_outer_tidal_terms(int value);
int get_ignore_tertiary(bool *value);
int set_ignore_tertiary(int value);
int get_include_inner_wind_terms(bool *value);
int set_include_inner_wind_terms(int value);
int get_include_outer_wind_terms(bool *value);
int set_include_outer_wind_terms(int value);
int get_include_magnetic_braking_terms(bool *value);
int set_include_magnetic_braking_terms(int value);
int get_include_spin_radius_mass_coupling_terms_star1(bool *value);
int set_include_spin_radius_mass_coupling_terms_star1(int value);
int get_include_spin_radius_mass_coupling_terms_star2(bool *value);
int set_include_spin_radius_mass_coupling_terms_star2(int value);
int get_include_spin_radius_mass_coupling_terms_star3(bool *value);
int set_include_spin_radius_mass_coupling_terms_star3(int value);
int get_include_inner_RLOF_terms(bool *value);
int set_include_inner_RLOF_terms(int value);
int get_include_outer_RLOF_terms(bool *value);
int set_include_outer_RLOF_terms(int value);
int get_include_linear_mass_change(bool *value);
int set_include_linear_mass_change(int value);
int get_include_linear_radius_change(bool *value);
int set_include_linear_radius_change(int value);
int get_relative_tolerance(double *value);
int set_relative_tolerance(double value);
int get_check_for_dynamical_stability_at_initialisation(bool *value);
int set_check_for_dynamical_stability_at_initialisation(int value);
int get_check_for_semisecular_regime_at_initialisation(bool *value);
int set_check_for_semisecular_regime_at_initialisation(int value);
int get_verbose(bool *value);
int set_verbose(int value);



/* addition Jan 2017 */
int get_check_for_semisecular_regime(bool *value);
int set_check_for_semisecular_regime(bool value);
int get_check_for_semisecular_regime_parameter(double *value);
int set_check_for_semisecular_regime_parameter(double value);
