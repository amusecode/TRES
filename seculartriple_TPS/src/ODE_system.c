/*	Worker code for SecularTriple, a secular triple gravitational dynamics code taking into account Newtonian, 1PN and 2.5PN terms	*/
/*  Also included: effects of tidal friction, wind mass loss & mass transfer */
/*	The relevant ODEs are solved consistently for each user supplied timestep using CVODE (Cohen & Hindmarsh 1996)	*/

#include "main_code.h"

int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data_f)
{
    
	UserData data;
	data = (UserData) data_f;
	
    /**********************
     * extract parameters *
     **********************/
    double global_time_step = data->global_time_step; // the global time-step
    
    int stellar_type1 = data->stellar_type1;
    int stellar_type2 = data->stellar_type2;
    int stellar_type3 = data->stellar_type3;
	double m1 = data->m1; // masses -- at the END of the global time-step
	double m2 = data->m2;				
	double m3 = data->m3;		
	double m1_convective_envelope = data->m1_convective_envelope;
	double m2_convective_envelope = data->m2_convective_envelope;
	double m3_convective_envelope = data->m3_convective_envelope;
    double R1 = data->R1; // radii -- at the END of the global time-step
    double R2 = data->R2;
    double R3 = data->R3;
    double R1_convective_envelope = data->R1_convective_envelope; // convective envelope radii
    double R2_convective_envelope = data->R2_convective_envelope;
    double R3_convective_envelope = data->R3_convective_envelope;

    double AMC_star1 = data->AMC_star1; // Apsidal Motion Constant
    double AMC_star2 = data->AMC_star2; // Apsidal Motion Constant
    double AMC_star3 = data->AMC_star3; // Apsidal Motion Constant  
    
//    printf("apsidal motion constant %g %g %g \n",AMC_star1, AMC_star2, AMC_star3);
    
    double luminosity_star1 = data->luminosity_star1;
    double luminosity_star2 = data->luminosity_star2;
    double luminosity_star3 = data->luminosity_star3;
    double gyration_radius_star1 = data->gyration_radius_star1; // gyration radius (NOT squared)     
    double gyration_radius_star2 = data->gyration_radius_star2; // gyration radius (NOT squared)     
    double gyration_radius_star3 = data->gyration_radius_star3; // gyration radius (NOT squared)             


    /* !!!!!!!!!!!!!!!!!!!!!!!! */
    /* HARDCODED gyration radii that overwrite the above user-specified gyration radii */
    //Gabriele
//   gyration_radius_star1 = set_crude_gyration_radii_based_on_stellar_structure(stellar_type1,m1);
//   gyration_radius_star2 = set_crude_gyration_radii_based_on_stellar_structure(stellar_type2,m2);
//   gyration_radius_star3 = set_crude_gyration_radii_based_on_stellar_structure(stellar_type3,m3);

    
    double moment_of_inertia_star1 = data->moment_of_inertia_star1;
    double moment_of_inertia_star2 = data->moment_of_inertia_star2;
    double moment_of_inertia_star3 = data->moment_of_inertia_star3;
    double moment_of_inertia_dot_star1 = data->moment_of_inertia_dot_star1;
    double moment_of_inertia_dot_star2 = data->moment_of_inertia_dot_star2;
    double moment_of_inertia_dot_star3 = data->moment_of_inertia_dot_star3;

    /* NOTE: k_div_T_tides_stari is not assumed to be constant during the integration */
//    double k_div_T_tides_star1 = data->k_div_T_tides_star1; // AMC divided by tidal dissipation time-scale
//    double k_div_T_tides_star2 = data->k_div_T_tides_star2;    
//    double k_div_T_tides_star3 = data->k_div_T_tides_star3;

    bool include_quadrupole_terms = data->include_quadrupole_terms;
    bool include_octupole_terms = data->include_octupole_terms;
    bool include_1PN_inner_terms = data->include_1PN_inner_terms;
    bool include_1PN_outer_terms = data->include_1PN_outer_terms;
    bool include_1PN_inner_outer_terms = data->include_1PN_inner_outer_terms;
    bool include_25PN_inner_terms = data->include_25PN_inner_terms;
    bool include_25PN_outer_terms = data->include_25PN_outer_terms;
    bool include_inner_tidal_terms = data->include_inner_tidal_terms;
    bool include_outer_tidal_terms = data->include_outer_tidal_terms;
    bool ignore_tertiary = data->ignore_tertiary;
    
    bool include_linear_mass_change = data->include_linear_mass_change;
    bool include_linear_radius_change = data->include_linear_radius_change;

    bool include_inner_wind_terms = data->include_inner_wind_terms;
    bool include_outer_wind_terms = data->include_outer_wind_terms;
    bool include_magnetic_braking_terms = data->include_magnetic_braking_terms;
    bool include_spin_radius_mass_coupling_terms_star1 = data->include_spin_radius_mass_coupling_terms_star1;
    bool include_spin_radius_mass_coupling_terms_star2 = data->include_spin_radius_mass_coupling_terms_star2;
    bool include_spin_radius_mass_coupling_terms_star3 = data->include_spin_radius_mass_coupling_terms_star3;
    bool include_inner_RLOF_terms = data->include_inner_RLOF_terms;
    bool include_outer_RLOF_terms = data->include_outer_RLOF_terms;
    bool star1_is_donor = data->star1_is_donor;
    bool star2_is_donor = data->star2_is_donor;
    bool star3_is_donor = data->star3_is_donor;
    double wind_mass_loss_rate_star1 = data->wind_mass_loss_rate_star1;
    double wind_mass_loss_rate_star2 = data->wind_mass_loss_rate_star2;
    double wind_mass_loss_rate_star3 = data->wind_mass_loss_rate_star3;
    double R1_dot = data->time_derivative_of_radius_star1;
    double R2_dot = data->time_derivative_of_radius_star2;
    double R3_dot = data->time_derivative_of_radius_star3;
    double inner_mass_transfer_rate = data->inner_mass_transfer_rate;
    double outer_mass_transfer_rate = data->outer_mass_transfer_rate;
    double inner_accretion_efficiency_wind_star1_to_star2 = data->inner_accretion_efficiency_wind_child1_to_child2;
    double inner_accretion_efficiency_wind_star2_to_star1 = data->inner_accretion_efficiency_wind_child2_to_child1;
    double outer_accretion_efficiency_wind_star1_to_star3 = data->outer_accretion_efficiency_wind_child2_to_child1; /* temporary */
    double outer_accretion_efficiency_wind_star2_to_star3 = data->outer_accretion_efficiency_wind_child2_to_child1; /* temporary */
    double outer_accretion_efficiency_wind_star3_to_inner_binary = data->outer_accretion_efficiency_wind_child1_to_child2;
    
    double inner_accretion_efficiency_mass_transfer = data->inner_accretion_efficiency_mass_transfer;
    double outer_accretion_efficiency_mass_transfer = data->outer_accretion_efficiency_mass_transfer;
    double inner_specific_AM_loss_mass_transfer = data->inner_specific_AM_loss_mass_transfer;
    double outer_specific_AM_loss_mass_transfer = data->outer_specific_AM_loss_mass_transfer;
    double inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2 = data->inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2;
    double inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1 = data->inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1;    

    double threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = data->threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero;
    double threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero = data->threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero;


    /******************
     * compute mdots *
     *****************/
    
     /* assumptions:
     * -- in inner binary, child1 is star1 & child2 is star2
     * -- wind mass loss rates are always 0 or negative
     * -- mass transfer rates are always 0 or negative
     */

    /* wind mass loss */
    double m1_dot = wind_mass_loss_rate_star1;
    double m2_dot = wind_mass_loss_rate_star2;
    double m3_dot = wind_mass_loss_rate_star3;

    m1_dot += -inner_accretion_efficiency_wind_star2_to_star1*wind_mass_loss_rate_star2; // wind from star2 to star1; minus sign because wind mass loss rates are always negative
    m2_dot += -inner_accretion_efficiency_wind_star1_to_star2*wind_mass_loss_rate_star1; // wind from star1 to star2; minus sign because wind mass loss rates are always negative    
    
    /* mass transfer */
    double effective_inner_mass_transfer_rate = inner_accretion_efficiency_mass_transfer*inner_mass_transfer_rate; // always 0 or negative
    if (star1_is_donor == TRUE)  // m1 transfers mass; m2 gains a fraction of this
    {
        m1_dot += inner_mass_transfer_rate;
        m2_dot += -effective_inner_mass_transfer_rate;
    }
    else if (star2_is_donor == TRUE) // m2 transfers mass; m1 gains a fraction of this
    {
        m1_dot += -effective_inner_mass_transfer_rate;
        m2_dot += inner_mass_transfer_rate;
    }

    /***************************************************************
     * compute time-dependent masses, radii and moments of inertia *
     ***************************************************************/
     
    /* masses and radii extracted above are the FINAL values
     * let q denote mi, Ri or Ii, then
     * q(t) = q_begin + q_dot*t, where t the time relative to the time for which q = q_begin
     * q_end = q_begin + q_dot*dt, where dt is the global time-step
     * hence q(t) = q_end + q_dot*(t-dt) */


    double t_dif = t-global_time_step; // i.e., t-dt
    if (include_linear_mass_change == TRUE)
    {
        m1 += m1_dot*t_dif; // note that before this statement, m1 corresponded to q_end
        m2 += m2_dot*t_dif;
        m3 += m3_dot*t_dif;
    }
    if (include_linear_radius_change == TRUE)
    {
        R1 += R1_dot*t_dif;
        R2 += R2_dot*t_dif;
        R3 += R3_dot*t_dif;
    }

    moment_of_inertia_star1 += moment_of_inertia_dot_star1*t_dif;
    moment_of_inertia_star2 += moment_of_inertia_dot_star2*t_dif;
    moment_of_inertia_star3 += moment_of_inertia_dot_star3*t_dif;

	/*	the ODE variables	
     *  g: argument of pericentre
     *  h: longitude of the ascending nodes
     */
     
	double x = Ith(yev,1); // log_10(1-e_in)
	double y = Ith(yev,2); // log_10(1-e_out)
	double g_in = Ith(yev,3);
	double g_out = Ith(yev,4);
    double h_in = Ith(yev,5);
    double h_out = Ith(yev,6);
	double a_in = Ith(yev,7);
	double a_out = Ith(yev,8);
	double cositot = Ith(yev,9);
	double spin_angular_frequency1 = Ith(yev,10);
	double spin_angular_frequency2 = Ith(yev,11);
	double spin_angular_frequency3 = Ith(yev,12);

	double e_in = 1.0 - pow(10.0,x);
	double e_out = 1.0 - pow(10.0,y);
    
    /*  track the positive changes in eccentricity to obtain
     *  the largest Kozai-Lidov amplitude
     */
    
    extern double e_in_prev;
    extern double tracker;
    extern double delta_e_in;
    extern double t_prev;

    double diff = e_in - e_in_prev;
    
    if (t >= t_prev) {
        if (diff >= -1.0e-5) {
            tracker += diff;
            if (tracker > delta_e_in) {
                delta_e_in = tracker;
            }
        }
        else {
            if (tracker > delta_e_in) {
                delta_e_in = tracker;
            }
            tracker = 0;
        }
    }
    else{
        tracker += diff;
    }
    t_prev = t;
    e_in_prev = e_in;
 
    /* in the case that the tertiary is not to be taken into account (ignore_tertiary == TRUE),
     * several quantities should still be set to some arbitrary, nonzero value in order to avoid nans at various instances
     * this does not affect the values of the quantities outside of this function,
     * nor the dots of these quantities computed below */
    if (ignore_tertiary == TRUE)
    {
        m3 = 1.0;
        a_out = 1.0e10;
        e_out = 0.1;
    }
    
    /* premature return values
     * currently not used */
    if (e_in < 0.0)
    {
        printf("e_in<=0 %g \n",e_in);
//        return 1;
    }
    if (e_out < 0.0)
    {
        printf("e_out<=0\n");
//        return 1;
    }
    if (data->stop_after_error_bool == TRUE)
    {
        printf("error occured -- stopping integration\n");
        return -1;
    }

    /********************************
     * preamble: derived quantities *
     * *****************************/
     
	/*	mass quantities	*/
	double m1_div_m2 = m1/m2;
	double m2_div_m1 = 1.0/m1_div_m2;    
    double m1_plus_m2 = m1+m2;
    double m1_plus_m2_plus_m3 = m1_plus_m2+m3;
    double m1_times_m2 = m1*m2;
    double m1_times_m2_times_m3 = m1_times_m2*m3;
    double m1_plus_m2_div_m3 = m1_plus_m2/m3;
    double m1_times_m1 = m1*m1;
    double m2_times_m2 = m2*m2;
        
	/*	eccentricity quantities	*/
	double e_in_p2 = e_in*e_in;
	double l_in_p2 = 1.0 - e_in_p2;
    double l_in = sqrt(l_in_p2);
    double l_in_p3 = l_in_p2*l_in;
    double l_in_p4 = l_in_p3*l_in;
    double l_in_p5 = l_in_p4*l_in;
    double l_in_p6 = l_in_p5*l_in;
    double l_in_p7 = l_in_p6*l_in;
    
	double e_out_p2 = e_out*e_out;
	double l_out_p2 = 1.0 - e_out_p2;
    double l_out = sqrt(l_out_p2);
    double l_out_p3 = l_out_p2*l_out;
    double l_out_p4 = l_out_p3*l_out;
    double l_out_p5 = l_out_p4*l_out;    
    double l_out_p6 = l_out_p5*l_out;
    double l_out_p7 = l_out_p6*l_out;

	/*	triple secular gravitational dynamics quantities */
    /* 2000ApJ...535..385F */
	double L_in = m1_times_m2*sqrt(CONST_G*a_in/(m1_plus_m2));
	double L_out = (m1_plus_m2)*m3*sqrt(CONST_G*a_out/(m1_plus_m2_plus_m3));
	double G_in = L_in*sqrt(l_in_p2);
	double G_out = L_out*sqrt(l_out_p2);
    double G_tot = sqrt( G_in*G_in + G_out*G_out + 2.0*G_in*G_out*cositot );

	double a_in_div_a_out = a_in/a_out;
    double C2,C3;
    if (include_quadrupole_terms == FALSE)
    {
        C2 = 0.0;
    }
    else
    {
        C2 = CONST_G*c_1div16*(m1_times_m2_times_m3/m1_plus_m2)*pow(l_out,-3.0)*a_in_div_a_out*a_in_div_a_out/a_out;
    }
    if (include_octupole_terms == FALSE)
    {
        C3 = 0.0;
    }
    else
    {
        C3 = -CONST_G*c_15div16*c_1div4*(m1_times_m2_times_m3/(m1_plus_m2*m1_plus_m2))*(m1-m2)*pow(l_out,-5.0)*a_in_div_a_out*a_in_div_a_out*a_in_div_a_out/a_out;
    }

	if (cositot > 1.0)
	{
		cositot = 2.0 - cositot;
	}

	if (cositot < -1.0)
	{
		cositot = -2.0 - cositot;
	}

	double cositot_p2 = cositot*cositot;
	double sinitot = sqrt(1.0 - cositot_p2); // NOTE: 0 < itot < PI, so sinitot > 0 always
	double sinitot_p2 = sinitot*sinitot;

	double sin_g_in = sin(g_in);
	double sin_2g_in = sin(2.0*g_in);
	double sin_g_out = sin(g_out);
	double cos_g_in = cos(g_in);
	double cos_2g_in = cos(2.0*g_in);
	double cos_g_out = cos(g_out);
	
	/*	required for octupole-order terms	*/
	double B = 2.0 + 5.0*e_in_p2 - 7.0*e_in_p2*cos_2g_in;
	double A = 4.0 + 3.0*e_in_p2 - c_5div2*B*sinitot_p2;
	double cosphi = -cos_g_in*cos_g_out - cositot*sin_g_in*sin_g_out;

    /* PN quantities */
    double f_25PN_e_in = 0.0, f_25PN_a_in = 0.0;
    if (include_25PN_inner_terms == TRUE)
    {
        f_25PN_e_in = f_25PN_e(e_in_p2);
        f_25PN_a_in = f_25PN_a(e_in_p2);
    }
    double f_25PN_e_out = 0.0, f_25PN_a_out = 0.0;
    if (include_25PN_outer_terms == TRUE)
    {
        f_25PN_e_out = f_25PN_e(e_out_p2);
        f_25PN_a_out = f_25PN_a(e_out_p2);
    }
    double f_1PN_in_out_m1_m2 = 0.0, f_1PN_in_out_L_in_tilde = 0.0, f_1PN_in_out_L_out_tilde =  0.0, f_1PN_in_out_LL = 0.0, f_1PN_in_out_e_in = 0.0, f_1PN_in_out_i = 0.0;
    if (include_1PN_inner_outer_terms == TRUE)
    {
        f_1PN_in_out_m1_m2 = m1_times_m1 + m1_times_m2 + m2_times_m2;
        f_1PN_in_out_L_in_tilde = L_in/(m1_times_m2/m1_plus_m2);
        f_1PN_in_out_L_out_tilde = L_out/(m1_plus_m2*m3/m1_plus_m2_plus_m3);
        f_1PN_in_out_LL = f_1PN_in_out_L_in_tilde*f_1PN_in_out_L_out_tilde*m1_plus_m2*(4.0*m1_plus_m2 + 3.0*m3);
        f_1PN_in_out_e_in = (2.0 - 5.0*e_in_p2)*m1_times_m1 + 3.0*(-2.0 + e_in_p2)*m1_times_m1 + (2.0 - 5.0*e_in_p2)*m2_times_m2;
        f_1PN_in_out_i = 3.0*a_in*m1_plus_m2_plus_m3*cositot*(f_1PN_in_out_e_in + 3.0*e_in_p2*f_1PN_in_out_m1_m2*cos_2g_in);
    }
    	
    /* tides quantities */
    double k_div_T_tides_star1 = compute_k_div_T_tides(stellar_type1,m1,m1_convective_envelope,m2,a_in,R1,R1_convective_envelope,luminosity_star1,spin_angular_frequency1,gyration_radius_star1, AMC_star1); // AMC divided by tidal dissipation time-scale
    double k_div_T_tides_star2 = compute_k_div_T_tides(stellar_type2,m2,m2_convective_envelope,m1,a_in,R2,R2_convective_envelope,luminosity_star2,spin_angular_frequency2,gyration_radius_star2, AMC_star2); // AMC divided by tidal dissipation time-scale
    double k_div_T_tides_star3 = compute_k_div_T_tides(stellar_type3,m3,m3_convective_envelope,m1+m2,a_out,R3,R3_convective_envelope,luminosity_star3,spin_angular_frequency3,gyration_radius_star3, AMC_star3); // AMC divided by tidal dissipation time-scale

//    printf("k_div_T %g %g %g \n",k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3);

	double R1_div_a_in = R1/a_in;
	double R1_div_a_in_p2 = R1_div_a_in*R1_div_a_in;
	double R1_div_a_in_p5 = R1_div_a_in_p2*R1_div_a_in_p2*R1_div_a_in;
	double R1_div_a_in_p6 = R1_div_a_in*R1_div_a_in_p5;
	double R2_div_a_in = R2/a_in;
	double R2_div_a_in_p2 = R2_div_a_in*R2_div_a_in;
	double R2_div_a_in_p5 = R2_div_a_in_p2*R2_div_a_in_p2*R2_div_a_in;
	double R2_div_a_in_p6 = R2_div_a_in*R2_div_a_in_p5;
	double R3_div_a_out = R3/a_out;
	double R3_div_a_out_p2 = R3_div_a_out*R3_div_a_out;
	double R3_div_a_out_p5 = R3_div_a_out_p2*R3_div_a_out_p2*R3_div_a_out;
	double R3_div_a_out_p6 = R3_div_a_out*R3_div_a_out_p5;
	
	double m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1 = m2_div_m1*R1_div_a_in_p6*k_div_T_tides_star1;
	double m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2 = m1_div_m2*R2_div_a_in_p6*k_div_T_tides_star2;
	double m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3 = m1_plus_m2_div_m3*R3_div_a_out_p6*k_div_T_tides_star3;
    
    double a_in_p2 = a_in*a_in;
    double a_in_p3 = a_in_p2*a_in;
    double a_in_p4 = a_in_p3*a_in;
	double n_in = sqrt(CONST_G*m1_plus_m2/a_in_p3); // mean orbital angular speed
    double spin_angular_frequency1_div_n_in = spin_angular_frequency1/n_in;
    double spin_angular_frequency2_div_n_in = spin_angular_frequency2/n_in;

    double a_out_p2 = a_out*a_out;
    double a_out_p3 = a_out_p2*a_out;
    double a_out_p4 = a_out_p3*a_out;
	double n_out = sqrt(CONST_G*m1_plus_m2_plus_m3/a_out_p3); // mean orbital angular speed
    double spin_angular_frequency3_div_n_out = spin_angular_frequency3/n_out;

    double f_tides1_in = f_tides1(e_in_p2);
    double f_tides2_in = f_tides2(e_in_p2);
    double f_tides3_in = f_tides3(e_in_p2);
    double f_tides4_in = f_tides4(e_in_p2);
    double f_tides5_in = f_tides5(e_in_p2);

    double f_tides1_out = f_tides1(e_out_p2);
    double f_tides2_out = f_tides2(e_out_p2);
    double f_tides3_out = f_tides3(e_out_p2);
    double f_tides4_out = f_tides4(e_out_p2);
    double f_tides5_out = f_tides5(e_out_p2);

    

    /************************************************
     * the calculations of the ODE right-hand-sides *
     * **********************************************/
     	

    /*******************************
     * e_in_dot                    *
     * *****************************/
    
	double e_in_dot_newtonian = 0.0;
    double e_in_dot_GR_1PN_in_out = 0.0;
    double e_in_dot_GR_25PN_in = 0.0;
    double e_in_dot_tides = 0.0;    

    /* Newtonian point particle -- up and including octupole order */    
    if (ignore_tertiary == FALSE)
    {
        e_in_dot_newtonian = C2*(l_in_p2/G_in)*(30.0*e_in*sinitot_p2*sin_2g_in) \
            + C3*e_out*(l_in_p2/G_in)*(35.0*cosphi*sinitot_p2*e_in_p2*sin_2g_in \
                - 10.0*cositot*sinitot_p2*cos_g_in*sin_g_out*l_in_p2 \
                - A*(sin_g_in*cos_g_out - cositot*cos_g_in*sin_g_out));
                
        printf("e_in_dot_new %g %g %g %g \n",C2*(l_in_p2/G_in)*(30.0*e_in*sinitot_p2*sin_2g_in), C3*e_out*(l_in_p2/G_in)*(35.0*cosphi*sinitot_p2*e_in_p2*sin_2g_in), - 10.0*cositot*sinitot_p2*cos_g_in*sin_g_out*l_in_p2*C3*e_out*(l_in_p2/G_in), - A*(sin_g_in*cos_g_out - cositot*cos_g_in*sin_g_out)*C3*e_out*(l_in_p2/G_in));
                
                
    }
    
    /* post-Newtonian point particle */
    if ((ignore_tertiary == FALSE) && (include_1PN_inner_outer_terms == TRUE))
    {
        e_in_dot_GR_1PN_in_out = -c_9div16*CONST_G*CONST_G*a_in*e_in*l_in*m1_times_m2*f_1PN_in_out_m1_m2*m3*sinitot_p2*sin_2g_in/( CONST_C_LIGHT_P2*a_out_p3*l_out_p3*L_in*m1_plus_m2*m1_plus_m2);
    }
    if (include_25PN_inner_terms == TRUE)
    {
        e_in_dot_GR_25PN_in = -c_304div15*CONST_G_P3*m1_times_m2*m1_plus_m2*(e_in/(CONST_C_LIGHT_P5*a_in_p4*l_in_p5))*f_25PN_e_in;
    }
    
    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        double e_in_dot_tides_star1 = -27.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency1_div_n_in);
        double e_in_dot_tides_star2 = -27.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency2_div_n_in);            

        if (e_in <= threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero)
        {
            e_in_dot_tides_star1 = e_in_dot_tides_star2 = 0.0;
//            printf("e_in_dot_tides %g\n",e_in_dot_tides_star1);
            
        }
        e_in_dot_tides = e_in_dot_tides_star1 + e_in_dot_tides_star2;
//        printf("e_in_dot_tides %g %g %g \n",e_in_dot_tides, e_in_dot_tides_star1, e_in_dot_tides_star2);
    }
    
    /* mass transfer */
    /* wind mass loss */

    /* combined */

    double e_in_dot = 0.0;

    if (e_in <= threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero)
	{
//        printf("setting e_in_dot zero %g %g %g %g %g \n",e_in,e_in_dot_newtonian,e_in_dot_GR_1PN_in_out,e_in_dot_GR_25PN_in,e_in_dot_tides);
	    e_in_dot = 0.0;
	}
    else
	{
	    e_in_dot = e_in_dot_newtonian + e_in_dot_GR_1PN_in_out + e_in_dot_GR_25PN_in + e_in_dot_tides;
    }
	Ith(ydot,1) = -1.0*pow(10.0,-x)*e_in_dot/log(10.0);

printf("e_in_dot %g %g %g %g %g %g\n",e_in,e_in_dot_newtonian,e_in_dot_GR_1PN_in_out,e_in_dot_GR_25PN_in,e_in_dot_tides, threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero);


    /*******************************
     * e_out_dot                       *
     * *****************************/

    double e_out_dot_newtonian = 0.0;
    double e_out_dot_GR_25PN_out = 0.0;
    double e_out_dot_tides = 0.0;
    
    if (ignore_tertiary == FALSE)
    {
        /* Newtonian point particle -- up and including octupole order */
        e_out_dot_newtonian = -C3*e_in*(l_out_p2/G_out)*(10.0*cositot*sinitot_p2*l_in_p2*sin_g_in*cos_g_out \
            + A*(cos_g_in*sin_g_out - cositot*sin_g_in*cos_g_out));
            
        /* post-Newtonian point particle */        
        if (include_25PN_outer_terms == TRUE)
        {
            e_out_dot_GR_25PN_out = -c_304div15*CONST_G_P3*m1_plus_m2*m3*m1_plus_m2_plus_m3*(e_out/(CONST_C_LIGHT_P5*a_out_p4*l_out_p5))*f_25PN_e_out;
        }
    
        /* tides */
        if (include_outer_tidal_terms == TRUE) // ad hoc/approximate: treats inner binary as point particle
        {
            double e_out_dot_tides_star3 = -27.0*(1.0+m1_plus_m2_div_m3)*m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3*R3_div_a_out_p2*e_out*pow(l_out,-13.0)*(f_tides3_out \
                - c_11div18*l_out_p3*f_tides4_out*spin_angular_frequency3_div_n_out);
            e_out_dot_tides = e_out_dot_tides_star3;
        }

        /* mass transfer */
        /* wind mass loss */
    }
    
    /* combined */
	double e_out_dot = e_out_dot_newtonian + e_out_dot_GR_25PN_out + e_out_dot_tides;
	Ith(ydot,2) = -1.0*pow(10.0,-y)*e_out_dot/log(10.0);
printf("e_out_dot %g %g %g %g \n",e_out,e_out_dot_newtonian,e_out_dot_GR_25PN_out,e_out_dot_tides);



    /*******************************
     * g_in_dot                    *
     * *****************************/

	double g_in_dot_newtonian = 0.0;
    double g_in_dot_GR_1PN_in = 0.0;
    double g_in_dot_GR_1PN_in_out = 0.0;
    double g_in_dot_tides = 0.0;    
     
    /* Newtonian point particle -- up and including octupole order */
    if (ignore_tertiary == FALSE)
    {
        g_in_dot_newtonian = 6.0*C2*((1.0/G_in)*(4.0*cositot_p2 + (5.0*cos_2g_in - 1.0)*(l_in_p2 - cositot_p2)) \
                + (cositot/G_out)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in))) \
            - C3*e_out*(e_in*((1.0/G_out) + (cositot/G_in))*(sin_g_in*sin_g_out*(10.0*(3.0*cositot_p2 \
                - 1.0)*(1.0 - e_in_p2) + A) - 5.0*B*cositot*cosphi) \
                - (l_in_p2/(e_in*G_in))*(sin_g_in*sin_g_out*10.0*cositot*sinitot_p2*(1.0 - 3.0*e_in_p2) \
                + cosphi*(3.0*A - 10.0*cositot_p2 + 2.0)));


//        printf("g_in_dot newtonian %g %g  %g \n", g_in_dot_newtonian,6.0*C2*((1.0/G_in)*(4.0*cositot_p2 + (5.0*cos_2g_in - 1.0)*(l_in_p2 - cositot_p2)) 
//                + (cositot/G_out)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in))) , C3*e_out*(e_in*((1.0/G_out) + (cositot/G_in))*(sin_g_in*sin_g_out*(10.0*(3.0*cositot_p2  - 1.0)*(1.0 - e_in_p2) + A) - 5.0*B*cositot*cosphi) - (l_in_p2/(e_in*G_in))*(sin_g_in*sin_g_out*10.0*cositot*sinitot_p2*(1.0 - 3.0*e_in_p2) + cosphi*(3.0*A - 10.0*cositot_p2 + 2.0))));


    }
    
    /* post-Newtonian point particle */  
    if (include_1PN_inner_terms == TRUE)
    {
        g_in_dot_GR_1PN_in = (3.0/(CONST_C_LIGHT_P2*a_in*l_in_p2))*pow(CONST_G*m1_plus_m2/a_in,3.0/2.0);
    }
    if ((ignore_tertiary == FALSE) && (include_1PN_inner_outer_terms == TRUE))
    {
        g_in_dot_GR_1PN_in_out = CONST_G*CONST_G*m1_times_m2_times_m3/(16.0*a_out_p3*CONST_C_LIGHT_P2*l_out_p3*m1_plus_m2*m1_plus_m2*m1_plus_m2_plus_m3)*( (a_in/G_in)*m1_plus_m2_plus_m3*( \
            l_in_p2*(5.0*m1_times_m2 - 3.0*m1_times_m2 + 5.0*m2_times_m2) - 9.0*f_1PN_in_out_m1_m2*( l_in_p2*cos_2g_in + 2.0*cositot_p2*sin_g_in*sin_g_in ) ) \
            + (1.0/G_out)*(-8.0*f_1PN_in_out_LL + f_1PN_in_out_i) );
    }
	
    /* tides/rotation */
    if (include_inner_tidal_terms == TRUE)
    {    
        double g_in_dot_tides_star1 = R1_div_a_in_p5*AMC_star1*(n_in/(l_in_p4))*(15.0*m2_div_m1*f_tides4_in/l_in_p6 \
		+ (1.0+m2_div_m1)*spin_angular_frequency1_div_n_in*spin_angular_frequency1_div_n_in);
        double g_in_dot_tides_star2 = R2_div_a_in_p5*AMC_star2*(n_in/(l_in_p4))*(15.0*m1_div_m2*f_tides4_in/l_in_p6 \
		+ (1.0+m1_div_m2)*spin_angular_frequency2_div_n_in*spin_angular_frequency2_div_n_in);
        g_in_dot_tides = g_in_dot_tides_star1 + g_in_dot_tides_star2;


//        printf("g_in_dot detail %g %g  %g %g  %g %g\n", g_in_dot_tides_star1, g_in_dot_tides_star2, 15.0*m2_div_m1*f_tides4_in/l_in_p6, (1.0+m2_div_m1)*spin_angular_frequency1_div_n_in*spin_angular_frequency1_div_n_in, 15.0*m1_div_m2*f_tides4_in/l_in_p6, (1.0+m1_div_m2)*spin_angular_frequency2_div_n_in*spin_angular_frequency2_div_n_in);



    }
    
    /* wind mass loss */
    /* mass transfer */    
    
    /* combined */  
    double g_in_dot = g_in_dot_newtonian + g_in_dot_GR_1PN_in + g_in_dot_GR_1PN_in_out + g_in_dot_tides;
	Ith(ydot,3) = g_in_dot;


printf("g_in_dot %g %g %g %g %g %g\n",g_in,g_in_dot, g_in_dot_newtonian, g_in_dot_GR_1PN_in, g_in_dot_GR_1PN_in_out, g_in_dot_tides);



    /*******************************
     * g_out_dot                   *
     * *****************************/
     
    double g_out_dot_newtonian = 0.0;
    double g_out_dot_GR_1PN_out = 0.0;
    double g_out_dot_GR_1PN_in_out = 0.0;
    double g_out_dot_tides = 0.0;
    
    if (ignore_tertiary == FALSE)
    {
        /* Newtonian point particle -- up and including octupole order */
        g_out_dot_newtonian = 3.0*C2*((2.0*cositot/G_in)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)) \
                + (1.0/G_out)*(4.0 + 6.0*e_in_p2 + (5.0*cositot_p2 - 3.0)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)))) \
            + C3*e_in*(sin_g_in*sin_g_out*(((4.0*e_out_p2 + 1.0)/(e_out*G_out))*10.0*cositot*sinitot_p2*l_in_p2 \
                - e_out*((1.0/G_in) + (cositot/G_out))*(A + 10.0*(3.0*cositot_p2 - 1.0)*l_in_p2)) \
                + cosphi*(5.0*B*cositot*e_out*((1.0/G_in) + (cositot/G_out)) + ((4.0*e_out_p2 + 1.0)/(e_out*G_out))*A));
    
        /* post-Newtonian point particle */  
        if (include_1PN_outer_terms == TRUE)
        {    
            g_out_dot_GR_1PN_out = (3.0/(CONST_C_LIGHT_P2*a_out*l_out_p2))*pow(CONST_G*m1_plus_m2_plus_m3/a_out,3.0/2.0);
        }
        if (include_1PN_inner_outer_terms == TRUE)
        {
            g_out_dot_GR_1PN_in_out = CONST_G*CONST_G*m1_times_m2_times_m3/(16.0*a_out_p3*CONST_C_LIGHT_P2*l_out_p3*m1_plus_m2*m1_plus_m2*m1_plus_m2_plus_m3)*( (8.0*f_1PN_in_out_LL - f_1PN_in_out_i) \
                - (1.0/(2.0*G_out))*( 2.0*cositot*(-8.0*f_1PN_in_out_LL + f_1PN_in_out_i) - 16.0*m1_plus_m2*m1_plus_m2*cositot*f_1PN_in_out_L_in_tilde*f_1PN_in_out_L_out_tilde*( \
                    16.0*m1_plus_m2*f_1PN_in_out_L_in_tilde*f_1PN_in_out_L_out_tilde*(7.0*m1_plus_m2 + 6.0*m3)*cositot \
                        + c_3div2*a_in*m1_plus_m2_plus_m3*(-f_1PN_in_out_e_in*(1.0 + 3.0*(2.0*cositot_p2 - 1.0)) + 18.0*e_in_p2*f_1PN_in_out_m1_m2*cos_2g_in*sinitot_p2) ) ) );
        }
    
        /* tides */
        if (include_outer_tidal_terms == TRUE)
        {    
            double g_out_dot_tides_star3 = R3_div_a_out_p5*AMC_star3*(n_out/(l_out_p4))*(15.0*m1_plus_m2_div_m3*f_tides4_out/l_out_p6 \
            + (1.0+m1_plus_m2_div_m3)*spin_angular_frequency3_div_n_out*spin_angular_frequency3_div_n_out);
            g_out_dot_tides = g_out_dot_tides_star3;
        }
        
        /* wind mass loss */
        /* mass transfer */    
    }
    
    /* combined */
    double g_out_dot = g_out_dot_newtonian + g_out_dot_GR_1PN_out + g_out_dot_GR_1PN_in_out + g_out_dot_tides;
	Ith(ydot,4) = g_out_dot;
	
//    printf("g_out_dot %g %g %g %g %g %g \n",g_out, g_out_dot, g_out_dot_newtonian, g_out_dot_GR_1PN_out, g_out_dot_GR_1PN_in_out, g_out_dot_tides);
    
    
    /*******************************
     * h_in_dot                    *
     * *****************************/
     
    double h_in_dot_newtonian = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */
    /* 2013MNRAS.431.2155N Eq. B8 */
    if (ignore_tertiary == FALSE)
    {
        h_in_dot_newtonian = -3.0*C2*(G_tot/(G_in*G_out*sinitot))*(2.0 + 3.0*e_in_p2 - 5.0*e_in_p2*cos_2g_in)*2.0*sinitot*cositot \
            - C3*e_in*e_out*(G_tot/(G_in*G_out))*(5.0*B*cositot*cosphi - A*sin_g_in*sin_g_out + 10.0*(1.0 - 3.0*cositot_p2)*l_in_p2*sin_g_in*sin_g_out);
    }
    Ith(ydot,5) = h_in_dot_newtonian;
    Ith(ydot,6) = h_in_dot_newtonian;
   
   
   
    /*******************************
     * a_in_dot                    *
     * *****************************/

	double a_in_dot_GR_25PN_in = 0.0;
    double a_in_dot_tides = 0.0;    
    double a_in_dot_wind=0.0;
    double a_in_dot_mass_transfer = 0.0;    
    
    /* post-Newtonian point particle */  
    if (include_25PN_inner_terms == TRUE)
    {
        a_in_dot_GR_25PN_in = -c_64div5*CONST_G_P3*m1_times_m2*m1_plus_m2*(1.0/(CONST_C_LIGHT_P5*a_in_p3*l_in_p7))*f_25PN_a_in;
    }

    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        double a_in_dot_tides_star1 = -6.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency1_div_n_in);
        double a_in_dot_tides_star2 = -6.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency2_div_n_in);            
        a_in_dot_tides = a_in_dot_tides_star1 + a_in_dot_tides_star2;
//        printf("a_in_dot_tides %g %g %g \n",a_in_dot_tides, a_in_dot_tides_star1, a_in_dot_tides_star2);
    }

    /* wind mass loss */
    if (include_inner_wind_terms == TRUE)
    {
        /* here, it is assumed that mass lost from star3 does not affect the inner orbit */
        double a_in_dot_wind_star1 = f_a_dot_mass_variations_fast_and_isotropic_wind(m1,wind_mass_loss_rate_star1,m2,a_in,inner_accretion_efficiency_wind_star1_to_star2);
        double a_in_dot_wind_star2 = f_a_dot_mass_variations_fast_and_isotropic_wind(m2,wind_mass_loss_rate_star2,m1,a_in,inner_accretion_efficiency_wind_star2_to_star1);
        a_in_dot_wind = a_in_dot_wind_star1 + a_in_dot_wind_star2;
    }
    
    /* mass transfer */
    if (include_inner_RLOF_terms == TRUE)
    {
        bool calculate = TRUE;
        double m_donor,m_accretor;
        if (star1_is_donor == TRUE)
        {
            m_donor = m1;
            m_accretor = m2;
        }
        else if (star2_is_donor == TRUE)
        {
            m_donor = m2;
            m_accretor = m1;
        }
        else
        {
            calculate = FALSE;
        }
        
        if (calculate == TRUE)
        {
            a_in_dot_mass_transfer = f_a_dot_mass_variations(m_donor,inner_mass_transfer_rate,m_accretor,a_in,inner_accretion_efficiency_mass_transfer,inner_specific_AM_loss_mass_transfer);
        }
    }
    
    /* combined */
    double a_in_dot = a_in_dot_GR_25PN_in + a_in_dot_tides + a_in_dot_wind + a_in_dot_mass_transfer;
	Ith(ydot,7) = a_in_dot;

printf("a_in_dot %g %g %g %g %g %g\n",a_in,a_in_dot, a_in_dot_GR_25PN_in, a_in_dot_tides, a_in_dot_wind, a_in_dot_mass_transfer);

    /*******************************
     * a_out_dot                   *
     * *****************************/

    double a_out_dot_GR_25PN_out = 0.0;
    double a_out_dot_tides = 0.0;
    double a_out_dot_wind = 0.0;
    double a_out_dot_mass_transfer = 0.0;
    
    if (ignore_tertiary == FALSE)
    {
        /* post-Newtonian point particle */  
        if (include_25PN_outer_terms == TRUE)
        {
            a_out_dot_GR_25PN_out = -c_64div5*CONST_G_P3*m1_plus_m2*m3*m1_plus_m2_plus_m3*(1.0/(CONST_C_LIGHT_P5*a_out_p3*l_out_p7))*f_25PN_a_out;
        }
    
        /* tides */
        if (include_outer_tidal_terms == TRUE)
        {
            double a_out_dot_tides_star3 = -6.0*(1.0+m1_plus_m2_div_m3)*m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3*R3_div_a_out_p2*a_out*pow(l_out,-15.0)*(f_tides1_out \
                - l_out_p3*f_tides2_out*spin_angular_frequency3_div_n_out);
            a_out_dot_tides = a_out_dot_tides_star3;
        }
        
        /* wind mass loss */
        if (include_outer_wind_terms == TRUE)
        {
            /* a_out_dot due to wind mass loss of star3, where part is accreted by the inner binary (treated as point mass) */
            double a_out_dot_wind_star3 = f_a_dot_mass_variations_fast_and_isotropic_wind(m3,wind_mass_loss_rate_star3,m1+m2,a_out,outer_accretion_efficiency_wind_star3_to_inner_binary);
    
            /* a_out_dot due to wind mass loss from the inner binary (treated as a point mass) to star3 */
            /* here, it is assumed that wind material from the inner binary can be accreted by star3, but not material due to nonconservative mass transfer in the inner binary */
            double m_dot_lost_from_inner_binary = \
                + wind_mass_loss_rate_star1*(1.0-inner_accretion_efficiency_wind_star1_to_star2) \
                + wind_mass_loss_rate_star2*(1.0-inner_accretion_efficiency_wind_star2_to_star1) \
                + inner_mass_transfer_rate*(1.0-inner_accretion_efficiency_mass_transfer);
            double m_dot_accreted_by_star3_from_inner_binary_winds = \
                + wind_mass_loss_rate_star1*(1.0-inner_accretion_efficiency_wind_star1_to_star2)*outer_accretion_efficiency_wind_star1_to_star3 \
                + wind_mass_loss_rate_star2*(1.0-inner_accretion_efficiency_wind_star2_to_star1)*outer_accretion_efficiency_wind_star2_to_star3;
    
            double a_out_dot_wind_inner_binary = 0.0;
            if (m_dot_lost_from_inner_binary < 0.0)
            {
                double effective_accretion_efficiency_wind_inner_binary_to_star3 = m_dot_accreted_by_star3_from_inner_binary_winds/m_dot_lost_from_inner_binary;
    
                a_out_dot_wind_inner_binary = f_a_dot_mass_variations_fast_and_isotropic_wind(m1+m2,m_dot_lost_from_inner_binary,m3,a_out,effective_accretion_efficiency_wind_inner_binary_to_star3);
            }
            a_out_dot_wind = a_out_dot_wind_star3 + a_out_dot_wind_inner_binary;
        }
        
        /* mass transfer */
        if (include_outer_RLOF_terms == TRUE)
        {
            bool calculate = TRUE;
            double m_donor,m_accretor;
            if (star3_is_donor == TRUE)
            {
                m_donor = m3;
                m_accretor = m1+m2;
                a_out_dot_mass_transfer = f_a_dot_mass_variations(m_donor,outer_mass_transfer_rate,m_accretor,a_out,outer_accretion_efficiency_mass_transfer,outer_specific_AM_loss_mass_transfer);
            }
        }
    }
    
    /* combined */
    double a_out_dot = a_out_dot_GR_25PN_out + a_out_dot_tides + a_out_dot_wind + a_out_dot_mass_transfer;
	Ith(ydot,8) = a_out_dot;



    /**********************************************
     * cositot_dot                                *
     * due to dynamical triple interaction only!  *
     * ********************************************/

    double cositot_dot = 0.0;
    if (ignore_tertiary == FALSE)
    {
        double G_in_dot = -G_in*e_in*(e_in_dot_newtonian+e_in_dot_GR_1PN_in_out)/l_in_p2;
        double G_out_dot = -G_out*e_out*e_out_dot_newtonian/l_out_p2;
        cositot_dot = (-1.0/(G_in*G_out))*(G_in_dot*(G_in + G_out*cositot) + G_out_dot*(G_out + G_in*cositot));
    }
    
	Ith(ydot,9) = cositot_dot;



    /************************************
     * spin_angular_frequency1_dot      *
     * **********************************/

	double spin_angular_frequency1_dot_tides = 0.0;
    double spin_angular_frequency1_dot_magnetic_braking = 0.0;
    double spin_angular_frequency1_dot_moment_of_inertia_plus_wind_changes = 0.0;
    double spin_angular_frequency1_dot_accretion_from_companions = 0.0;
    
    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        spin_angular_frequency1_dot_tides = 3.0*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*(m2_div_m1/(gyration_radius_star1*gyration_radius_star1)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency1_div_n_in);
    }

    /* magnetic braking */
    if (include_magnetic_braking_terms == TRUE)
    {
        spin_angular_frequency1_dot_magnetic_braking = spin_angular_frequency_dot_magnetic_braking(spin_angular_frequency1,m1,wind_mass_loss_rate_star1,gyration_radius_star1,R1,R1_dot);
    }
    
    /* changes of moment inertia, wind & AM accretion */
    if (include_spin_radius_mass_coupling_terms_star1 == TRUE)
    {
        /* changes in the spin frequency due to 1) changes in the moment of inertia and 2) loss of spin AM in the wind */
        spin_angular_frequency1_dot_moment_of_inertia_plus_wind_changes = spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes(spin_angular_frequency1,m1,R1,moment_of_inertia_star1,moment_of_inertia_dot_star1, wind_mass_loss_rate_star1,threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero);
        
        /* changes in the spin frequency due to spin AM transferred from the companion -- here, a `spin AM accretion efficiency' of 1 is assumed */
        /* also, it is assumed that star1 can only receive spin AM from star 2, not star3 */
        double m_dot_accreted_by_star1_from_wind_of_star2 = -inner_accretion_efficiency_wind_star2_to_star1*wind_mass_loss_rate_star2;
        spin_angular_frequency1_dot_accretion_from_companions = spin_angular_frequency_dot_accretion_from_companion(moment_of_inertia_star1, m_dot_accreted_by_star1_from_wind_of_star2, spin_angular_frequency2, R2);
//        printf("dots %g %g %g\n",spin_angular_frequency1,spin_angular_frequency1_dot_moment_of_inertia_plus_wind_changes,m_dot_accreted_by_star1_from_wind_of_star2);
    }

    double spin_angular_frequency1_dot = spin_angular_frequency1_dot_tides + spin_angular_frequency1_dot_magnetic_braking + spin_angular_frequency1_dot_moment_of_inertia_plus_wind_changes + spin_angular_frequency1_dot_accretion_from_companions;
	Ith(ydot,10) = spin_angular_frequency1_dot;
    


    /************************************
     * spin_angular_frequency2_dot      *
     * **********************************/

	double spin_angular_frequency2_dot_tides = 0.0;
    double spin_angular_frequency2_dot_magnetic_braking = 0.0;
    double spin_angular_frequency2_dot_moment_of_inertia_plus_wind_changes = 0.0;
    double spin_angular_frequency2_dot_accretion_from_companions = 0.0;
    
    /* tides */
    if (include_inner_tidal_terms == TRUE)
    {
        spin_angular_frequency2_dot_tides = 3.0*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*(m1_div_m2/(gyration_radius_star2*gyration_radius_star2)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency2_div_n_in);
    }
    
    /* magnetic braking */
    if (include_magnetic_braking_terms == TRUE)
    {
        spin_angular_frequency2_dot_magnetic_braking = spin_angular_frequency_dot_magnetic_braking(spin_angular_frequency2,m2,wind_mass_loss_rate_star2,gyration_radius_star2,R2,R2_dot);
    }

    /* changes of moment inertia, wind & AM accretion */
    if (include_spin_radius_mass_coupling_terms_star2 == TRUE)
    {
        /* changes in the spin frequency due to 1) changes in the moment of inertia and 2) loss of spin AM in the wind */
        spin_angular_frequency2_dot_moment_of_inertia_plus_wind_changes = spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes(spin_angular_frequency2,m2,R2,moment_of_inertia_star2,moment_of_inertia_dot_star2, wind_mass_loss_rate_star2,threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero);

        /* changes in the spin frequency due to spin AM transferred from the companion -- here, a `spin AM accretion efficiency' of 1 is assumed */
        /* also, it is assumed that star2 can only receive spin AM from star 1, not star3 */
        double m_dot_accreted_by_star2_from_wind_of_star1 = -inner_accretion_efficiency_wind_star1_to_star2*wind_mass_loss_rate_star1;
        spin_angular_frequency2_dot_accretion_from_companions = spin_angular_frequency_dot_accretion_from_companion(moment_of_inertia_star2, m_dot_accreted_by_star2_from_wind_of_star1, spin_angular_frequency1, R1);
    }
    
    double spin_angular_frequency2_dot = spin_angular_frequency2_dot_tides + spin_angular_frequency2_dot_magnetic_braking + spin_angular_frequency2_dot_moment_of_inertia_plus_wind_changes + spin_angular_frequency2_dot_accretion_from_companions;
	Ith(ydot,11) = spin_angular_frequency2_dot;


    /************************************
     * spin_angular_frequency3_dot      *
     * **********************************/

	double spin_angular_frequency3_dot_tides = 0.0;
    double spin_angular_frequency3_dot_magnetic_braking = 0.0;
    double spin_angular_frequency3_dot_moment_of_inertia_plus_wind_changes = 0.0;
    double spin_angular_frequency3_dot_accretion_from_companions = 0.0;

    if (ignore_tertiary == FALSE)
    {
        /* tides */
        if (include_outer_tidal_terms == TRUE)
        {
            spin_angular_frequency3_dot_tides = 3.0*m1_plus_m2_div_m3_times_R3_div_a_out_p6_times_k_div_T_tides_star3*(m1_plus_m2_div_m3/(gyration_radius_star3*gyration_radius_star3)) \
                *(n_out/(l_out_p6*l_out_p6))*(f_tides2_out - l_out_p3*f_tides5_out*spin_angular_frequency3_div_n_out);
        }
        
        /* magnetic braking */
        if (include_magnetic_braking_terms == TRUE)
        {
            spin_angular_frequency3_dot_magnetic_braking = spin_angular_frequency_dot_magnetic_braking(spin_angular_frequency3,m3,wind_mass_loss_rate_star3,gyration_radius_star3,R3,R3_dot);
        }
        
        /* changes of moment inertia, wind & AM accretion */
        if (include_spin_radius_mass_coupling_terms_star3 == TRUE)
        {
            /* changes in the spin frequency due to 1) changes in the moment of inertia and 2) loss of spin AM in the wind */        
            spin_angular_frequency3_dot_moment_of_inertia_plus_wind_changes = spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes(spin_angular_frequency3,m3,R3,moment_of_inertia_star3,moment_of_inertia_dot_star3, wind_mass_loss_rate_star3,threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero);
    
            /* changes in the spin frequency due to spin AM transferred from the companions in the inner orbit -- here, a `spin AM accretion efficiency' of 1 is assumed */
            double m_dot_accreted_by_star3_from_wind_of_star1 = -(1.0-inner_accretion_efficiency_wind_star1_to_star2)*outer_accretion_efficiency_wind_star1_to_star3*wind_mass_loss_rate_star1;
            double m_dot_accreted_by_star3_from_wind_of_star2 = -(1.0-inner_accretion_efficiency_wind_star2_to_star1)*outer_accretion_efficiency_wind_star2_to_star3*wind_mass_loss_rate_star2;
            
            spin_angular_frequency3_dot_accretion_from_companions += spin_angular_frequency_dot_accretion_from_companion(moment_of_inertia_star3, m_dot_accreted_by_star3_from_wind_of_star1, spin_angular_frequency1, R1);
            spin_angular_frequency3_dot_accretion_from_companions += spin_angular_frequency_dot_accretion_from_companion(moment_of_inertia_star3, m_dot_accreted_by_star3_from_wind_of_star2, spin_angular_frequency2, R2);
            
        }
    }
    
    double spin_angular_frequency3_dot = spin_angular_frequency3_dot_tides + spin_angular_frequency3_dot_magnetic_braking + spin_angular_frequency3_dot_moment_of_inertia_plus_wind_changes + spin_angular_frequency3_dot_accretion_from_companions;
	Ith(ydot,12) = spin_angular_frequency3_dot;
    
//#ifdef VERBOSE
//    if (fabs(spin_angular_frequency3_dot) > 1.0e4)
//    {
//        printf("dots -- e_in %g e_out %g g_in %g g_out %g a_in %g a_out %g spin_1 %g spin_2 %g spin_3 %g\n", e_in_dot,e_out_dot,g_in_dot,g_out_dot,a_in_dot,a_out_dot,spin_angular_frequency1_dot,spin_angular_frequency2_dot,spin_angular_frequency3_dot);
//        printf("dots -- spin_3 parts %g %g %g %g\n",spin_angular_frequency3_dot_tides,spin_angular_frequency3_dot_magnetic_braking,spin_angular_frequency3_dot_moment_of_inertia_plus_wind_changes,spin_angular_frequency3_dot_accretion_from_companions);
//        printf("arguments for spin_angular_frequency3_dot_moment_of_inertia_plus_wind_changes %g %g %g %g %g %g\n",spin_angular_frequency3,m3,R3,moment_of_inertia_star3,moment_of_inertia_dot_star3, wind_mass_loss_rate_star3);
//    }
//endif

	return 0;
}

/***************************************************
 * separate smaller functions for right-hand sides *
 ***************************************************/

/* tides (1981A&A....99..126H) */
double f_tides1(double e_p2)
{
    return 1.0 + e_p2*(c_31div2 + e_p2*(c_255div8 + e_p2*(c_185div16 + e_p2*c_25div64)));
}
double f_tides2(double e_p2)
{
    return 1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16));
}
double f_tides3(double e_p2)
{
    return 1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64));
}
double f_tides4(double e_p2)
{
    return 1.0 + e_p2*(c_3div2 + e_p2*c_1div8);
}
double f_tides5(double e_p2)
{
    return 1.0 + e_p2*(3.0 + e_p2*c_3div8);
}
double f_25PN_e(double e_p2)
{
    return 1.0 + e_p2*c_121div304;
}
double f_25PN_a(double e_p2)
{
    return 1.0 + e_p2*(c_73div24 + e_p2*c_37div96);
}

/* magnetic braking */
double spin_angular_frequency_dot_magnetic_braking(double spin_angular_frequency, double mass, double wind_mass_loss_rate, double gyration_radius, double radius, double radius_dot)
{
//    return spin_angular_frequency*( (wind_mass_loss_rate/mass)*(-1.0 + c_2div3/(gyration_radius*gyration_radius)) - 2.0*(radius_dot/radius) );
    return spin_angular_frequency*(wind_mass_loss_rate/mass)*c_2div3/(gyration_radius*gyration_radius);
}

/* effect of mass & radius changes on stellar spin -- assuming that the spin angular momentum remains constant */
double spin_angular_frequency_dot_mass_radius_changes(double spin_angular_frequency, double m, double m_dot, double R, double R_dot)
{
    return -spin_angular_frequency*(m_dot/m + 2.0*R_dot/R);
}

/* effect of moment of inertia changes & wind mass loss on stellar spin */
double spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes(double spin_angular_frequency, double mass, double radius, double moment_of_inertia, double moment_of_inertia_dot, double m_dot_wind, double threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero)
{
    double spin_angular_frequency_dot = spin_angular_frequency*(c_2div3*m_dot_wind*radius*radius/moment_of_inertia - moment_of_inertia_dot/moment_of_inertia);
    if ((spin_angular_frequency <= threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero) && (spin_angular_frequency_dot < 0.0))
    {
//        printf("threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero %g %g \n",spin_angular_frequency,threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero);
        spin_angular_frequency_dot = 0.0;
    }
    
    return spin_angular_frequency_dot;
}

/* effect wind mass accretion from companion on spin */
double spin_angular_frequency_dot_accretion_from_companion(double moment_of_inertia, double m_dot_accretion, double companion_spin_angular_frequency, double companion_radius)
{
    double spin_angular_momentum_wind_accretion_efficiency_companion_to_primary = 1.0; /* working assumption */
    
    return companion_spin_angular_frequency*c_2div3*spin_angular_momentum_wind_accretion_efficiency_companion_to_primary*m_dot_accretion*companion_radius*companion_radius/moment_of_inertia;
}

/* effect of mass variations on orbit */
double f_a_dot_mass_variations(double m_donor, double m_donor_dot, double m_accretor, double a, double beta, double gamma)
{
    /* Lectore notes on binary evolution by Onno Pols -- chapter 7 Eq. 7.14 */
    return -2.0*a*(m_donor_dot/m_donor)*(1.0 - beta*(m_donor/m_accretor) - (1.0-beta)*(gamma + c_1div2)*m_donor/(m_donor+m_accretor) );
}

double f_a_dot_mass_variations_fast_and_isotropic_wind(double m_donor, double m_donor_dot, double m_accretor, double a, double beta)
{
    /* Lectore notes on binary evolution by Onno Pols -- chapter 7 Eq. 7.14 */
    double m_accretor_dot = -beta*m_donor_dot;
    double m_total = m_donor + m_accretor;
    double m_total_dot = (1.0 - beta)*m_donor_dot;
    return a*(-m_total_dot/m_total - 2.0*m_accretor_dot/m_accretor - 2.0*beta*m_donor_dot/m_donor);
}

/* root finding functions */
int froot_delaunay(realtype t, N_Vector yev, realtype *gout, void *data_f)
{
	UserData data;
	data = (UserData) data_f;

    /* by default checks are not done (reduces computational cost) */
    /* setting gout[i] to >0 means no root will be found */
    gout[0] = 1e10;
    gout[1] = 1e10;
    gout[2] = 1e10;
    gout[3] = 1e10;
    gout[4] = 1e10;
    gout[5] = 1e10;
    gout[6] = 1e10;

    bool check_for_dynamical_stability = data->check_for_dynamical_stability;
    bool check_for_inner_collision = data->check_for_inner_collision;
    bool check_for_outer_collision = data->check_for_outer_collision;
    bool check_for_inner_RLOF = data->check_for_inner_RLOF;
    bool check_for_outer_RLOF = data->check_for_outer_RLOF;
    
    
    /* addition Jan 2017 */
    bool check_for_semisecular_regime = data->check_for_semisecular_regime;

    double m1 = data->m1;
    double m2 = data->m2;				
    double m3 = data->m3;		
    
    double R1 = data->R1;
    double R2 = data->R2;
    double R3 = data->R3;

    double a_in = Ith(yev,7);
    double a_out = Ith(yev,8);
    double e_in = 1.0 - pow(10.0,Ith(yev,1)); /* Current e1 */
    double e_out = 1.0 - pow(10.0,Ith(yev,2)); /* Current e2 */

    double rp_in = a_in*(1.0 - e_in);
    double rp_out = a_out*(1.0 - e_out);
    
    int roche_radius_specification = data->roche_radius_specification;
    int stability_limit_specification = data->stability_limit_specification;
    
    if (check_for_dynamical_stability == TRUE)
    {
        /*	check for dynamical stability */
        double itot = acos(Ith(yev,9));
        double a_out_div_a_in_crit = a_out_div_a_in_dynamical_stability(m1,m2,m3,a_in,a_out,e_in,e_out,itot,stability_limit_specification);
        gout[0] = a_out/a_in - a_out_div_a_in_crit;///1.5;
        printf("dyn stab %g %g %g %g \n",a_out, a_in, a_out_div_a_in_crit, gout[0]);        
    }

    if (check_for_inner_collision == TRUE)
    {
        /*	check for collision at periastron (inner binary)	*/
        gout[1] = rp_in - (R1 + R2);
        printf("collision %g %g %g %g %g\n",a_in, rp_in, R1, R2, gout[1]);
        
    }
    
    if (check_for_outer_collision == TRUE)
    {
        /*	check for "collision" at periastron (outer binary)	*/
        /*  the "radius" of the inner binary is set to the inner apocenter distance */
        double ra_in = a_in*(1.0 + e_in);
        gout[2] = rp_out - (R3 + ra_in);
    }
    
    if (check_for_inner_RLOF == TRUE)
    {
        double spin_angular_frequency1 = Ith(yev,10);
        double spin_angular_frequency2 = Ith(yev,11);
        double spin_angular_frequency_inner_orbit_periapse = sqrt( CONST_G*(m1+m2)*(1.0+e_in)/(rp_in*rp_in*rp_in) );
        double f1 = spin_angular_frequency1/spin_angular_frequency_inner_orbit_periapse;
        double f2 = spin_angular_frequency2/spin_angular_frequency_inner_orbit_periapse;
        
        double roche_radius_pericenter_inner_star1 = roche_radius(rp_in, m1/m2, e_in, f1, roche_radius_specification);
        gout[3] = R1 - roche_radius_pericenter_inner_star1;
        printf("RLOF %g %g %g \n",R1, roche_radius_pericenter_inner_star1, gout[3]);

        double roche_radius_pericenter_inner_star2 = roche_radius(rp_in, m2/m1, e_in, f2, roche_radius_specification);

        gout[4] = R2 - roche_radius_pericenter_inner_star2;
        printf("RLOF2 %g %g %g \n",R2, roche_radius_pericenter_inner_star2, gout[4]);
        
    }
    if (check_for_outer_RLOF == TRUE)
    {

        double spin_angular_frequency3 = Ith(yev,12);
        double spin_angular_frequency_outer_orbit_periapse = sqrt( CONST_G*(m1+m2+m3)*(1.0+e_out)/(rp_out*rp_out*rp_out) );        
        double f3 = spin_angular_frequency3/spin_angular_frequency_outer_orbit_periapse;
        
        double roche_radius_pericenter_outer_star3 = roche_radius(rp_out, m3/(m1+m2), e_out, f3, roche_radius_specification);
        gout[5] = R3 - roche_radius_pericenter_outer_star3;
    }            
    
    
    
    /* addition Jan/Feb 2017 */
    if (check_for_semisecular_regime == TRUE)
    {
        double check_for_semisecular_regime_parameter = data->check_for_semisecular_regime_parameter;
        double a_out_div_a_in_crit = a_out_div_a_in_semisecular_regime(m1,m2,m3,e_in,e_out,check_for_semisecular_regime_parameter);
//        gout[6] = fabs(a_out/a_in - a_out_div_a_in_crit;
        gout[6] = a_out/a_in - a_out_div_a_in_crit;

            printf("sil1 %.16g %.16g %g %g %.16g\n",a_out_div_a_in_crit, a_out/a_in, a_out, a_in, t);
//            printf("sil2 %g %g %g \n",t_dif,t, global_time_step);

    }
   
	return 0;
}

double a_out_div_a_in_semisecular_regime(double m1, double m2, double m3, double e_in, double e_out, double check_for_semisecular_regime_parameter)
{
    /* also used in interface.py */
    /* 2014ApJ...781...45A Eq. (18) */
            
    double a_out_div_a_in_crit = (1.0/(1.0-e_out))*pow( check_for_semisecular_regime_parameter*5.0*M_PI*(m3/(m1+m2))*(1.0/sqrt(1.0-e_in)), 1.0/3.0);
    return a_out_div_a_in_crit;
}

double a_out_div_a_in_dynamical_stability(double m1, double m2, double m3, double a_in, double a_out, double e_in, double e_out, double itot, int stability_limit_specification)
{
    /* wrapper used in interface.py */
    if (stability_limit_specification == 1){
        return a_out_div_a_in_dynamical_stability_petrovich_15_simple(e_in,e_out);}
    else if (stability_limit_specification == 2){
        return a_out_div_a_in_dynamical_stability_petrovich_15(m1,m2,m3,e_in,e_out,a_in,a_out);}
    else if (stability_limit_specification == 3){
        return a_out_div_a_in_dynamical_stability_holman_stype_98(m1,m2,m3,e_out);}        
    else if (stability_limit_specification == 4){
        return a_out_div_a_in_dynamical_stability_holman_ptype_98(m1,m2,m3,e_in);}         
    else if (stability_limit_specification == 5){
        return a_out_div_a_in_dynamical_stability_vynatheya(m1,m2,m3,e_in,e_out,itot);}         
    else if (stability_limit_specification == 6){
        return a_out_div_a_in_dynamical_stability_tory(m1,m2,m3,e_in,e_out,itot);}         
    else {
        return a_out_div_a_in_dynamical_stability_mardling_aarseth_01(m1,m2,m3,e_out,itot);}

}

double a_out_div_a_in_dynamical_stability_mardling_aarseth_01(double m1, double m2, double m3, double e_out, double itot)
{
    /* Mardling & Aarseth criterion (2001MNRAS.321..398M) 
     * including the `ad hoc' inclination factor
     * itot is assumed to be in radians! */
     
    double q_out = m3/(m1+m2);
    double a_out_div_a_in_crit = (1.0/(1.0-e_out))*2.8*pow( (1.0+q_out)*(1.0+e_out)/sqrt(1.0-e_out),c_2div5)*(1.0 - 0.3*itot/M_PI);
    
    return a_out_div_a_in_crit;
}



double a_out_div_a_in_dynamical_stability_petrovich_15_simple(double e_in,  double e_out)
{
    /* Petrovich criterion (2015ApJ...808..120P) 
       for star + 2 planets */     
//     printf("Petrovich \n");

    double a_out_div_a_in_crit = 1.83 * (1+e_in)/(1-e_out);        
    return a_out_div_a_in_crit;
}


double a_out_div_a_in_dynamical_stability_petrovich_15(double m1, double m2, double m3,double e_in,  double e_out, double a_in,  double a_out)
{
    /* Petrovich criterion (2015ApJ...808..120P) 
       for star + 2 planets */ 
//     printf("Petrovich \n");

    double mu_in = min(m1,m2)/max(m1,m2);
    double mu_out = m3/max(m1,m2);
    double max_mu = max(mu_in, mu_out);
    double a_out_div_a_in_crit = (1+e_in)/(1-e_out) * (2.4*pow(max_mu,1./3.)*pow(a_out/a_in, 0.5)+1.15);        
    return a_out_div_a_in_crit;
}

double a_out_div_a_in_dynamical_stability_holman_stype_98(double m1, double m2, double m3, double e_out)
{
    /* Holman criterion ( 1999AJ....117..621H) 
        for planet orbiting a star (s+p)+s 
        updated by Quarles et al. 2020 2020AJ....159...80Q for i=0, 30,45 & 180
        */     
     printf("Holman S-type \n");

    double mu = m3/(max(m1,m2)+m3);
    double a_out_div_a_in_crit = 1./(0.464-0.38*mu-0.631*e_out+0.586*mu*e_out+0.15*e_out*e_out-0.198*mu*e_out*e_out);
    return a_out_div_a_in_crit;
}

double a_out_div_a_in_dynamical_stability_holman_ptype_98(double m1, double m2, double m3, double e_in)
{
    /* Holman criterion ( 1999AJ....117..621H) 
         for planet orbiting a binary (s+s)+p 
         updated by Quarles et al. 2018 2018ApJ...856..150Q 
         */     
     printf("Holman P-type \n");

    double mu = min(m1,m2)/(m1+m2);
    double a_out_div_a_in_crit = 1.6+5.1*e_in-2.22*e_in*e_in+4.12*mu-4.27*e_in*mu-5.09*mu*mu+4.61*e_in*e_in*mu*mu;
    return a_out_div_a_in_crit;
}

double a_out_div_a_in_dynamical_stability_vynatheya(double m1, double m2, double m3, double e_in, double e_out, double itot)
{
    /* Vynatheya criterion (2022MNRAS.516.4146V) */     
     printf("Vynatheya \n");

    double q_out = m3/(m1+m2);
    double cositot = cos(itot);
    double e_in_max = pow(1-5./3.*cositot*cositot,0.5);
    double e_in_av = max(e_in, 0.5*e_in_max*e_in_max);
    double Ycrit = 2.4 * pow((1+q_out)/(1+e_in_av)/pow(1-e_out,0.5), 0.4) * ((1-0.2*e_in_av+e_out)/8.*(cositot -1)+1);
    double a_out_div_a_in_crit = Ycrit * (1+e_in_av)/(1-e_out);
//    double a_out_div_a_in_crit = Ycrit * (1+e_in??)/(1-e_out);
    return a_out_div_a_in_crit;
}

double a_out_div_a_in_dynamical_stability_tory(double m1, double m2, double m3, double e_in, double e_out, double itot)
{
    /* Tory criterion (2022PASA...39...62T) */     
     printf("Tory \n");

    double q_out = (m1+m2)/m3;
    double f = pow(10, -0.6+0.04*q_out) * pow(q_out, 0.32+0.1*q_out);
    double g = -0.4*cos(itot)+1.4;
    if (itot>M_PI/3.){
        g = -0.1773*pow(itot,4) + 1.1211*pow(itot,3) -1.9149*itot*itot + 0.5022*itot + 1.6222;
    }
    double logh = -1.*itot*pow(q_out,1.3)/1500;
    
    double a_out_div_a_in_crit = 1./(f*g*pow(10,logh)) / (1-e_out);
    return a_out_div_a_in_crit;
}



double roche_radius(double rp, double q, double e, double f, int roche_radius_specification)
{   
    
    if (roche_radius_specification == 1){ //printf("roche radius of sepinsky \n");
        return roche_radius_pericenter_sepinsky(rp, q, e, f);}
    else if (roche_radius_specification == 2){ //printf("classical roche radius of eggleton for circular binaries \n");  
        return roche_radius_pericenter_eggleton(rp/(1.-e), q);} 
    else { //printf("roche radius of eggleton with eccentricity factor \n");
    return roche_radius_pericenter_eggleton(rp, q);}

}



double roche_radius_pericenter_eggleton(double rp, double q)
{
    /* 2007ApJ...660.1624S Eqs. (45) */    
    /* q is defined as m_primary/m_secondary */
    double q_pow_one_third = pow(q,c_1div3);
    double q_pow_two_third = q_pow_one_third*q_pow_one_third;
    return rp*0.49*q_pow_two_third/(0.6*q_pow_two_third + log(1.0 + q_pow_one_third));
}

double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f)
{
    /* 2007ApJ...660.1624S Eqs. (47)-(52) */
    double log_q = log10(q);
    double A = f*f*(1.0 + e); // assumes pericenter
    double log_A = log10(A);

    double R_L_pericenter_eggleton = roche_radius_pericenter_eggleton(rp,q);
    double ratio = 0.0; // this is R_L divided by R_L_pericenter_eggleton

    if (log_q < 0.0)
    {
        if (log_A <= -0.1)
        {
            double c = 0.5*(1.0+A) + log_q;
            ratio = 1.0 + 0.11*(1.0-A) - 0.05*(1.0-A)*exp(-c*c);
        }
        if ((log_A > -0.1) && (log_A < 0.2))
        {
            double g_0 = 0.9978 - 0.1229*log_A - 0.1273*log_A*log_A;
            double g_1 = 0.001 + 0.02556*log_A;
            double g_2 = 0.0004 + 0.0021*log_A;
            ratio = g_0 + g_1*log_q * g_2*log_q*log_q;
        }
        if (log_A >= 0.2)
        {
            double num_0 = 6.3014*pow(log_A,1.3643);
            double den_0 = exp(2.3644*pow(log_A,0.70748)) - 1.4413*exp(-0.0000184*pow(log_A,-4.5693));
            double i_0 = num_0/den_0;

            double den_1 = 0.0015*exp(8.84*pow(log_A,0.282)) + 15.78;
            double i_1 = log_A/den_1;

            double num_2 = 1.0 + 0.036*exp(8.01*pow(log_A,0.879));
            double den_2 = 0.105*exp(7.91*pow(log_A,0.879));
            double i_2 = num_2/den_2;

            double den_3 = 1.38*exp(-0.035*pow(log_A,0.76)) + 23.0*exp(-2.89*pow(log_A,0.76));
            double i_3 = 0.991/den_3;

            double c = log_q + i_3;
            ratio = i_0 + i_1*exp(-i_2*c*c);
        }
    }
    if (log_q >= 0.0)
    {
        if (log_A <= -0.1)
        {
            ratio = 1.226 - 0.21*A - 0.15*(1.0-A)*exp( (0.25*A - 0.3)*pow(log_q,1.55) );
        }
        if ((log_A > -0.1) && (log_A < 0.2))
        {
            double log_A_p2 = log_A*log_A;
            double h_0 = 1.0071 - 0.0907*log_A - 0.0495*log_A_p2;
            double h_1 = -0.004 - 0.163*log_A - 0.214*log_A_p2;
            double h_2 = 0.00022 - 0.0108*log_A - 0.02718*log_A_p2;
            ratio = h_0 + h_1*log_q + h_2*log_q*log_q;
        }
        if (log_A >= 0.2)
        {
            double num_0 = 1.895*pow(log_A,0.837);
            double den_0 = exp(1.636*pow(log_A,0.789)) - 1.0;
            double j_0 = num_0/den_0;

            double num_1 = 4.3*pow(log_A,0.98);
            double den_1 = exp(2.5*pow(log_A,0.66)) + 4.7;
            double j_1 = num_1/den_1;

            double den_2 = 8.8*exp(-2.95*pow(log_A,0.76)) + 1.64*exp(-0.03*pow(log_A,0.76));
            double j_2 = 1.0/den_2;

//            double j_3 = 0.256*exp(-1.33*pow(log_A,2.9))*( 5.5*exp(1.33*pow(log_A,2.9)) + 1.0 );
            double j_3 = 0.256*(5.5 + exp(-1.33*pow(log_A,2.9)));

            ratio = j_0 + j_1*exp(-j_2*pow(log_q,j_3));
            
//            printf("log_A %g\n",log_A);
//            printf("1 %g %g %g \n",num_0,den_0,j_0);
//            printf("2 %g %g %g \n",num_1,den_1,j_1);            
//            printf("2 %g %g %g \n",den_2,j_2,j_3);            
//            printf("ratio %g %g %g \n",ratio);            
        }
    }
    if (ratio == 0.0)
    {
        printf("unrecoverable error occurred in function roche_radius_pericenter_sepinsky in ODE_system.c\n");
        printf("log_q %g log_A %g ratio %g\n",log_q,log_A,ratio);
        printf("rp %g q %g e %g f %g\n",rp,q,e,f);
        exit(-1);
    }
    
    return ratio*R_L_pericenter_eggleton;
}
