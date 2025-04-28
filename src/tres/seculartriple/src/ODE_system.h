int fev_triad(realtype t, N_Vector yev, N_Vector ydot, void *data);
int fev_delaunay(realtype t, N_Vector yev, N_Vector ydot, void *data);
double f_tides1(double e_p2);
double f_tides2(double e_p2);
double f_tides3(double e_p2);
double f_tides4(double e_p2);
double f_tides5(double e_p2);
double f_25PN_e(double e_p2);
double f_25PN_a(double e_p2);

double spin_angular_frequency_dot_magnetic_braking(double spin_angular_frequency, double mass, double wind_mass_loss_rate, double gyration_radius, double radius, double radius_dot);
double spin_angular_frequency_dot_mass_radius_changes(double spin_angular_frequency, double m, double m_dot, double R, double R_dot);
double spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes(double spin_angular_frequency, double mass, double radius, double moment_of_inertia, double moment_of_inertia_dot, double m_dot_wind, double threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero);
double spin_angular_frequency_dot_accretion_from_companion(double moment_of_inertia, double m_dot_accretion, double companion_spin_angular_frequency, double companion_radius);

double f_a_dot_mass_variations(double m_donor, double m_donor_dot, double m_accretor, double a, double beta, double gamma);
double f_a_dot_mass_variations_fast_and_isotropic_wind(double m_donor, double m_donor_dot, double m_accretor, double a, double beta);
int froot_delaunay(realtype t, N_Vector yev, realtype *gout, void *data);

double a_out_div_a_in_semisecular_regime(double m1, double m2, double m3, double e_in, double e_out, double check_for_semisecular_regime_parameter);
double a_out_div_a_in_dynamical_stability(double m1, double m2, double m3, double a_in,double a_out,double e_in,double e_out, double itot, int stability_limit_specification);
double a_out_div_a_in_dynamical_stability_mardling_aarseth_01(double m1, double m2, double m3, double e_out, double itot);
double a_out_div_a_in_dynamical_stability_petrovich_15_simple(double e_in, double e_out);
double a_out_div_a_in_dynamical_stability_petrovich_15(double m1, double m2, double m3,double e_in,  double e_out, double a_in,  double a_out);
double a_out_div_a_in_dynamical_stability_holman_stype_98(double m1, double m2, double m3, double e_out);
double a_out_div_a_in_dynamical_stability_holman_ptype_98(double m1, double m2, double m3, double e_in);
double a_out_div_a_in_dynamical_stability_vynatheya(double m1, double m2, double m3, double e_in, double e_out, double itot);
double a_out_div_a_in_dynamical_stability_tory(double m1, double m2, double m3, double e_in, double e_out, double itot);


double roche_radius_pericenter_eggleton(double rp, double q);
double roche_radius_pericenter_sepinsky(double rp, double q, double e, double f);
double roche_radius(double rp, double q, double e, double f, int roche_radius_specification);

