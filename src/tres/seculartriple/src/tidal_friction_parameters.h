bool check_for_radiative_damping(int stellar_type, double mass);
bool check_for_convective_damping(int stellar_type);
double set_crude_gyration_radii_based_on_stellar_structure(int stellar_type, double mass);

double compute_k_div_T_tides
(
    int stellar_type,
    double mass,
    double convective_envelope_mass,
    double companion_mass,
    double semimajor_axis,
    double radius,
    double convective_envelope_radius,
    double luminosity,
    double spin_angular_frequency,
    double gyration_radius,
    double amc
);
