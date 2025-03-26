int compute_effect_of_SN_on_orbital_vectors
(
    double m1, double m2, double m3,
    double e1_vec_x, double e1_vec_y, double e1_vec_z,
    double e2_vec_x, double e2_vec_y, double e2_vec_z,
    double h1_vec_x, double h1_vec_y, double h1_vec_z,
    double h2_vec_x, double h2_vec_y, double h2_vec_z,
    double Vkick1_vec_x, double Vkick1_vec_y, double Vkick1_vec_z,
    double Vkick2_vec_x, double Vkick2_vec_y, double Vkick2_vec_z,
    double Vkick3_vec_x, double Vkick3_vec_y, double Vkick3_vec_z,
    double delta_m1, double delta_m2, double delta_m3,
    double f1, double f2,
    double *V1_prime_x, double *V1_prime_y, double *V1_prime_z,
    double *V2_prime_x, double *V2_prime_y, double *V2_prime_z, 
    double *V3_prime_x, double *V3_prime_y, double *V3_prime_z,
    double *e1_vec_x_p, double *e1_vec_y_p, double *e1_vec_z_p,
    double *e2_vec_x_p, double *e2_vec_y_p, double *e2_vec_z_p,
    double *h1_vec_x_p, double *h1_vec_y_p, double *h1_vec_z_p,
    double *h2_vec_x_p, double *h2_vec_y_p, double *h2_vec_z_p,
    double *cos_phi1, double *cos_phi2, double *R, double *v_sys,
    double *r1dotv1, double *r1dotvk, double *v1dotvk,
    double *r1dotvsys, double *v1dotvsys, double *r2dotv2, double *r2dotvsys, double *v2dotvsys,
    double *r1dotr2, double *r1dotv2
);
int compute_eccentricity_vector(double total_mass, double r[3], double v[3], double e_vec[3]);
int compute_h_vector(double mu, double r[3], double v[3], double h_vec[3]);

int compute_orbital_vectors_from_orbital_elements(double m1, double m2, double m3,
    double a1, double a2, double e1, double e2,
    double INCL_rel,
    double AP1, double AP2, double LAN1, double LAN2,
    double *e1_vec_x, double *e1_vec_y, double *e1_vec_z,
    double *e2_vec_x, double *e2_vec_y, double *e2_vec_z,
    double *h1_vec_x, double *h1_vec_y, double *h1_vec_z,
    double *h2_vec_x, double *h2_vec_y, double *h2_vec_z);
int compute_orbital_elements_from_orbital_vectors(double m1, double m2, double m3, 
    double e1_vec_x,double e1_vec_y,double e1_vec_z,
    double e2_vec_x,double e2_vec_y,double e2_vec_z,
    double h1_vec_x,double h1_vec_y,double h1_vec_z,
    double h2_vec_x,double h2_vec_y,double h2_vec_z,
    double *a1, double *a2, double *e1, double *e2,
    double *INCL_rel,
    double *AP1, double *AP2, double *LAN1, double *LAN2);

void cross3(double a[3], double b[3], double result[3]);
double norm3(double v[3]);
double norm3_squared(double v[3]);
double dot3(double a[3], double b[3]);
