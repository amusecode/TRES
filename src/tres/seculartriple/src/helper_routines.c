#include "main_code.h"

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
)
{
    double e1_vec[3] = {e1_vec_x,e1_vec_y,e1_vec_z};
    double e2_vec[3] = {e2_vec_x,e2_vec_y,e2_vec_z};
    double h1_vec[3] = {h1_vec_x,h1_vec_y,h1_vec_z};
    double h2_vec[3] = {h2_vec_x,h2_vec_y,h2_vec_z};
 
    double Vkick1_vec[3] = {Vkick1_vec_x,Vkick1_vec_y,Vkick1_vec_z};
    double Vkick2_vec[3] = {Vkick2_vec_x,Vkick2_vec_y,Vkick2_vec_z};
    double Vkick3_vec[3] = {Vkick3_vec_x,Vkick3_vec_y,Vkick3_vec_z};
    
    double e1_vec_unit[3],e2_vec_unit[3];
    double q1_vec[3],q2_vec[3];
    double q1_vec_unit[3],q2_vec_unit[3];
    
    double h1 = norm3(h1_vec);
    double h2 = norm3(h2_vec);
    double e1 = norm3(e1_vec);
    double e2 = norm3(e2_vec);

    double jsq1 = 1.0 - e1*e1;
    double jsq2 = 1.0 - e2*e2;
    
    double a1 = h1*h1*(m1+m2)/( CONST_G*m1*m1*m2*m2*jsq1 );
    double a2 = h2*h2*(m1+m2+m3)/( CONST_G*(m1+m2)*(m1+m2)*m3*m3*jsq2 );

    cross3(h1_vec,e1_vec,q1_vec);
    cross3(h2_vec,e2_vec,q2_vec);
    double q1 = norm3(q1_vec);
    double q2 = norm3(q2_vec);

    for (int i=0; i<3; i++)
    {
        e1_vec_unit[i] = e1_vec[i]/e1;
        e2_vec_unit[i] = e2_vec[i]/e2;
        q1_vec_unit[i] = q1_vec[i]/q1;
        q2_vec_unit[i] = q2_vec[i]/q2;
    }

    double cos_f1 = cos(f1);
    double sin_f1 = sin(f1);
    double cos_f2 = cos(f2);
    double sin_f2 = sin(f2);
    
    double r1 = a1*jsq1/(1.0 + e1*cos_f1);
    double r2 = a2*jsq2/(1.0 + e2*cos_f2);
    double v1 = sqrt(CONST_G*(m1+m2)/(a1*jsq1));
    double v2 = sqrt(CONST_G*(m1+m2+m3)/(a2*jsq2));
    double r1_vec[3],r2_vec[3];
    double v1_vec[3],v2_vec[3];

    /* '_p' denotes 'prime', i.e. after the SN */
    double r1_vec_p[3],r2_vec_p[3];
    double v1_vec_p[3],v2_vec_p[3];
    
    /* post-SN masses */
    /* note sign convention */
    double m1_p = m1 - delta_m1;
    double m2_p = m2 - delta_m2;
    double m3_p = m3 - delta_m3;
    
    double R1_vec[3],R2_vec[3],R3_vec[3];
    double V1_vec[3],V2_vec[3],V3_vec[3];
    double V1_vec_p[3],V2_vec_p[3],V3_vec_p[3];
    double r1_CM_vec[3],r1_CM_vec_p[3];
    double v1_CM_vec[3],v1_CM_vec_p[3];
    
    /* without loss of generality, set the initial CM of the triple to the origin */
    double r2_CM_vec[3] = {0.0,0.0,0.0};
    double v2_CM_vec[3] = {0.0,0.0,0.0};

    double v_sys_vec[3];
    
    for (int i=0; i<3; i++)
    {
        /* pre-SN */
        r1_vec[i] = r1*( cos_f1*e1_vec_unit[i] + sin_f1*q1_vec_unit[i] );
        r2_vec[i] = r2*( cos_f2*e2_vec_unit[i] + sin_f2*q2_vec_unit[i] );
        v1_vec[i] = v1*( -sin_f1*e1_vec_unit[i] + (e1+cos_f1)*q1_vec_unit[i] );
        v2_vec[i] = v2*( -sin_f2*e2_vec_unit[i] + (e2+cos_f2)*q2_vec_unit[i] );
    
        r1_CM_vec[i] = r2_CM_vec[i] + (m3/(m1+m2+m3))*r2_vec[i];
        v1_CM_vec[i] = v2_CM_vec[i] + (m3/(m1+m2+m3))*v2_vec[i];
        R3_vec[i] = r2_CM_vec[i] - ((m1+m2)/(m1+m2+m3))*r2_vec[i];
        V3_vec[i] = v2_CM_vec[i] - ((m1+m2)/(m1+m2+m3))*v2_vec[i];

        R1_vec[i] = r1_CM_vec[i] + (m2/(m1+m2))*r1_vec[i];
        V1_vec[i] = v1_CM_vec[i] + (m2/(m1+m2))*v1_vec[i];
        R2_vec[i] = r1_CM_vec[i] - (m1/(m1+m2))*r1_vec[i];
        V2_vec[i] = v1_CM_vec[i] - (m1/(m1+m2))*v1_vec[i];

        /* post-SN (positions are assumed to be unchanged) */
        V1_vec_p[i] = V1_vec[i] + Vkick1_vec[i];
        V2_vec_p[i] = V2_vec[i] + Vkick2_vec[i];
        V3_vec_p[i] = V3_vec[i] + Vkick3_vec[i];
        
        r1_CM_vec_p[i] = ((m1_p*R1_vec[i] + m2_p*R2_vec[i])/(m1_p+m2_p)); /* inner CM position changes because of mass loss */
        v1_CM_vec_p[i] = ((m1_p*V1_vec_p[i] + m2_p*V2_vec_p[i])/(m1_p+m2_p)); /* inner CM velocity changes because of mass loss and/or kick velocities */
        
        r1_vec_p[i] = r1_vec[i]; /* inner relative position vector does not change */
        r2_vec_p[i] = r1_CM_vec_p[i] - R3_vec[i]; /* RELATIVE position vector of outer binary changes because CM position of inner binary changes */
        v1_vec_p[i] = V1_vec_p[i] - V2_vec_p[i]; /* only changed by kicks */
        v2_vec_p[i] = v1_CM_vec_p[i] - V3_vec_p[i];
        
        //temp[i] = (delta_m1/(m1+m2-delta_m1))*(v1_CM_vec[i] - V1_vec[i]);
        //temp[i] = v1_CM_vec_p[i] - v1_CM_vec[i];
        //v_sys_vec[i] = v2_vec_p[i] - v2_vec[i];
        v_sys_vec[i] = v1_CM_vec_p[i] - v1_CM_vec[i];
        
    }

    *cos_phi1 = dot3(Vkick1_vec,v1_vec)/( norm3(Vkick1_vec)*norm3(v1_vec) );
    //*cos_phi2 = dot3(v1_CM_vec_p,v2_vec)/(norm3(v1_CM_vec_p)*norm3(v2_vec));
    *cos_phi2 = dot3(v_sys_vec,v2_vec)/(norm3(v_sys_vec)*norm3(v2_vec));
    
    *R = norm3(r2_vec_p);
    //*v_sys = norm3(v1_CM_vec_p);
    *v_sys = norm3(v_sys_vec);
    
    *r1dotv1 = dot3(r1_vec,v1_vec);
    *r1dotvk = dot3(r1_vec,Vkick1_vec);
    *v1dotvk = dot3(v1_vec,Vkick1_vec);

    *r1dotvsys = dot3(r1_vec,v_sys_vec);
    *v1dotvsys = dot3(v1_vec,v_sys_vec);

    *r2dotvsys = dot3(r2_vec,v_sys_vec);
    *v2dotvsys = dot3(v2_vec,v_sys_vec);
    
    *r1dotr2 = dot3(r1_vec,r2_vec);
    *r1dotv2 = dot3(r1_vec,v2_vec);
    
    *r2dotv2 = dot3(r2_vec,v2_vec);
//    printf("code V0 %g\n",norm3(v2_vec));
//    printf("R0 %g\n",norm3(r2_vec));
//    printf("R %g\n",norm3(r2_vec_p));
    
#ifdef IGNORE
    //printf("test %g\n",norm3(v1_CM_vec)/norm3(V1_vec)*(m1+m2));
    //printf("expr1 %g %g %g\n",v1_CM_vec_p[0]-v1_CM_vec[0],v1_CM_vec_p[1]-v1_CM_vec[1],v1_CM_vec_p[2]-v1_CM_vec[2]);
    //printf("expr2 %g %g %g\n",temp[0],temp[1],temp[2]);

    printf("pre rv1 %g %g %g %g %g %g\n",r1_vec[0],r1_vec[1],r1_vec[2],v1_vec[0],v1_vec[1],v1_vec[2]);
    printf("pre rv2 %g %g %g %g %g %g\n",r2_vec[0],r2_vec[1],r2_vec[2],v2_vec[0],v2_vec[1],v2_vec[2]);
    printf("pre rvcm1 %g %g %g %g %g %g\n",r1_CM_vec[0],r1_CM_vec[1],r1_CM_vec[2],v1_CM_vec[0],v1_CM_vec[1],v1_CM_vec[2]);
    printf("post rv1 %g %g %g %g %g %g\n",r1_vec_p[0],r1_vec_p[1],r1_vec_p[2],v1_vec_p[0],v1_vec_p[1],v1_vec_p[2]);
    printf("post rv2 %g %g %g %g %g %g\n",r2_vec_p[0],r2_vec_p[1],r2_vec_p[2],v2_vec_p[0],v2_vec_p[1],v2_vec_p[2]);
    printf("post rvcm1 %g %g %g %g %g %g\n",r1_CM_vec_p[0],r1_CM_vec_p[1],r1_CM_vec_p[2],v1_CM_vec_p[0],v1_CM_vec_p[1],v1_CM_vec_p[2]);
#endif

    double e1_vec_p[3],h1_vec_p[3];
    double e2_vec_p[3],h2_vec_p[3];
        
    compute_h_vector(m1_p*m2_p/(m1_p+m2_p),             r1_vec_p, v1_vec_p, h1_vec_p);
    compute_h_vector((m1_p+m2_p)*m3_p/(m1_p+m2_p+m3_p), r2_vec_p, v2_vec_p, h2_vec_p);

    compute_eccentricity_vector(m1_p+m2_p,      r1_vec_p, v1_vec_p, e1_vec_p);
    compute_eccentricity_vector(m1_p+m2_p+m3_p, r2_vec_p, v2_vec_p, e2_vec_p);

//    printf("test e 1 %g %g %g\n",e1_vec[0],e1_vec[1],e1_vec[2]);
//    printf("test e 2 %g %g %g\n",e1_vec_p[0],e1_vec_p[1],e1_vec_p[2]);
//    printf("test h 1 %g %g %g\n",h1_vec[0],h1_vec[1],h1_vec[2]);
//    printf("test h 2 %g %g %g\n",h1_vec_p[0],h1_vec_p[1],h1_vec_p[2]);

    *V1_prime_x = V1_vec_p[0];
    *V1_prime_y = V1_vec_p[1];
    *V1_prime_z = V1_vec_p[2];

    *V2_prime_x = V2_vec_p[0];
    *V2_prime_y = V2_vec_p[1];
    *V2_prime_z = V2_vec_p[2];

    *V3_prime_x = V3_vec_p[0];
    *V3_prime_y = V3_vec_p[1];
    *V3_prime_z = V3_vec_p[2];

    *e1_vec_x_p = e1_vec_p[0];
    *e1_vec_y_p = e1_vec_p[1];
    *e1_vec_z_p = e1_vec_p[2];
    *e2_vec_x_p = e2_vec_p[0];
    *e2_vec_y_p = e2_vec_p[1];
    *e2_vec_z_p = e2_vec_p[2];

    *h1_vec_x_p = h1_vec_p[0];
    *h1_vec_y_p = h1_vec_p[1];
    *h1_vec_z_p = h1_vec_p[2];
    *h2_vec_x_p = h2_vec_p[0];
    *h2_vec_y_p = h2_vec_p[1];
    *h2_vec_z_p = h2_vec_p[2];

    return 0;
}
int compute_eccentricity_vector(double total_mass, double r[3], double v[3], double e_vec[3])
{
    double v_dot_v = dot3(v,v);
    double r_dot_v = dot3(r,v);
    double r_norm = norm3(r);
    for (int i=0; i<3; i++)
    {
        e_vec[i] = (r[i]*v_dot_v - v[i]*r_dot_v)/(CONST_G*total_mass) - r[i]/r_norm;
    }
    return 0;
}
int compute_h_vector(double mu, double r[3], double v[3], double h_vec[3])
{
    cross3(r,v,h_vec);
    for (int i=0; i<3; i++)
    {
        h_vec[i] *= mu;
    }
    return 0;
}

int compute_orbital_vectors_from_orbital_elements(double m1, double m2, double m3,
    double a1, double a2, double e1, double e2,
    double INCL_rel,
    double AP1, double AP2, double LAN1, double LAN2,
    double *e1_vec_x, double *e1_vec_y, double *e1_vec_z,
    double *e2_vec_x, double *e2_vec_y, double *e2_vec_z,
    double *h1_vec_x, double *h1_vec_y, double *h1_vec_z,
    double *h2_vec_x, double *h2_vec_y, double *h2_vec_z)
{
    double tiny_double = 1.0e-5;
    if (INCL_rel < tiny_double)
    {
        INCL_rel = tiny_double;
    }

    if (INCL_rel > M_PI-tiny_double)
    {
        INCL_rel = M_PI-tiny_double;
    }

    double cos_INCL_rel = cos(INCL_rel);
    double cos_AP1 = cos(AP1);
    double cos_AP2 = cos(AP2);
    double sin_AP1 = sin(AP1);
    double sin_AP2 = sin(AP2);
    double cos_LAN1 = cos(LAN1);
//    double cos_LAN2 = cos(LAN2);
    double sin_LAN1 = sin(LAN1);
//    double sin_LAN2 = sin(LAN2);
    double sin_LAN2 = -sin_LAN1;
    double cos_LAN2 = -cos_LAN1;
    
    double L1 = m1*m2*sqrt(CONST_G*a1/(m1+m2));
    double L2 = (m1+m2)*m3*sqrt(CONST_G*a2/(m1+m2+m3));
    double G1 = L1*sqrt(1.0 - e1*e1);
    double G2 = L2*sqrt(1.0 - e2*e2);
    double Gtot = sqrt(G1*G1 + G2*G2 + 2.0*G1*G2*cos_INCL_rel);
    double cos_INCL1 = (Gtot*Gtot + G1*G1 - G2*G2)/(2.0*Gtot*G1);
    double cos_INCL2 = (Gtot*Gtot - G1*G1 + G2*G2)/(2.0*Gtot*G2);
    double INCL1 = acos(cos_INCL1);
    double INCL2 = acos(cos_INCL2);
    double sin_INCL1 = sin(INCL1);
    double sin_INCL2 = sin(INCL2);
//    sin_INCL1 = sqrt(1.0 - cos_INCL1*cos_INCL1);
//    sin_INCL2 = sqrt(1.0 - cos_INCL2*cos_INCL2);
    double h1 = G1;
    double h2 = G2;

    *e1_vec_x = e1*(cos_LAN1*cos_AP1 - sin_LAN1*sin_AP1*cos_INCL1);
    *e1_vec_y = e1*(sin_LAN1*cos_AP1 + cos_LAN1*sin_AP1*cos_INCL1);
    *e1_vec_z = e1*(sin_AP1*sin_INCL1);
    
    *h1_vec_x = h1*(sin_LAN1*sin_INCL1);
    *h1_vec_y = h1*(-cos_LAN1*sin_INCL1);
    *h1_vec_z = h1*(cos_INCL1);

    *e2_vec_x = e2*(cos_LAN2*cos_AP2 - sin_LAN2*sin_AP2*cos_INCL2);
    *e2_vec_y = e2*(sin_LAN2*cos_AP2 + cos_LAN2*sin_AP2*cos_INCL2);
    *e2_vec_z = e2*(sin_AP2*sin_INCL2);
    
    *h2_vec_x = h2*(sin_LAN2*sin_INCL2);
    *h2_vec_y = h2*(-cos_LAN2*sin_INCL2);
    *h2_vec_z = h2*(cos_INCL2);

//    double h1_vec[3] = {*h1_vec_x,*h1_vec_y,*h1_vec_z};
//    double h2_vec[3] = {*h2_vec_x,*h2_vec_y,*h2_vec_z};
//    double te = dot3(h1_vec,h2_vec)/(h1*h2);
//    printf("INCL_rel_1 %g\n",acos(te));
    return 0;
}

int compute_orbital_elements_from_orbital_vectors
(
    double m1, double m2, double m3, 
    double e1_vec_x,double e1_vec_y,double e1_vec_z,
    double e2_vec_x,double e2_vec_y,double e2_vec_z,
    double h1_vec_x,double h1_vec_y,double h1_vec_z,
    double h2_vec_x,double h2_vec_y,double h2_vec_z,
    double *a1, double *a2, double *e1, double *e2,
    double *INCL_rel,
    double *AP1, double *AP2, double *LAN1, double *LAN2)
{
    double e1_vec[3] = {e1_vec_x,e1_vec_y,e1_vec_z};
    double e2_vec[3] = {e2_vec_x,e2_vec_y,e2_vec_z};
    double e1_sq = norm3_squared(e1_vec);
    double e2_sq = norm3_squared(e2_vec);    
    *e1 = sqrt(e1_sq);
    *e2 = sqrt(e2_sq);

    double h1_vec[3] = {h1_vec_x,h1_vec_y,h1_vec_z};
    double h2_vec[3] = {h2_vec_x,h2_vec_y,h2_vec_z};
    double h1_sq = norm3_squared(h1_vec);
    double h2_sq = norm3_squared(h2_vec);
    
    *a1 = h1_sq*(m1+m2)/( CONST_G*m1*m1*m2*m2*(1.0 - e1_sq) );
    *a2 = h2_sq*(m1+m2+m3)/( CONST_G*(m1+m2)*(m1+m2)*m3*m3*(1.0 - e2_sq) );
    
    double h1 = sqrt(h1_sq);
    double h2 = sqrt(h2_sq);
    
//    double X_vec[3] = {1.0,0.0,0.0};
//    double Y_vec[3] = {0.0,1.0,0.0};
//    double Z_vec[3] = {0.0,0.0,1.0};

//    double cos_INCL1 = dot3(h1_vec,z_vec)/h1;
//    double cos_INCL2 = dot3(h2_vec,z_vec)/h2;
    double cos_INCL_rel = dot3(h1_vec,h2_vec)/(h1*h2);
    *INCL_rel = acos(cos_INCL_rel);

    double htot_vec[3];
    double X_vec[3],Y_vec[3],Z_vec[3];
    for (int i=0; i<3; i++)
    {
        htot_vec[i] = h1_vec[i] + h2_vec[i];
    }
    double htot_vec_norm = norm3(htot_vec);
    for (int i=0; i<3; i++)
    {
        Z_vec[i] = htot_vec[i]/htot_vec_norm;
    }
    double f = 1.0/sqrt( Z_vec[0]*Z_vec[0] + Z_vec[2]*Z_vec[2] );
    X_vec[0] = Z_vec[2]*f;
    X_vec[1] = 0.0;
    X_vec[2] = -Z_vec[0]*f;
    cross3(Z_vec,X_vec,Y_vec);
//    printf("X %g %g %g \n",X_vec[0],X_vec[1],X_vec[2]);
//    printf("Y %g %g %g \n",Y_vec[0],Y_vec[1],Y_vec[2]);
//    printf("Z %g %g %g \n",Z_vec[0],Z_vec[1],Z_vec[2]);
    
    double LAN1_vec[3],LAN2_vec[3];
    double LAN1_vec_unit[3],LAN2_vec_unit[3];
    cross3(Z_vec,h1_vec,LAN1_vec);

    double LAN1_vec_norm = norm3(LAN1_vec);
    double e1_vec_unit[3], e2_vec_unit[3];
    double h1_vec_unit[3], h2_vec_unit[3];
    for (int i=0; i<3; i++)
    {
        LAN1_vec_unit[i] = LAN1_vec[i]/LAN1_vec_norm;
        LAN2_vec_unit[i] = -LAN1_vec_unit[i];
        e1_vec_unit[i] = e1_vec[i]/(*e1);
        e2_vec_unit[i] = e2_vec[i]/(*e2);
        h1_vec_unit[i] = h1_vec[i]/h1;
        h2_vec_unit[i] = h2_vec[i]/h2;
    }

    double sin_LAN1 = dot3(LAN1_vec,Y_vec);
    double sin_LAN2 = -sin_LAN1;
//    double sin_LAN2 = dot3(LAN2_vec,y_vec);    
    double cos_LAN1 = dot3(LAN1_vec,X_vec);
//    double cos_LAN2 = dot3(LAN2_vec,x_vec);    
    double cos_LAN2 = -cos_LAN1;

    double e1_cross_h1[3], e2_cross_h2[3];
    cross3(e1_vec_unit,h1_vec_unit,e1_cross_h1);
    cross3(e2_vec_unit,h2_vec_unit,e2_cross_h2);
    double sin_AP1 = dot3(LAN1_vec_unit,e1_cross_h1);
    double sin_AP2 = dot3(LAN2_vec_unit,e2_cross_h2);
    double cos_AP1 = dot3(LAN1_vec_unit,e1_vec_unit);
    double cos_AP2 = dot3(LAN2_vec_unit,e2_vec_unit);

    *LAN1 = atan2(sin_LAN1,cos_LAN1);
    *LAN2 = atan2(sin_LAN2,cos_LAN2);
    *AP1 = atan2(sin_AP1,cos_AP1);
    *AP2 = atan2(sin_AP2,cos_AP2);

    return 0;
}


void cross3(double a[3], double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}
double norm3(double v[3])
{
    double result = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return result;
}
double norm3_squared(double v[3])
{
    double result = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return result;
}
double dot3(double a[3], double b[3])
{
    double result = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    return result;
}
