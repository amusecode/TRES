/*
Routine that calculates the ratio k/T for tidal friction, where k is the apsidal motion constant and T the tidal friction time scale as they appear in Hut (1981; 1981A&A....99..126H).
The entire routine has been copied directly from c-code that is part of binary_c -- need to talk to Rob Izzard about this before making the code public

!!! Needs to be thoroughly checked/tested !!!

Note that Izzard's "E2-prescription" is not yet implemented.

Adrian Hamers 29-09-2014
*/

//#define USE_RADIATIVE_DAMPING (((stardata->star[star_number].stellar_type==MAIN_SEQUENCE)&&(MORE_OR_EQUAL(stardata->star[star_number].mass,1.25)))||(stardata->star[star_number].stellar_type==CHeB)||(stardata->star[star_number].stellar_type==HeMS))

//#define USE_CONVECTIVE_DAMPING (stardata->star[star_number].stellar_type<HeWD)

#include "main_code.h"

double const MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN = 1.2; // in future: make user-adjustable
int const MAIN_SEQUENCE = 1;
int const CHeB = 4;
int const HeMS = 7;
int const HeWD = 10;

int const NS = 13; 
int const BH = 14;
int const PREMS = 17;
int const PLANET = 18;
int const BD = 19;



double set_crude_gyration_radii_based_on_stellar_structure(int stellar_type, double mass)
{ 
    //printf("tides in c");
    double gyration_radius;
    
    bool USE_RADIATIVE_DAMPING = check_for_radiative_damping(stellar_type,mass);
    bool USE_CONVECTIVE_DAMPING = check_for_convective_damping(stellar_type);
    
    if ((USE_RADIATIVE_DAMPING==TRUE) || (USE_CONVECTIVE_DAMPING == TRUE))
    {
        gyration_radius = 0.1;
    }
    else
    {
        gyration_radius = 0.21;
    }
    return gyration_radius;
}

bool check_for_radiative_damping(int stellar_type, double mass)
{
    if ((stellar_type == MAIN_SEQUENCE) && (mass/CONST_MSUN >= MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN))
    {
        return TRUE;
    }
    else if ((stellar_type == CHeB) || (stellar_type == HeMS))
    {
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}   
bool check_for_convective_damping(int stellar_type)
{
    if ((stellar_type < HeWD) || (stellar_type == PREMS))
    {
        return TRUE;
    }
    else
    {
        return FALSE;
    }
}

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
)
{
    bool USE_RADIATIVE_DAMPING = check_for_radiative_damping(stellar_type,mass);
    bool USE_CONVECTIVE_DAMPING = check_for_convective_damping(stellar_type);
    double k_div_T_tides;
    
        
    if (USE_RADIATIVE_DAMPING == TRUE) // radiative damping
    {
        double E2 = 1.592e-09*pow(mass/CONST_MSUN,2.84); // Hurley prescription; Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329
        
#ifdef IGNORE
        // Izzard's prescription not yet implemented
		else if(stardata->preferences->E2_prescription==E2_IZZARD)
		{
			if(stardata->star[star_number].stellar_type<HERTZSPRUNG_GAP)
			{
				double fburn=1.0;
				double x=stardata->star[star_number].aj/stardata->star[star_number].tms;
				if(x<fburn)
				{
					/* log mass and Z */
					double logm=log10(stardata->star[star_number].mass);
					double logz=log10(stardata->common.metallicity);
 
					/* fits for am and E20 */
					double am = 0.15*sin(3.2*logm)+0.31*logm;
					double E20 = -1.23250e+01+1.04550e+01*logm-4.15420e-01*logz-7.18650e+00*logm*logm+1.97940e+00*logm*logm*logm;
			  		E20=pow(10.0,E20);

					/* calc log(E2/E20) */
					E2 = -pow(x+am,4.0)*pow(MAX(1,x/0.95),30.0);
			  
					/* hence E2 */
					E2 = E20 * pow(10.0,E2);
			  
					/* verbosity */
					/*
					if(star->starnum==1)
					{
						printf("E2 kw=%d I=%g (fburn=%g x=%g E20=%g) H=%g\n",
						star->stellar_type,
						E2,
						fburn,x,
						E20,E2_Hurley);
					}
					*/
				}
				else
				{
					/* no conv core */
					E2=0.0;
				}
			}
			else
			{
				E2=0.0;
			}
		}
#endif
        
        k_div_T_tides = E2*pow(1.0 + companion_mass/mass,5.0/6.0)*radius*sqrt(CONST_G*mass/(pow(semimajor_axis,5.0)));
//        printf("rad damping %g %g \n",mass, k_div_T_tides);
        return k_div_T_tides*1.;
        
//        kTradiative_damping = 1.9782e+04*sqrt((mass.value_in(units.MSun)*(radius.value_in(units.RSun))**2))/((semimajor_axis.value_in(units.AU))**5))*E2*pow(1.0+companion_mass/mass,5.0/6.0)
//		kTradiative_damping = kTradiative_damping/YEAR_LENGTH_IN_SECONDS;	/* This converts (k/T) to units of s^-1 */
    }
    else if (USE_CONVECTIVE_DAMPING == TRUE) // convective damping
    {
        double P_orb = 2.0*M_PI*sqrt((semimajor_axis*semimajor_axis*semimajor_axis)/(CONST_G*(mass + companion_mass)));
//        printf("a %g\n",semimajor_axis);
        double P_spin,P_tid;
        
        if (spin_angular_frequency == 0.0)
        {
            P_tid = P_orb;
        }
        else
        {
            P_spin = 2.0*M_PI/spin_angular_frequency;
            P_tid = 1.0/( 1e-10 + fabs( 1.0/P_orb - 1.0/P_spin) );

        double tau_convective = pow( (convective_envelope_mass*convective_envelope_radius*(radius - (1.0/2.0)*convective_envelope_radius))/(3.0*luminosity), 1.0/3.0);
//	print 'tau',envelope_mass,envelope_mass*envelope_radius*(radius - (1.0/2.0)*envelope_radius)/(3.0*luminosity)

        double f_convective = pow(P_tid/(2.0*tau_convective),2.0);
        f_convective = min(1.0,f_convective);

        k_div_T_tides = (2.0/21.0)*(f_convective/tau_convective)*(convective_envelope_mass/mass);
        
        if ((convective_envelope_mass <= 0.0) || (convective_envelope_radius <= 0.0))
        {
            k_div_T_tides = 0.0;
        }
        
//        printf("test par conv %g %g %g %g %g \n",mass,radius,convective_envelope_mass,convective_envelope_radius,spin_angular_frequency);
//        printf("test conv %g %g %g %g %g \n",P_orb,tau_convective,P_tid,P_spin,f_convective);
//        printf("convective damping %g %g %g %g \n",mass, k_div_T_tides, f_convective, tau_convective);
//        return 3.0e+5;//t_viscous = 5yr based on Fabrycky & Tremaine for 1 Msun star
//        printf("conv damping %g %g \n",mass, k_div_T_tides);

        return k_div_T_tides;//*100;//*10000;
        }
    }
    
    else if (stellar_type == NS or stellar_type == BH) //no tides for NS or BH
    {
        k_div_T_tides = 0;   
//        printf("ns/bh tidess %g %i \n", k_div_T_tides, stellar_type);   
        return k_div_T_tides;      
    }
    else if (stellar_type == PLANET or stellar_type == BD)  //based on Fabrycky & Tremaine 2007, appendix
    { 
        double T_viscous = 0.001/1e6; // in Myr
//        printf("planet/bd tidess %g %g %i \n", k, T_viscous, stellar_type);   
        return amc/T_viscous;      
    }
    else // degenerate damping -- 1984MNRAS.207..433C
    {

//        double seconds_in_year = 365.25*24.0*3600.0;
//        double tau_degenerate = 1.3e7*seconds_in_year;
        double tau_degenerate = 1.3e7/1e6;// in Myr
        k_div_T_tides = (1.0/(3.0*tau_degenerate))*gyration_radius*gyration_radius*pow((luminosity/CONST_L_SUN)/(mass/CONST_MSUN),5.0/7.0);
        
//        printf("deg damping %g  \n",k_div_T_tides);
//        printf("deg damping %g %g %g %g \n",k_div_T_tides, gyration_radius, luminosity/CONST_L_SUN,mass/CONST_MSUN);
        return k_div_T_tides;
    }
}
