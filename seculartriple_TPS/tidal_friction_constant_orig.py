#not used
"""
Routine that calculates the ratio k/T for tidal friction, where k is the apsidal motion constant and T the tidal friction time scale as they appear in Hut (1981; 1981A&A....99..126H).
The entire routine has been copied directly from c-code that is part of binary_c and converted to Python.
Needs to be checked/tested. Note that Izzard's "E2-prescription" is not yet implemented.

Adrian Hamers 16-06-2014
"""

#define USE_RADIATIVE_DAMPING (((stardata->star[star_number].stellar_type==MAIN_SEQUENCE)&&(MORE_OR_EQUAL(stardata->star[star_number].mass,1.25)))||(stardata->star[star_number].stellar_type==CHeB)||(stardata->star[star_number].stellar_type==HeMS))

#define USE_CONVECTIVE_DAMPING (stardata->star[star_number].stellar_type<HeWD)

from amuse.units import units,constants
import numpy

MAIN_SEQUENCE=1|units.stellar_type
CHeB=4|units.stellar_type
HeMS=7|units.stellar_type
HeWD=10|units.stellar_type

def check_for_radiative_damping(stellar_type,mass):
    if (stellar_type==MAIN_SEQUENCE and mass.value_in(units.MSun) >= 1.25):
        return True
    elif (stellar_type == CHeB or stellar_type == HeMS):
        return True
    else:
        return False
   
def check_for_convective_damping(stellar_type):    
    if stellar_type < HeWD:
        return True
    else:
        return False

def tidal_friction_constant(stellar_type,mass,companion_mass,semimajor_axis,radius,envelope_mass,envelope_radius,luminosity,spin_angular_frequency,gyration_radius_star1):
    
    USE_RADIATIVE_DAMPING = check_for_radiative_damping(stellar_type,mass)
    USE_CONVECTIVE_DAMPING = check_for_convective_damping(stellar_type)
    
    if USE_RADIATIVE_DAMPING == True: ### radiative damping ###
        E2 = 1.592e-09*pow(mass.value_in(units.MSun),2.84) ### Hurley prescription; Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329. ###
        
        ### Izzard's prescription not yet implemented ###
        """
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
        """   
        
        k_div_T_tides = E2*pow(1.0 + companion_mass/mass,5.0/6.0)*radius*numpy.sqrt(constants.G*mass/(semimajor_axis**5))
        return k_div_T_tides
        
#        kTradiative_damping = 1.9782e+04*sqrt((mass.value_in(units.MSun)*(radius.value_in(units.RSun))**2))/((semimajor_axis.value_in(units.AU))**5))*E2*pow(1.0+companion_mass/mass,5.0/6.0)
#		kTradiative_damping = kTradiative_damping/YEAR_LENGTH_IN_SECONDS;	/* This converts (k/T) to units of s^-1 */
             
    elif USE_CONVECTIVE_DAMPING == True: ### convective damping ###
        P_orb = 2.0*numpy.pi*numpy.sqrt((semimajor_axis**3)/(constants.G*(mass + companion_mass)))
        if spin_angular_frequency.value_in(1.0/units.s) == 0.0:
            P_tid = P_orb
        else:
            P_spin = 2.0*numpy.pi/spin_angular_frequency
            P_tid_s = 1.0/( 1e-10 + numpy.fabs( 1.0/(P_orb.value_in(units.s)) - 1.0/(P_spin.value_in(units.s)) ) )
            P_tid = P_tid_s | units.s

        tau_convective = pow( (envelope_mass*envelope_radius*(radius - (1.0/2.0)*envelope_radius))/(3.0*luminosity), 1.0/3.0)
	print 'tau',envelope_mass,envelope_mass*envelope_radius*(radius - (1.0/2.0)*envelope_radius)/(3.0*luminosity)

	#print 'tau convective',tau_convective
        f_convective = (P_tid/(2.0*tau_convective))**2

        f_convective = numpy.amin([1.0,f_convective])
        
        k_div_T_tides = (2.0/21.0)*(f_convective/tau_convective)*(envelope_mass/mass)
        print 'convective damping', k_div_T_tides
        return k_div_T_tides

    else: ### degenerate damping -- 1984MNRAS.207..433C ###
        tau_degenerate = 1.3e7 | units.yr 
        k_div_T_tides = (1.0/(3.0*tau_degenerate))*gyration_radius_star1**2*pow(luminosity.value_in(units.LSun)/mass.value_in(units.MSun),5.0/7.0)
        print 'degenerate damping', k_div_T_tides
        return k_div_T_tides
