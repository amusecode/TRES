#not used - file should be removed 
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


MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN = 1.2 # in future: make user-adjustable

MAIN_SEQUENCE=1|units.stellar_type
CHeB=4|units.stellar_type
HeMS=7|units.stellar_type
HeWD=10|units.stellar_type

NS=13|units.stellar_type
BH=14|units.stellar_type
PREMS=17|units.stellar_type
PLANET=18|units.stellar_type
BD=19|units.stellar_type



def set_gyration_radius(stellar_type, mass):
    USE_RADIATIVE_DAMPING = check_for_radiative_damping(stellar_type,mass)
    USE_CONVECTIVE_DAMPING = check_for_convective_damping(stellar_type)
    
    if (USE_RADIATIVE_DAMPING or USE_CONVECTIVE_DAMPING):
        gyration_radius = 0.1
    else:
        gyration_radius = 0.21
    return gyration_radius



def check_for_radiative_damping(stellar_type,mass):
    if (stellar_type==MAIN_SEQUENCE and mass.value_in(units.MSun) >= MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN):
        return True
    elif (stellar_type == CHeB or stellar_type == HeMS):
        return True
    else:
        return False
   
def check_for_convective_damping(stellar_type):    
    if stellar_type < HeWD or stellar_type == PREMS:
        return True
    else:
        return False

def tidal_friction_constant(stellar_type,mass,companion_mass,semimajor_axis,radius,convective_envelope_mass,convective_envelope_radius,luminosity,spin_angular_frequency, gyration_radius, amc):
    
    USE_RADIATIVE_DAMPING = check_for_radiative_damping(stellar_type,mass)
    USE_CONVECTIVE_DAMPING = check_for_convective_damping(stellar_type)
    
    if USE_RADIATIVE_DAMPING: ### radiative damping ###
        E2 = 1.592e-09*pow(mass.value_in(units.MSun),2.84) ### Hurley prescription; Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329. ###
                
        k_div_T_tides = E2*pow(1.0 + companion_mass/mass,5.0/6.0)*radius*numpy.sqrt(constants.G*mass/(semimajor_axis**5))
#        print('radiative damping', mass, k_div_T_tides)
        return k_div_T_tides
        
#        kTradiative_damping = 1.9782e+04*sqrt((mass.value_in(units.MSun)*(radius.value_in(units.RSun))**2))/((semimajor_axis.value_in(units.AU))**5))*E2*pow(1.0+companion_mass/mass,5.0/6.0)
#		kTradiative_damping = kTradiative_damping/YEAR_LENGTH_IN_SECONDS;	/* This converts (k/T) to units of s^-1 */
             
    elif USE_CONVECTIVE_DAMPING: ### convective damping ###
        P_orb = 2.0*numpy.pi*numpy.sqrt((semimajor_axis**3)/(constants.G*(mass + companion_mass)))
        if spin_angular_frequency.value_in(1.0/units.s) == 0.0:
            P_tid = P_orb
        else:
            P_spin = 2.0*numpy.pi/spin_angular_frequency
            P_tid_s = 1.0/( 1e-10 + numpy.fabs( 1.0/(P_orb.value_in(units.s)) - 1.0/(P_spin.value_in(units.s)) ) )
            P_tid = P_tid_s | units.s

        tau_convective = pow( (convective_envelope_mass*convective_envelope_radius*(radius - (1.0/2.0)*convective_envelope_radius))/(3.0*luminosity), 1.0/3.0)
#	print 'tau',envelope_mass,envelope_mass*envelope_radius*(radius - (1.0/2.0)*envelope_radius)/(3.0*luminosity)

	#print 'tau convective',tau_convective
        f_convective = (P_tid/(2.0*tau_convective))**2

        f_convective = numpy.amin([1.0,f_convective])
        
        k_div_T_tides = (2.0/21.0)*(f_convective/tau_convective)*(convective_envelope_mass/mass)
#        print('convective damping', mass, k_div_T_tides)
        return k_div_T_tides

    elif stellar_type==NS or stellar_type==BH:
        ### no tides for NS or BH
        k_div_T_tides = 0      
#        print('ns/bh tides', k_div_T_tides, stellar_type)
        return k_div_T_tides      
    elif stellar_type==PLANET or stellar_type==BD: #based on Fabrycky & Tremaine 2007, appendix
        T_viscous = 0.001|units.yr
        return amc/T_viscous
    else: ### degenerate damping -- 1984MNRAS.207..433C ###
        tau_degenerate = 1.3e7 | units.yr 
        k_div_T_tides = (1.0/(3.0*tau_degenerate))*gyration_radius**2*pow(luminosity.value_in(units.LSun)/mass.value_in(units.MSun),5.0/7.0)
#        print('degenerate damping', k_div_T_tides)
#        print('degenerate damping ', k_div_T_tides, gyration_radius_star1, luminosity.value_in(units.LSun),mass.value_in(units.MSun))
        return k_div_T_tides
