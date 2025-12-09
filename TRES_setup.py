from amuse.community.seba.interface import SeBa
from amuse.datamodel import Particles
from amuse.units import units
from seculartriple_TPS.interface import SecularTriple

from TRES_options import max_mass, absolute_min_mass, options_mesa, REPORT_USER_WARNINGS

from interactions import *
from tidal_friction_constant import *

import numpy as np

#-------
#to initialize the triple object
def make_stars(mass_list):
    stars = Particles(len(mass_list))
    stars.is_star = True
    stars.is_donor = False
    
    if max(mass_list) > max_mass:
        print('error: masses not in allowed range')
        print(mass_list, 'should be below:', max_mass)
        print('max_mass settable in TRES_options.py')
        return stars, False

    if min(mass_list) <= absolute_min_mass:
        print('error: masses not in allowed range')
        print(mass_list, 'should be at least above:', absolute_min_mass)
        print('absolute_min_mass settable in TRES_options.py')
        print('substellar objects can be included through EXCLUDE_SSO in TRES_options.py')
        return stars, False

    if mass_list[0] < mass_list[1]:
        spare = mass_list[0]
        mass_list[0] = mass_list[1]
        mass_list[1] = spare

    stars.mass = mass_list
    stars.initial_mass = mass_list
    
    return stars, True

def make_bins(stars, semimajor_axis_list,
        eccentricity_list,
        relative_inclination_list,
        argument_of_pericenter_list,
        longitude_of_ascending_node_list):

    bins = Particles(len(semimajor_axis_list))
    bins.is_star = False
    bins.is_mt_stable = True
    bins.part_dt_mt = 1.
    bins.bin_type = bin_type['unknown'] 
    
    
    if min(semimajor_axis_list) < 0.|units.RSun:
        print('error: semimajor_axis not in allowed range', semimajor_axis_list)
        return bins, False, eccentricity_list
    if (min(eccentricity_list) < 0.) or (max(eccentricity_list) >= 1.):
        print('error: eccentricity not in allowed range', eccentricity_list)
        return bins, False, eccentricity_list
    for i in range(len(eccentricity_list)):
        if eccentricity_list[i] < minimum_eccentricity:
            eccentricity_list[i] = minimum_eccentricity
 
    if (min(relative_inclination_list) < 0.) or (max(relative_inclination_list) > np.pi):
        print('error: relative_inclination not in allowed range', relative_inclination_list)
        return bins, False, eccentricity_list
    if (min(argument_of_pericenter_list) < -1.*np.pi) or (max(argument_of_pericenter_list) > np.pi):
        print('error: argument_of_pericenter not in allowed range', argument_of_pericenter_list)
        return bins, False, eccentricity_list
    if (min(longitude_of_ascending_node_list) < -1.*np.pi) or (max(longitude_of_ascending_node_list) > np.pi):
        print('error: longitude_of_ascending_node not in allowed range', longitude_of_ascending_node_list)
        return bins, False, eccentricity_list
      
    if len(semimajor_axis_list) == 2: 
        longitude_of_ascending_node_list.append(longitude_of_ascending_node_list[0] - np.pi)
        if semimajor_axis_list[0] >= semimajor_axis_list[1]: 
            print('error input parameters, should be:')
            print('inner_semimajor_axis < outer_semimajor_axis' )
            return bins, False, eccentricity_list
    if len(semimajor_axis_list) > 2: 
        print('quadruples not implemented yet in make_bins')
        return bins, False, eccentricity_list
        
       
    bins[0].child1 = stars[0]
    bins[0].child2 = stars[1]
    if len(bins) > 1: #adjust for quads
        bins[1].child1 = stars[2]
        bins[1].child2 = bins[0]
        
    for i in range(len(bins)):
        bins[i].child1.parent = bins[i]
        bins[i].child2.parent = bins[i]
        bins[i].semimajor_axis = semimajor_axis_list[i]
        bins[i].eccentricity = eccentricity_list[i]
        bins[i].argument_of_pericenter = argument_of_pericenter_list[i]
        bins[i].longitude_of_ascending_node = longitude_of_ascending_node_list[i]

    
    bins.mass_transfer_rate = 0.0 | units.MSun/units.yr
    bins.accretion_efficiency_mass_transfer = 1.0
    bins.accretion_efficiency_wind_child1_to_child2 = 0.0
    bins.accretion_efficiency_wind_child2_to_child1 = 0.0
    # binary evolutionary settings
    bins.specific_AM_loss_mass_transfer = 2.5

    return bins, True, eccentricity_list
    
    

def make_particle_sets(
        mass_list, semimajor_axis_list, eccentricity_list,  
        relative_inclination_list, argument_of_pericenter_list, 
        longitude_of_ascending_node_list):

        stars, correct_params_mass = make_stars(mass_list)
                
        bins, correct_params_orbit, eccentricity_list = make_bins(
            stars, semimajor_axis_list, eccentricity_list,
            relative_inclination_list,
            argument_of_pericenter_list, longitude_of_ascending_node_list)

        correct_params = True
        if not correct_params_mass or not correct_params_orbit:
            correct_params = False

        return stars, bins, correct_params

#-------
#setup community codes

def setup_stellar_code(stellar_code, stars):
    stellar_code.particles.add_particles(stars)

    if stellar_code.__module__.split(".")[-2]=="mesa_r15140":
        options_mesa(stellar_code)

    return stellar_code


def setup_secular_code(triple, secular_code, stop_at_semisecular_regime):
    triple_set = triple.as_set()
    triple_time = triple_set.time
    secular_code.triples.add_particles(triple_set)
    secular_code.parameters.verbose = False
#    secular_code.parameters.verbose = True #AB Set to true when debugging the secular code

    #needed for initialisation in some circumstances
    secular_code.model_time = triple_time

    secular_code.parameters.equations_of_motion_specification = 0
    secular_code.parameters.roche_radius_specification = 0
    #0: eccentric eggleton, 1: sepinsky, 2: classical circular eggleton
    secular_code.parameters.stability_limit_specification = 0
    #for stars 0, 5-6, for exoplanets 1-4
    #0: mardling & aarseth 2001, 1:petrovich et al. 2015 simple, 2:petrovich et al. 2015
    #3: holman et al. 98 s-type, 4: holman et al. 98 p-type,
    #5: vynatheya+ 22
    #6: tory+ 22

    secular_code.parameters.ignore_tertiary = False

    secular_code.parameters.include_quadrupole_terms = True
    secular_code.parameters.include_octupole_terms = True
    secular_code.parameters.include_inner_wind_terms = True
    secular_code.parameters.include_outer_wind_terms = True
    secular_code.parameters.include_inner_RLOF_terms = True
    secular_code.parameters.include_outer_RLOF_terms = True
    secular_code.parameters.include_magnetic_braking_terms = False # not tested

    secular_code.parameters.include_inner_tidal_terms = True
    secular_code.parameters.include_outer_tidal_terms = True

    secular_code.parameters.include_tertiary_tidal_terms_circ = False
    secular_code.parameters.include_tertiary_tidal_terms = False
    if (secular_code.parameters.include_tertiary_tidal_terms_circ and secular_code.parameters.include_tertiary_tidal_terms):
        if REPORT_USER_WARNINGS:
            print('Both circular and eccentric tertiary tides are switched on. Please pick one or the other.' )
        sys.exit('Both circular and eccentric tertiary tides are switched on. Please pick one or the other.')


    secular_code.parameters.include_1PN_inner_terms = True
    secular_code.parameters.include_1PN_outer_terms = True
    secular_code.parameters.include_1PN_inner_outer_terms = False ### warning: probably broken
    secular_code.parameters.include_25PN_inner_terms = True
    secular_code.parameters.include_25PN_outer_terms = True

    secular_code.parameters.check_for_dynamical_stability = True
    secular_code.parameters.check_for_dynamical_stability_at_initialisation = True

    secular_code.parameters.check_for_semisecular_regime = stop_at_semisecular_regime
    secular_code.parameters.check_for_semisecular_regime_at_initialisation = stop_at_semisecular_regime

    secular_code.parameters.check_for_inner_collision = True
    secular_code.parameters.check_for_outer_collision = True

    secular_code.parameters.check_for_inner_RLOF = True
    secular_code.parameters.check_for_outer_RLOF = True

    secular_code.parameters.include_spin_radius_mass_coupling_terms_star1 = True
    secular_code.parameters.include_spin_radius_mass_coupling_terms_star2 = True
    secular_code.parameters.include_spin_radius_mass_coupling_terms_star3 = True

     # accuracy of secular code
#        secular_code.parameters.input_precision = 1.0e-10#1.0e-5
#        secular_code.parameters.relative_tolerance = 1.0e-10
#        secular_code.parameters.threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = 1.0e-12
    secular_code.parameters.threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero = 1.0e-7|units.Myr**-1

    secular_code.parameters.include_linear_mass_change = True #needed for Jspin conservation
    secular_code.parameters.include_linear_radius_change = True #needed for Jspin conservation

#    channel_from_secular = secular_code.triples.new_channel_to(triple_set)
#    channel_to_secular = triple_set.new_channel_to(secular_code.triples)

    return secular_code

#-------
