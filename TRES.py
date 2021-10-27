# to do 
# min teken in mean anomaly

## Triple:      Triple evolution
##              computes the evolution of a given triple
##              given any initial conditions (M, m, l, A, a, E, e, i, G, g, O, o, T, z).
 
from amuse.community.seba.interface import SeBa
from interactions import *
from tidal_friction_constant import *

import os, sys
import time
from seculartriple_TPS.interface import SecularTriple
from amuse.units import units, constants
from amuse.datamodel import Particles
from amuse.support.console import set_printing_strategy
from amuse.io import write_set_to_file
from amuse.units import quantities
from scipy.stats import maxwell
from scipy import optimize
from math import sqrt, isnan
import numpy as np

REPORT_USER_WARNINGS = True

REPORT_DEBUG = False
REPORT_DT = False 
REPORT_SN_EVOLUTION = False
REPORT_TRIPLE_EVOLUTION = False 

no_stellar_evolution = False

#constants
time_step_factor_stable_mt = 0.01 #1% mass loss during mass transfer
# lowering this to 0.005 makes the code twice as slow
time_step_factor_ecc = 0.01
#Rl_fraction = 0.8
# 0.01 -> error in the semi-major axis of about 0.5%
maximum_wind_mass_loss_factor = 0.01 
error_dm = 0.05
#maximum_radius_change_factor = 0.005
error_dr = 0.05 #0.01
minimum_time_step = 1.e-9 |units.Myr
# min_mass = 0.08 |units.MSun # for stars
# absolute_min_mass = 0.008|units.MSun # AMUSE can't handle planets -> also for secondaries and tertiaries
max_mass = 100 |units.MSun
maximum_time_step_factor = 100.
maximum_time_step_factor_after_stable_mt = 5. 
time_step_factor_find_RLOF = 0.5
#Rl_fraction = 0.9#1.0-10.*error_dr # ratio or star radius over Roche lobe at which time step is decreased
                              # radius grows maximally by error_dr
time_step_factor_kozai = 0.025 # 0.2*0.1, 0.2-> for error in kozai timescale, 0.1-> 10 steps per cycle
kozai_type_factor = 10.

kanonical_neutron_star_mass = 1.4|units.MSun
fall_back_mass = 41 |units.MSun

stellar_types_SN_remnants = [13,14,15]|units.stellar_type # remnant types created through a supernova
stellar_types_remnants = [7,8,9,10,11,12,13,14,15]|units.stellar_type
stellar_types_dr = [2,4,7,8,9,10,11,12,13,14,15]|units.stellar_type #stars which go through a instantaneous radius change at formation; hertzsprung gap stars (small envelope perturbation) + horizontal branch stars + remnants


class Triple_Class:
    #-------
    #setup stellar system
    def __init__(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, 
            metallicity, tend, number, maximum_radius_change_factor, 
            stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer,
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime,  
            stop_at_SN, SN_kick_distr, stop_at_CPU_time, max_CPU_time,
            file_name, file_type, dir_plots):
        
        self.set_stopping_conditions(stop_at_mass_transfer, stop_at_init_mass_transfer,stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer,
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime,  stop_at_SN, stop_at_CPU_time)
            
        if inner_primary_mass < inner_secondary_mass:
            spare = inner_primary_mass
            inner_primary_mass = inner_secondary_mass
            inner_secondary_mass = spare     

        correct_params, inner_eccentricity, outer_eccentricity = self.test_initial_parameters(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis, inner_eccentricity, outer_eccentricity,
            relative_inclination, inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node)   
        if correct_params == False:
            self.triple = Particles(1)
            self.triple.correct_params = correct_params 
            return
            
        outer_longitude_of_ascending_node = inner_longitude_of_ascending_node - np.pi
                        
        stars = self.make_stars(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis)
        bins = self.make_bins(stars, inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node)

            
        self.instantaneous_evolution = False # no secular evolution        
        self.maximum_radius_change_factor = maximum_radius_change_factor
        self.fixed_timestep = -1|units.Myr

        self.tend = tend 
        self.previous_time = 0.0|units.yr
        self.previous_dt = 0.0|units.yr
        self.file_name = file_name
        self.file_type = file_type
        self.SN_kick_distr = SN_kick_distr
        self.max_CPU_time = max_CPU_time

        self.triple = bins[1]
        self.triple.time = 0.0|units.yr
        self.triple.relative_inclination = relative_inclination 
        self.triple.is_star = False #maybe not necessary?
        self.triple.dynamical_instability = False 
        self.triple.number = number 
        self.triple.error_flag_secular = 0
        self.triple.CPU_time = 0.0
            
        self.setup_stellar_code(metallicity, stars)
        self.setup_secular_code(self.triple.as_set())      
        self.initial_angular_frequency() 
                        
        self.triple.dynamical_instability_at_initialisation = False
        self.triple.semisecular_regime_at_initialisation = False
        self.triple.mass_transfer_at_initialisation = False
        self.triple.correct_params = correct_params 
           
        self.secular_code.check_for_dynamical_stability()
        if stop_at_dynamical_instability == True and self.secular_code.triples[0].dynamical_instability == True:
            self.triple.dynamical_instability_at_initialisation = True
            self.triple.dynamical_instability = True
            return 

        self.secular_code.check_for_semisecular_regime()
        if stop_at_semisecular_regime == True and self.secular_code.triples[0].semisecular_regime == True:
            self.triple.semisecular_regime_at_initialisation = True
            self.triple.semisecular_regime = True        
            return
        
#        if stop_at_semisecular_regime == True:                   
#            self.secular_code.check_for_semisecular_regime()
#            if self.secular_code.triples[0].semisecular_regime == True:
#                self.triple.semisecular_regime_at_initialisation = True
#                self.triple.semisecular_regime = True
#            return             

        self.check_RLOF() 
        if self.has_tertiary_donor() and (self.stop_at_outer_mass_transfer or self.stop_at_mass_transfer or self.stop_at_init_mass_transfer): 
            self.triple.mass_transfer_at_initialisation = True
            self.triple.bin_type = bin_type['rlof']
            return
        if self.has_donor() and (self.stop_at_mass_transfer or self.stop_at_init_mass_transfer):
            self.triple.mass_transfer_at_initialisation = True
            if self.is_binary(self.triple.child2):
                self.triple.child2.bin_type = bin_type['rlof'] 
            elif self.is_binary(self.triple.child1):
                self.triple.child1.bin_type = bin_type['rlof']    
            else:
                print('currently not implemented')
                exit(-1)    
            return
        
        self.triple.kozai_type = self.get_kozai_type()
        self.update_stellar_parameters() 
        self.update_time_derivative_of_radius()
        self.update_previous_stellar_parameters()
        

    def set_stopping_conditions(self, stop_at_mass_transfer,stop_at_init_mass_transfer,stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer,
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime, stop_at_SN, stop_at_CPU_time):

        if stop_at_disintegrated == False:
            print('stop_at_disintegrated = False not possible yet. After the disintegration of the triple, further evolution can be done with SeBa directly. ') 
            exit(1)
        if stop_at_outer_mass_transfer == False:
            print('stop_at_outer_mass_transfer = False not possible yet. Methodology is as of yet non-existent.' )
            exit(1)
        if stop_at_outer_collision == False:
            print('stop_at_outer_collision = False not possible. Non-hierarchical triples can not be simulated using the secular equations as used in TRES. Further evolution should be done by other means, e.g. one of the N-body codes implemented in AMUSE.' )
            exit(1)
        if stop_at_dynamical_instability == False:
            print('stop_at_dynamical_instability = False not possible. Unstable triples can not be simulated using the secular equations as used in TRES. Further evolution should be done by other means, e.g. one of the N-body codes implemented in AMUSE.') 
            exit(1)

                            
        self.stop_at_mass_transfer = stop_at_mass_transfer            
        self.stop_at_init_mass_transfer = stop_at_init_mass_transfer
        self.stop_at_outer_mass_transfer = stop_at_outer_mass_transfer            

        self.stop_at_stable_mass_transfer =  stop_at_stable_mass_transfer
        self.stop_at_eccentric_stable_mass_transfer = stop_at_eccentric_stable_mass_transfer
        self.stop_at_unstable_mass_transfer = stop_at_unstable_mass_transfer
        self.stop_at_eccentric_unstable_mass_transfer = stop_at_eccentric_unstable_mass_transfer

        self.stop_at_merger = stop_at_merger            
        self.stop_at_disintegrated = stop_at_disintegrated            
        self.stop_at_inner_collision = stop_at_inner_collision            
        self.stop_at_outer_collision = stop_at_outer_collision            

        self.stop_at_dynamical_instability = stop_at_dynamical_instability            
        self.stop_at_semisecular_regime = stop_at_semisecular_regime
        self.stop_at_SN = stop_at_SN
        self.stop_at_CPU_time = stop_at_CPU_time
    
    def make_stars(self, inner_primary_mass, inner_secondary_mass, outer_mass, inner_semimajor_axis, outer_semimajor_axis):
        stars = Particles(3)
        stars.is_star = True
        stars.is_donor = False

        stars[0].mass = inner_primary_mass
        stars[1].mass = inner_secondary_mass
        stars[2].mass = outer_mass

        stars[0].initial_mass = inner_primary_mass
        stars[1].initial_mass = inner_secondary_mass
        stars[2].initial_mass = outer_mass

        return stars 
         
    def make_bins(self, stars, inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node):     

        bins = Particles(2)
        bins.is_star = False
        bins.is_mt_stable = True
        bins.part_dt_mt = 1.
        bins.bin_type = bin_type['unknown'] #Unknown

        bins[0].child1 = stars[0]
        bins[0].child2 = stars[1]
        bins[0].child1.parent = bins[0]
        bins[0].child2.parent = bins[0]

        bins[0].semimajor_axis = inner_semimajor_axis
        bins[0].eccentricity = inner_eccentricity
        bins[0].argument_of_pericenter = inner_argument_of_pericenter
        bins[0].longitude_of_ascending_node = inner_longitude_of_ascending_node
        
        bins[0].mass_transfer_rate = 0.0 | units.MSun/units.yr
        bins[0].accretion_efficiency_mass_transfer = 1.0
        bins[0].accretion_efficiency_wind_child1_to_child2 = 0.0
        bins[0].accretion_efficiency_wind_child2_to_child1 = 0.0

        bins[1].child1 = stars[2]
        bins[1].child2 = bins[0]
        bins[1].child1.parent = bins[1]
        bins[1].child2.parent = bins[1]
        
        bins[1].semimajor_axis = outer_semimajor_axis
        bins[1].eccentricity = outer_eccentricity
        bins[1].argument_of_pericenter = outer_argument_of_pericenter                
        bins[1].longitude_of_ascending_node = outer_longitude_of_ascending_node
        
        bins[1].mass_transfer_rate = 0.0 | units.MSun/units.yr        
        bins[1].accretion_efficiency_mass_transfer = 1.0
        bins[1].accretion_efficiency_wind_child1_to_child2 = 0.0
        bins[1].accretion_efficiency_wind_child2_to_child1 = 0.0

        # binary evolutionary settings
        bins[0].specific_AM_loss_mass_transfer = 2.5 
        bins[1].specific_AM_loss_mass_transfer = 2.5

        return bins
    #-------
            
    #-------
    #setup community codes
    def setup_stellar_code(self, metallicity, stars):
        self.stellar_code = SeBa()
#        self.stellar_code = SeBa(redirection='none')

        #stopping conditions:
#        self.stellar_code.stopping_conditions.supernova_detection.enable()        

        self.stellar_code.parameters.metallicity = metallicity
        self.stellar_code.particles.add_particles(stars)
        self.channel_from_stellar = self.stellar_code.particles.new_channel_to(stars)
        self.channel_to_stellar = stars.new_channel_to(self.stellar_code.particles)
#        self.channel_from_stellar.copy()
        self.channel_from_stellar.copy_attributes(["age", "mass", "core_mass", "radius", "core_radius", "convective_envelope_radius",  "convective_envelope_mass", "stellar_type", "luminosity", "wind_mass_loss_rate",  "temperature"])  #"gyration_radius_sq"

    def setup_secular_code(self, triple_set):
        self.secular_code = SecularTriple()
#        self.secular_code = SecularTriple(redirection='none')
        self.secular_code.triples.add_particles(triple_set)
        self.secular_code.parameters.verbose = False
#        self.secular_code.parameters.verbose = True
        
        self.secular_code.parameters.equations_of_motion_specification = 0
        self.secular_code.parameters.include_quadrupole_terms = True
        self.secular_code.parameters.include_octupole_terms = True        
        self.secular_code.parameters.include_inner_wind_terms = True
        self.secular_code.parameters.include_outer_wind_terms = True
        self.secular_code.parameters.include_inner_RLOF_terms = True
        self.secular_code.parameters.include_outer_RLOF_terms = True
        self.secular_code.parameters.include_magnetic_braking_terms = False # not tested

        self.secular_code.parameters.include_inner_tidal_terms = True
        self.secular_code.parameters.include_outer_tidal_terms = True
        
        self.secular_code.parameters.include_1PN_inner_terms = True
        self.secular_code.parameters.include_1PN_outer_terms = True
        self.secular_code.parameters.include_1PN_inner_outer_terms = False ### warning: probably broken
        self.secular_code.parameters.include_25PN_inner_terms = True
        self.secular_code.parameters.include_25PN_outer_terms = True

        self.secular_code.parameters.check_for_dynamical_stability = True
        self.secular_code.parameters.check_for_dynamical_stability_at_initialisation = True

        self.secular_code.parameters.check_for_semisecular_regime = self.stop_at_semisecular_regime
        self.secular_code.parameters.check_for_semisecular_regime_at_initialisation = self.stop_at_semisecular_regime
        
        self.secular_code.parameters.check_for_inner_collision = True
        self.secular_code.parameters.check_for_outer_collision = True

        self.secular_code.parameters.check_for_inner_RLOF = True 
        self.secular_code.parameters.check_for_outer_RLOF = True 
        
        self.secular_code.parameters.include_spin_radius_mass_coupling_terms_star1 = True
        self.secular_code.parameters.include_spin_radius_mass_coupling_terms_star2 = True
        self.secular_code.parameters.include_spin_radius_mass_coupling_terms_star3 = True
        
         # accuracy of secular code
#        self.secular_code.parameters.input_precision = 1.0e-10#1.0e-5
#        self.secular_code.parameters.relative_tolerance = 1.0e-10
#        self.secular_code.parameters.threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = 1.0e-12
        self.secular_code.parameters.threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero = 1.0e-7|units.Myr**-1


        self.secular_code.parameters.include_linear_mass_change = True #needed for Jspin conservation
        self.secular_code.parameters.include_linear_radius_change = True #needed for Jspin conservation

        self.channel_from_secular = self.secular_code.triples.new_channel_to(triple_set)
        self.channel_to_secular = triple_set.new_channel_to(self.secular_code.triples)
    #-------

    #-------
    def test_initial_parameters(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node):
            

#         if max(inner_primary_mass, outer_mass) > max_mass:  
#             print(inner_primary_mass, inner_secondary_mass, outer_mass)
#             print('should be within:', min_mass, '-', max_mass)
#             print('error: masses not in allowed range')
#             return False, 0,0
                        
#         if min(inner_secondary_mass, outer_mass) <= absolute_min_mass:  
#             print(inner_primary_mass, inner_secondary_mass, outer_mass)
#             print('should be at least above:', absolute_min_mass)
#             print('error: masses not in allowed range')
#             return False, 0,0

        if inner_semimajor_axis >= outer_semimajor_axis:
            print('error input parameters, should be:')
            print('inner_semimajor_axis < outer_semimajor_axis' )
            return False, 0,0
        if (inner_semimajor_axis < 0.|units.RSun):
            print('error: inner separation not in allowed range')
            return False, 0,0
        if (outer_semimajor_axis < 0.|units.RSun):
            print('error: outer separation not in allowed range')
            return False, 0,0
    
        if (inner_eccentricity < 0.) or (inner_eccentricity > 1.):
            print('error: inner eccentricity not in allowed range')
            return False, 0,0
        if (outer_eccentricity < 0.) or (outer_eccentricity > 1.):
            print('error: outer eccentricity not in allowed range')
            return False, 0,0
        if (inner_eccentricity < minimum_eccentricity):
            inner_eccentricity = minimum_eccentricity
        if (outer_eccentricity < minimum_eccentricity):
            outer_eccentricity = minimum_eccentricity
    
        if (relative_inclination < 0.) or (relative_inclination > np.pi):
            print('error: relative inclination not in allowed range')
            return False, 0,0
    
        if (inner_argument_of_pericenter < -1.*np.pi) or (inner_argument_of_pericenter > np.pi):
            print('error: inner argument of pericenter not in allowed range')
            return False, 0,0
        if (outer_argument_of_pericenter < -1.*np.pi) or (outer_argument_of_pericenter > np.pi):
            print('error: outer argument of pericenter not in allowed range')
            return False, 0,0
    
        if (inner_longitude_of_ascending_node < -1.*np.pi) or (inner_longitude_of_ascending_node > np.pi):
            print('error: inner longitude of ascending node not in allowed range')
            return False, 0,0
            
        return True, inner_eccentricity, outer_eccentricity 
        
        
    def initial_angular_frequency(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        self.previous_time = self.triple.time
        if stellar_system.is_star:
            stellar_system.spin_angular_frequency = lang_spin_angular_frequency(stellar_system)
#            stellar_system.spin_angular_frequency = corotating_spin_angular_frequency_binary(stellar_system.parent.semimajor_axis, self.get_mass(stellar_system.parent.child1), self.get_mass(stellar_system.parent.child2))
        else:
            self.initial_angular_frequency(stellar_system.child1)        
            self.initial_angular_frequency(stellar_system.child2)
                               
    #-------

    #-------
    def refresh_memory(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        self.triple.time = self.previous_time
        if stellar_system.is_star:
            stellar_system.mass = stellar_system.previous_mass 
            stellar_system.radius = stellar_system.previous_radius           
            stellar_system.stellar_type = stellar_system.previous_stellar_type
            stellar_system.moment_of_inertia_of_star = stellar_system.previous_moment_of_inertia_of_star
            stellar_system.time_derivative_of_radius = stellar_system.previous_time_derivative_of_radius
            stellar_system.spin_angular_frequency = stellar_system.previous_spin_angular_frequency
            
        else:
            self.refresh_memory(stellar_system.child1)        
            self.refresh_memory(stellar_system.child2)
            stellar_system.mass = stellar_system.previous_mass
            stellar_system.bin_type = stellar_system.previous_bin_type
            
            if self.is_triple(stellar_system):
                stellar_system.kozai_type = stellar_system.previous_kozai_type

    def update_previous_stellar_parameters(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        self.previous_time = self.triple.time
        if stellar_system.is_star:
            stellar_system.previous_mass = self.get_mass(stellar_system)      
            stellar_system.previous_radius = stellar_system.radius
            stellar_system.previous_stellar_type = stellar_system.stellar_type
            stellar_system.previous_moment_of_inertia_of_star = stellar_system.moment_of_inertia_of_star
            stellar_system.previous_spin_angular_frequency = stellar_system.spin_angular_frequency 
            if self.triple.time == quantities.zero: #initialization
               stellar_system.previous_time_derivative_of_radius = 0.0 | units.RSun/units.yr
            else:
               stellar_system.previous_time_derivative_of_radius = stellar_system.time_derivative_of_radius
        else:
            self.update_previous_stellar_parameters(stellar_system.child1)        
            self.update_previous_stellar_parameters(stellar_system.child2)
            stellar_system.previous_mass = self.get_mass(stellar_system) 
            stellar_system.previous_bin_type = stellar_system.bin_type
            
            if self.is_triple(stellar_system):
                stellar_system.previous_kozai_type = stellar_system.kozai_type
                
    #-------

    #-------
    def update_time_derivative_of_radius(self, stellar_system = None):
        #update time_derivative_of_radius for effect of wind on spin
        #radius change due to stellar evolution, not mass transfer
        if stellar_system == None:
            stellar_system = self.triple
                
        time_step = self.triple.time - self.previous_time

        if self.triple.time == quantities.zero:
            #initialization
            self.triple.child2.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.triple.child2.child2.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.triple.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
        else:     
            if stellar_system.is_star:
                stellar_system.time_derivative_of_radius = (stellar_system.radius - stellar_system.previous_radius)/time_step
            else:
                self.update_time_derivative_of_radius(stellar_system.child1)        
                self.update_time_derivative_of_radius(stellar_system.child2)
    #-------
    
    #-------
    def update_time_derivative_of_moment_of_inertia(self, stellar_system = None):
        #update time_derivative_of_radius for effect of changing Jspin
        if stellar_system == None:
            stellar_system = self.triple
                
        time_step = self.triple.time - self.previous_time

        if self.triple.time == quantities.zero:
            #initialization
            self.triple.child2.child1.time_derivative_of_moment_of_inertia = 0.0 | units.MSun*units.RSun**2/units.yr
            self.triple.child2.child2.time_derivative_of_moment_of_inertia = 0.0 | units.MSun*units.RSun**2/units.yr
            self.triple.child1.time_derivative_of_moment_of_inertia = 0.0 | units.RSun/units.yr
        else:     
            if stellar_system.is_star:
                stellar_system.time_derivative_of_moment_of_inertia = (stellar_system.moment_of_inertia_of_star - stellar_system.previous_moment_of_inertia_of_star)/time_step
            else:
                self.update_time_derivative_of_moment_of_inertia(stellar_system.child1)        
                self.update_time_derivative_of_moment_of_inertia(stellar_system.child2)
    #-------

    #-------
    def update_stellar_parameters(self, stellar_system = None):
        # for the convective envelope mass:
        # the prescription of Hurley, Pols & Tout 2000 is implemented in SeBa, however note that the prescription in BSE is different
        # for the convective envelope radius:
        # the prescription of Hurley, Tout & Pols 2002 is implemented in SeBa, however note that the prescription in BSE is different
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            
            stellar_system.gyration_radius = 0.
            stellar_system.apsidal_motion_constant = self.apsidal_motion_constant(stellar_system) 
            if stellar_system.core_radius > stellar_system.radius:
                #can happen very late on the agb before WD formation
                stellar_system.core_radius = stellar_system.radius                
            stellar_system.moment_of_inertia_of_star = self.moment_of_inertia(stellar_system)

            if stellar_system.convective_envelope_radius < 0|units.RSun:
                print('convective_envelope_radius < 0')
                exit(1)
            if stellar_system.convective_envelope_radius == 0|units.RSun:
                stellar_system.convective_envelope_mass = 1.e-10 |units.MSun    
                stellar_system.convective_envelope_radius = 1.e-10 |units.RSun   
                 
        else:
            self.update_stellar_parameters(stellar_system.child1)        
            self.update_stellar_parameters(stellar_system.child2)
            if self.is_triple(stellar_system):
                stellar_system.kozai_type = self.get_kozai_type()           
    #-------

    #-------
    # useful functions general
    
    #whether or not a stellar system consists of just two stars
    def is_binary(self, stellar_system=None):
        if stellar_system == None:
            stellar_system = self.triple    
    
        if not stellar_system.is_star and stellar_system.child1.is_star and stellar_system.child2.is_star:
            return True
        else:
            return False

    def is_triple(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if not stellar_system.is_star:
            if stellar_system.child1.is_star and self.is_binary(stellar_system.child2):
                return True
            elif stellar_system.child2.is_star and self.is_binary(stellar_system.child1):
                return True

        return False

    def has_donor(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            if stellar_system.is_donor:
                return True
        else:
            if self.has_donor(stellar_system.child1) or self.has_donor(stellar_system.child2):
                return True                        
            
        return False            


#    def has_contact_system(self, stellar_system = None):
#        if stellar_system == None:
#            stellar_system = self.triple
#            
#        if stellar_system.is_star:
#            return False
#        elif self.is_binary(stellar_system):
#            if stellar_system.child1.is_donor and stellar_system.child2.is_donor:
#                return True
#        else:
#            if self.has_contact_system(stellar_system.child1):
#                return True
#            if self.has_contact_system(stellar_system.child2):
#                return True
#            
#        return False            

# if a merger is currently taking place, not if a merger has happened in the past
    def has_merger(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        else:
            if self.has_merger(stellar_system.child1):
                return True
            if self.has_merger(stellar_system.child2):
                return True
            if stellar_system.bin_type == bin_type['merger']:  
                return True    

        return False            

# if a disruption is currently taking place, not if a disruption has happened in the past
    def has_disintegrated(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        else:
            if self.has_disintegrated(stellar_system.child1):
                return True
            if self.has_disintegrated(stellar_system.child2):
                return True
            if stellar_system.bin_type == bin_type['disintegrated']:  
                return True    
            
        return False            

#doesn't work well, as it uses bin_types that are set later -> use has_tertiary_donor
# if a mass transfer in the outer binary of the triple is currently taking place, not if a mass transfer has happened in the past
#    def has_outer_mass_transfer(self, stellar_system = None): 
#        if stellar_system == None:
#            stellar_system = self.triple
#            
#        if stellar_system.is_star:
#            return False
#        elif self.is_binary(stellar_system):
#            return False
#        else:
#            if self.has_outer_mass_transfer(stellar_system.child1):
#                return True
#            if self.has_outer_mass_transfer(stellar_system.child2):
#                return True
#            if stellar_system.bin_type != bin_type['unknown'] and stellar_system.bin_type != bin_type['detached']:  
#                return True    
#            
#        return False            


# if a mass transfer in the outer binary of the triple is currently taking place, not if a mass transfer has happened in the past
    def has_tertiary_donor(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        elif self.is_binary(stellar_system):
            return False
        else:
            if self.has_tertiary_donor(stellar_system.child1): 
                return True
            if self.has_tertiary_donor(stellar_system.child2): 
                return True
            if stellar_system.child1.is_star and stellar_system.child1.is_donor and not stellar_system.child2.is_star:
                return True 
            if stellar_system.child2.is_star and stellar_system.child2.is_donor and not stellar_system.child1.is_star:
                return True             
        return False            
       
    def contains_SN_remnant(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type in stellar_types_SN_remnants:
                return True
        else:
            if self.contains_SN_remnant(stellar_system.child1) or self.contains_SN_remnant(stellar_system.child2):
                return True                        
            
        return False            
    

    def has_stellar_type_changed_into_SN_remnant(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_SN_remnants:
                return True
        else:
            if self.has_stellar_type_changed_into_SN_remnant(stellar_system.child1) or self.has_stellar_type_changed_into_SN_remnant(stellar_system.child2):
                return True                        
            
        return False            


    def has_stellar_type_changed(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type:
                return True
        else:
            if self.has_stellar_type_changed(stellar_system.child1) or self.has_stellar_type_changed(stellar_system.child2):
                return True                        
            
        return False            
    
    

    def has_kozai_type_changed(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if self.is_triple(stellar_system):
            if stellar_system.kozai_type != stellar_system.previous_kozai_type:
                return True
        else: #binaries and single stars do not have a kozai timescale
            return False            

#obsolete?
#    def is_system_stable(self, stellar_system = None):
#        if stellar_system == None:
#            stellar_system = self.triple
#            
#        if stellar_system.is_star:
#            return True
#        elif self.is_binary(stellar_system):
#            return stellar_system.is_mt_stable
#        else:
#            if stellar_system.is_mt_stable and self.is_system_stable(stellar_system.child1) and self.is_system_stable(stellar_system.child2):
#                return True                        
#            
#        return False                    

    def get_mass(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            return stellar_system.mass
        else:
            M1 = self.get_mass(stellar_system.child1)        
            M2 = self.get_mass(stellar_system.child2)
            return M1 + M2 
    #-------

    #-------
    # useful functions general
            
    def orbital_period(self, bs):
        if not bs.is_star:
            Porb = 2*np.pi * np.sqrt(bs.semimajor_axis**3/constants.G / self.get_mass(bs))
            return Porb
        else:
            print('orbital_period: single star does not have a period')
            exit(-1)

    def orbital_angular_momentum(self, bs):
        if not bs.is_star:
            M = self.get_mass(bs.child1)
            m = self.get_mass(bs.child2)
            a = bs.semimajor_axis
            e = bs.eccentricity
            J = M*m * np.sqrt(constants.G*a*(1-e**2)/(M+m))
        
            if REPORT_BINARY_EVOLUTION:
                print('Jorb:', M, m, a, e, J)
                
            return J
        else:
            print('orbital_angular_momentum: single star does not have an orbit')
            exit(-1)
    
    def spin_angular_momentum(self, ss):
        if ss.is_star:
            return ss.moment_of_inertia_of_star * ss.spin_angular_frequency
        else:
            print('spin_angular_momentum: structure stellar system unknown')        
            exit(2)
            
    def apsidal_motion_constant(self, star):
    
        if star.stellar_type in [13]|units.stellar_type: #ns
            #based on Brooke & Olle 1955, for n=1 polytrope
            return 0.260
    
        elif star.stellar_type in [13, 14]|units.stellar_type: #bh
            # Hamers et al. 2013
            return 0.
        elif star.stellar_type in [1,7,10,11,12]|units.stellar_type:#ms, he-ms, wd
            #based on Brooke & Olle 1955, for n=3 polytrope
            return 0.0144            

        elif star.stellar_type in [0,2,3,4,5,6,8,9,17]|units.stellar_type:#low-mass ms, hg, gb, cheb, agb, he-g, pre-ms
            #based on Brooke & Olle 1955, for n=3 polytrope
#            return 0.143 
            #based on Claret & Gimenez 1992, 96, 225 the value should be smaller, try:
            return 0.05
        elif star.stellar_type in [18]|units.stellar_type:#planet
            #based on Brooke & Olle 1955, for n=3 polytrope
            return 0.0144            
        elif star.stellar_type in [19]|units.stellar_type:#bd
            #based on Brooke & Olle 1955, for n=3 polytrope
            return 0.0144
        else:
            print('apsidal motion constant: stellar_type unknown')
            print(star.stellar_type)
            exit(2)
            
    #Hurley, Pols, Tout 2000
    def moment_of_inertia(self, star):
        if star.is_star:
            k2 = 0.1
            k3 = 0.21
        
            if star.stellar_type in stellar_types_remnants:
                I = k3*(star.mass)*star.radius**2 
            else:            
                I = k2*(star.mass - star.core_mass)*star.radius**2 + k3*star.core_mass*star.core_radius**2
                
            return I                   
        else:
            print('moment_of_inertia: structure stellar system unknown')        
            exit(2)


    def octupole_parameter(self):
        if self.is_triple():
           if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
           else: 
                star = self.triple.child2
                bin = self.triple.child1

           return (self.get_mass(bin.child1)-self.get_mass(bin.child2))/self.get_mass(bin) * bin.semimajor_axis/self.triple.semimajor_axis * self.triple.eccentricity/(1-self.triple.eccentricity**2)           

        else:
            print('Octupole parameter needs triple system')
            return np.nan   




    def kozai_timescale(self):
        if self.is_triple():
           alpha_kozai = 1.
           if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
           else: 
                star = self.triple.child2
                bin = self.triple.child1

            
           P_in = self.orbital_period(bin) #period inner binary 
           P_out = self.orbital_period(self.triple)#period outer binary 
           return alpha_kozai * P_out**2 / P_in * (self.get_mass(self.triple) / self.get_mass(star)) * (1-self.triple.eccentricity**2)**1.5       

        else:
            print('Kozai timescale needs triple system')
            return np.nan   
            
    def get_kozai_type(self):
        if self.is_triple():
            if self.secular_code.parameters.ignore_tertiary == True:
                return False

            t_kozai = self.kozai_timescale()
            if t_kozai >kozai_type_factor * self.tend:
                return False
            
            t_ev = self.get_min_stellar_evolution_timescale_of_system()
            if t_kozai < kozai_type_factor * t_ev:
                return True
            else:
                return False
        else:
           print('Kozai type needs triple system')
           exit(1)   

    def get_min_stellar_evolution_timescale_of_system(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            return stellar_evolution_timescale(stellar_system)
        else:
            t1 = self.get_min_stellar_evolution_timescale_of_system(stellar_system.child1)        
            t2 = self.get_min_stellar_evolution_timescale_of_system(stellar_system.child2)
            return max(t1, t2)


    def check_RLOF(self):
        if self.triple.is_star:
            return
        elif self.is_binary():
            Rl1 = roche_radius(self, self.child1)
            Rl2 = roche_radius(self, self.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print('Roche lobe radii:', Rl1, Rl2)
                print('Stellar radii:', self.triple.child1.radius, self.triple.child2.radius)
            
            self.triple.child1.is_donor = False
            self.triple.child2.is_donor = False
            
            if self.triple.child1.radius >= Rl1:
                self.triple.child1.is_donor = True               
            if self.triple.child2.radius >= Rl2:
                self.triple.child2.is_donor = True                             
 
        elif self.is_triple() and self.secular_code.parameters.ignore_tertiary == True:
            # for distrupted binary
            if self.triple.child1.is_star:
                bin = self.triple.child2
            else:
                bin = self.triple.child1
        
            Rl1 = roche_radius(bin, bin.child1, self)
            Rl2 = roche_radius(bin, bin.child2, self)
            if REPORT_TRIPLE_EVOLUTION:
                print('Roche lobe radii:', Rl1, Rl2)
                print('Stellar radii:', bin.child1.radius, bin.child2.radius)
            
            bin.child1.is_donor = False
            bin.child2.is_donor = False
            
            if bin.child1.radius >= Rl1:
                bin.child1.is_donor = True               
            if bin.child2.radius >= Rl2:
                bin.child2.is_donor = True                             
 
 
        elif self.is_triple():
            if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
            else:
                star = self.triple.child2
                bin = self.triple.child1

            #assumping secular code always returns inner binary first
            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
    
            if REPORT_TRIPLE_EVOLUTION:
                print('Roche lobe radii:', Rl1, Rl2, Rl3)
                print('Stellar radii:', bin.child1.radius, bin.child2.radius, star.radius)
                print('binary Roche lobe radii:', roche_radius(bin, bin.child1, self), roche_radius(bin, bin.child2, self), roche_radius(self.triple, star, self))
                print('eccentric binary Roche lobe radii:', roche_radius(bin, bin.child1, self)* (1-bin.eccentricity), roche_radius(bin, bin.child2, self)* (1-bin.eccentricity), roche_radius(self.triple, star, self)*(1-self.triple.eccentricity))
                print('Masses:', bin.child1.mass, bin.child2.mass, star.mass)
                print('Semi:', bin.semimajor_axis, self.triple.semimajor_axis)
                print('Ecc:', bin.eccentricity, self.triple.eccentricity)
                print('Stellar type:', bin.child1.stellar_type, bin.child2.stellar_type, star.stellar_type)
                print('Spin:', bin.child1.spin_angular_frequency, bin.child2.spin_angular_frequency, star.spin_angular_frequency)
                

            bin.child1.is_donor = False
            bin.child2.is_donor = False
            star.is_donor = False
                
            if bin.child1.radius >= Rl1:
                bin.child1.is_donor = True
            if bin.child2.radius >= Rl2:
                bin.child2.is_donor = True
            if star.radius >= Rl3:
                star.is_donor = True

            if star.is_donor and (bin.child1.is_donor or bin.child2.is_donor):
                print('RLOF in inner and outer binary')
                print(Rl1, bin.child1.radius, Rl2, bin.child2.radius)
                print(RL3, star.radius)
                exit(1)                   
                
        else:
            print('check_RLOF: structure stellar system unknown')
            exit(2)    
                     
#            
#    def determine_partial_timestep_stable_mass_transfer(self, stellar_system = None):
#        if stellar_system == None:
#            stellar_system = self.triple
#
#        if stellar_system.is_star:
#            return np.inf |units.Myr 
#        else:
#            dt1 = self.determine_partial_timestep_stable_mass_transfer(stellar_system.child1)        
#            dt2 = self.determine_partial_timestep_stable_mass_transfer(stellar_system.child2)
#            dt =  stellar_system.part_dt_mt
#            return min(dt, min(dt1, dt2))
                   
    #-------

    #-------
    # useful functions for printing
    def print_star(self, star):
        if star.is_star:
            print('star:')
            print(star.age, )
            print(star.stellar_type, )
            print(star.mass, )
            print(star.radius, )
            print(star.core_mass, )
            print(star.core_radius,)
            print(star.convective_envelope_mass,)
            print(star.convective_envelope_radius,)
            print(star.luminosity,)
            print(star.temperature,)
            print(star.wind_mass_loss_rate,)
            print(star.spin_angular_frequency,)
            print(star.is_donor)
            print('\t')
        else:
            print('print_star needs a star')
            exit(2)
    
    
    def print_binary(self, binary):
        if not binary.is_star:
            print(self.get_mass(binary), )
            print(binary.semimajor_axis,) 
            print(binary.eccentricity, )
            print(binary.argument_of_pericenter, )
            print(binary.longitude_of_ascending_node,)
            print(binary.mass_transfer_rate,)
#            print(binary.mass_transfer_timescale,)
            print(binary.accretion_efficiency_mass_transfer,)
            print(binary.accretion_efficiency_wind_child1_to_child2,)
            print(binary.accretion_efficiency_wind_child2_to_child1,)
            print(binary.specific_AM_loss_mass_transfer,)
            print(binary.is_mt_stable)
            print('\t')
        else:
            print('print_binary needs a binary')        
            exit(2)
    
    
    def print_stellar_system(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            print(stellar_system.number,)
            print(stellar_system.relative_inclination)

        if stellar_system.is_star:
            self.print_star(stellar_system)
        else:
            print('binary star: ')
            self.print_binary(stellar_system)
            self.print_stellar_system(stellar_system.child1)
            self.print_stellar_system(stellar_system.child2)
        print('\t')
    #-------
    #-------
    #don't change this unless you know what you're doing
    def remove_parents(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        parents = []
        if stellar_system.is_star:
            try:
                p = stellar_system.parent
                stellar_system.parent = 0
                return p
            except AttributeError: #when there is no parent
                return 0
        else:
            parents.append(self.remove_parents(stellar_system.child1))
            parents.append(self.remove_parents(stellar_system.child2))

            p = stellar_system.parent
#            except AttributeError: #when there is no parent=
            if p != None:
                stellar_system.parent = 0
                parents.append(p)                                    

            return parents
            
    #don't change this unless you know what you're doing
    def set_parents(self, parents, stellar_system=None):            
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
                stellar_system.parent = parents
        else:
            self.set_parents(parents[0], stellar_system.child1)
            self.set_parents(parents[1], stellar_system.child2)
            if len(parents) == 3:
                stellar_system.parent = parents[2]
            elif len(parents) != 2:
                print('set_parents: structure stellar system unknown') 
                exit(2)
 
    def save_snapshot(self):
        file_name = self.file_name

        if self.file_type == 'txt':
            print(self.file_name,self.file_type)
            parents = self.remove_parents()
            write_set_to_file(self.triple.as_set(), self.file_name, self.file_type) 
            self.set_parents(parents)

        else:
            write_set_to_file(self.triple.as_set(), self.file_name, self.file_type, version='2.0', append_to_file=True) 

    #some minor parameters are missing:
#        self.instantaneous_evolution = False 
#        self.tend = tend #...
#        self.triple.time = 0.0|units.yr
#        self.previous_time = 0.0|units.yr
#        self.file_name
#        self.file_type

    #-------
        
    #-------
    # determining time steps
    def determine_time_step_wind(self, stellar_system = None):
    #note: returned value can be inf when the wind_mass_loss_rate = 0
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr

            if stellar_system.wind_mass_loss_rate * -1. > quantities.zero:
                dt = maximum_wind_mass_loss_factor * (stellar_system.mass - stellar_system.core_mass)/stellar_system.wind_mass_loss_rate*-1.  
            if REPORT_DT:
                print("Dt_wind_star = ", dt)
            return dt 
        else:
            dt1 = self.determine_time_step_wind(stellar_system.child1)        
            dt2 = self.determine_time_step_wind(stellar_system.child2)
            if REPORT_DT:
                print("Dt_wind_binary = ", dt1, dt2)
            return min(dt1, dt2) 
    


    def determine_time_step_radius_change(self, stellar_system = None):
        #note: returned value can be inf when the change in radius <= 0
        #radius is only necessary for tides

        if not self.secular_code.parameters.include_inner_tidal_terms and not self.secular_code.parameters.include_outer_tidal_terms:
            return np.inf |units.Myr

        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            if stellar_system.time_derivative_of_radius == quantities.zero:
                dt = np.inf |units.Myr
            else:     
                if stellar_system.previous_time_derivative_of_radius == quantities.zero:
                    growth_factor = 0.1
                elif stellar_system.previous_time_derivative_of_radius * stellar_system.time_derivative_of_radius < quantities.zero:
                    growth_factor = 0.01
                else:
                    growth_factor = stellar_system.previous_time_derivative_of_radius/stellar_system.time_derivative_of_radius 
                    if growth_factor > 1:
                        growth_factor = 1./growth_factor
                
                if stellar_system.stellar_type > 1|units.stellar_type and stellar_system.time_derivative_of_radius < stellar_system.previous_time_derivative_of_radius: #not a MS star and Rdot < Rdot_prev
                    growth_factor = 1. 

                dt = abs(growth_factor * self.maximum_radius_change_factor*stellar_system.radius / stellar_system.time_derivative_of_radius)


            if REPORT_DT:
                print("Dt_radius_change_star = ", dt)
            return dt 
        else:
            dt1 = self.determine_time_step_radius_change(stellar_system.child1)        
            dt2 = self.determine_time_step_radius_change(stellar_system.child2)
            if REPORT_DT:
                print("Dt_radius_change_binary = ", dt1, dt2)
            return min(dt1, dt2) 

    def determine_time_step_kozai(self):
    #note: returned value can be inf when the system is a binary or single star            
        if self.is_triple():
            dt = self.kozai_timescale()*time_step_factor_kozai
            if REPORT_DT:
                print("Dt_kozai = ", dt)
        else:
            dt = np.inf |units.Myr
        
        P_out = self.orbital_period(self.triple)
        P_in = self.orbital_period(self.triple.child2)        
        return dt

 
    def determine_time_step_stable_mt(self, stellar_system = None):
    #note: returned value can be inf when the mass_transfer_rate = 0 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.is_donor:
                dt = abs(time_step_factor_stable_mt*stellar_system.mass/stellar_system.parent.mass_transfer_rate)
            if REPORT_DT:
                print("Dt_mt_star = ", dt, time_step_factor_stable_mt, stellar_system.mass, stellar_system.parent.mass_transfer_rate )
            return dt
        else:
            dt1 = self.determine_time_step_stable_mt(stellar_system.child1)        
            dt2 = self.determine_time_step_stable_mt(stellar_system.child2)
            if REPORT_DT:
                print("Dt_mt_binary = ", dt1, dt2)
            return min(dt1, dt2) 
         
         
    #based on secular code      
    def e_dot_tides(self, star,m_comp, semi, eccentricity):
        e_p2 = eccentricity**2
        l = np.sqrt(1.0 - e_p2)        
        f_tides3 = 1.0 + e_p2*(15./4. + e_p2*(15./8. + e_p2*5./64.))
        f_tides4 = 1.0 + e_p2*(3./2 + e_p2/8.)
    
        R_div_a = star.radius/semi        
        spin = star.spin_angular_frequency
        k_div_T_tides = tidal_friction_constant(star.stellar_type, star.mass, m_comp, semi, star.radius, star.convective_envelope_mass, star.convective_envelope_radius, star.luminosity, spin)
        n = corotating_spin_angular_frequency_binary(semi, star.mass, m_comp) # mean orbital angular speed

        e_dot = -27.0*(1.0+m_comp/star.mass)*(m_comp/star.mass) * R_div_a**6 * k_div_T_tides * R_div_a**2 *eccentricity* (l**-13.0) * (f_tides3 - 11.0/18.0*(l**3)*f_tides4*(spin/n))
#        print(star.mass, m_comp, semi, eccentricity, k_div_T_tides, e_dot)
        return e_dot
         
         
         
    def determine_time_step_tides(self, stellar_system = None):
#        print("determine_time_step_tides")
        if stellar_system == None:
            stellar_system = self.triple


        if self.triple.is_star:
            return False, 0        
        elif self.is_binary():
            Rl1 = roche_radius(self, self.child1)
            Rl2 = roche_radius(self, self.child2)
            if self.triple.child1.radius >= Rl1 or self.triple.child2.radius >= Rl2:
                return abs(time_step_factor_stable_mt*min(self.triple.child1.mass, self.triple.child2.mass)/self.triple.mass_transfer_rate)
            else: 
                de_dt_child1 = abs(self.e_dot_tides(self.triple.child1, self.triple.child2.mass, self.triple.semimajor_axis, self.triple.eccentricity))
                de_dt_child2 = abs(self.e_dot_tides(self.triple.child2, self.triple.child1.mass, self.triple.semimajor_axis, self.triple.eccentricity))
                return time_step_factor_ecc * self.triple.eccentricity / max(de_dt_child1, de_dt_child2)
        elif self.is_triple():
        
            if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
            else:
                star = self.triple.child2
                bin = self.triple.child1

            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
            if star.radius >= Rl3:
                exit(2)  
            elif bin.child1.radius >= Rl1 or bin.child2.radius >= Rl2:
                return abs(time_step_factor_stable_mt*min(bin.child1.mass, bin.child2.mass)/self.triple.mass_transfer_rate)
            else:
                de_dt_star = abs(self.e_dot_tides(star, self.get_mass(bin), self.triple.semimajor_axis, self.triple.eccentricity))
                de_dt_bin_child1 = abs(self.e_dot_tides(bin.child1, bin.child2.mass, bin.semimajor_axis, bin.eccentricity))
                de_dt_bin_child2 = abs(self.e_dot_tides(bin.child2, bin.child1.mass, bin.semimajor_axis, bin.eccentricity))
                dt_star = time_step_factor_ecc * self.triple.eccentricity / de_dt_star
                dt_bin = time_step_factor_ecc * bin.eccentricity / max(de_dt_bin_child1, de_dt_bin_child2)
                return min(dt_star, dt_bin)

        else: 
            print('determine_time_step_tides: structure stellar system unknown')        
            exit(2)    

    
         
         
#    def close_to_RLOF(self):
#
#        if self.triple.is_star:
#            return False, 0        
#        elif self.is_binary():
#            Rl1 = roche_radius(self, self.child1)
#            Rl2 = roche_radius(self, self.child2)
#            if self.triple.child1.radius >= Rl1 or self.triple.child2.radius >= Rl2:
#                print('close to rlof already rlof')
#                exit(0)
#            elif self.triple.child1.radius >= Rl_fraction*Rl1 and self.triple.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof = maximum(self.triple.child1.radius/Rl1, self.triple.child2.radius/Rl2)
#                print('dt_close_to_mt', time_step_factor_stable_mt, self.triple.child1.mass, self.triple.child2.mass, self.triple.mass_transfer_rate)
#                dt = abs(time_step_factor_stable_mt*min(self.triple.child1.mass, self.triple.child2.mass)/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt  
#            elif self.triple.child1.radius >= Rl_fraction*Rl1:
#                ratio_rad_rlof = self.triple.child1.radius/Rl1
#                print('dt_close_to_mt', time_step_factor_stable_mt, self.triple.child1.mass, self.triple.mass_transfer_rate)
#                dt = abs(time_step_factor_stable_mt*self.triple.child1.mass/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt  
#            elif self.triple.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof = self.triple.child2.radius/Rl2
#                print('dt_close_to_mt', time_step_factor_stable_mt, self.triple.child2.mass, self.triple.mass_transfer_rate)
#                dt = abs(time_step_factor_stable_mt*self.triple.child2.mass/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt              
#            else:
#                return False, 0, 0    
#        elif self.is_triple():
#            if self.triple.child1.is_star:
#                star = self.triple.child1
#                bin = self.triple.child2
#            else:
#                star = self.triple.child2
#                bin = self.triple.child1
#
#            #assumping secular code always returns inner binary first
#            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
#            if bin.child1.radius >= Rl1 or bin.child2.radius >= Rl2 or star.radius >= Rl3:
#                print('close to rlof already rlof')
#                exit(0)
#
#            dt=[]
#            ratio_rad_rlof=[]
#            if star.radius >= Rl_fraction*Rl3:
#                ratio_rad_rlof.append(star.radius / Rl3)
#                print('dt_close_to_mt', time_step_factor_stable_mt, star.mass, star.parent.mass_transfer_rate)
#                dt.append(abs(time_step_factor_stable_mt*star.mass/star.parent.mass_transfer_rate))
#            if bin.child1.radius >= Rl_fraction*Rl1:
#                ratio_rad_rlof.append(bin.child1.radius / Rl1)
#                print('dt_close_to_mt', time_step_factor_stable_mt, bin.child1.mass, bin.mass_transfer_rate)
#                dt.append(abs(time_step_factor_stable_mt*bin.child1.mass/bin.mass_transfer_rate))
#            if bin.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof.append(bin.child2.radius / Rl2)
#                print('dt_close_to_mt', time_step_factor_stable_mt, bin.child2.mass, bin.mass_transfer_rate)
#                dt.append(abs(time_step_factor_stable_mt*bin.child2.mass/bin.mass_transfer_rate))
#                    
#            if len(dt) > 0:
#                print('close to rlof', ratio_rad_rlof, dt)
#                return True, max(ratio_rad_rlof), min(dt)
#            else:
#                return False, 0, 0
#        else:
#            print('close_to_RLOF: structure stellar system unknown')        
#            exit(2)    
    
        
    def determine_time_step(self):         
        if REPORT_DT:
            print("Dt = ", self.stellar_code.particles.time_step, self.tend, self.previous_time)

        #maximum time_step            
        time_step_max = self.tend - self.triple.time        
                               
        # time_step of stellar evolution
        time_step_stellar_code = self.stellar_code.particles.time_step

        # small time_step after type change
        # done automatically by SeBa when due to wind mass losses, but not during RLOF
        if self.has_stellar_type_changed():
            time_step_stellar_code.append(minimum_time_step)                    
                    
        # small time_step during heavy wind mass losses
        time_step_wind =  self.determine_time_step_wind()

        # small time_step during phases of fast growth (note: not yet during fast shrinkage)
        time_step_radius_change = self.determine_time_step_radius_change()             
        
#        #during stable mass transfer   
#        if self.has_donor():
##            print(time_step, self.determine_time_step_stable_mt())
#            time_step = min(time_step, self.determine_time_step_stable_mt())
#            if REPORT_DT or REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
#                print('donor time:',  self.determine_time_step_stable_mt())

        #tides - ADD ALSO DURING MT
        time_step_tides = np.inf |units.Myr 
	#interesting alternative, slows down code 
#         if self.secular_code.parameters.include_inner_tidal_terms or self.secular_code.parameters.include_outer_tidal_terms:    
#             time_step_tides = self.determine_time_step_tides()  	
                
        if REPORT_DT or REPORT_DEBUG:
            print('time:', time_step_max, time_step_stellar_code, time_step_wind, time_step_radius_change, time_step_tides)



        if not self.has_donor():
            if REPORT_DT:
                print('no RLOF' )
            time_step = min(time_step_radius_change, min(time_step_wind, min( min(time_step_stellar_code), time_step_max)))
            if time_step_tides < time_step:
                Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
                if REPORT_DT:
                    print(time_step, time_step_tides)
                    print(Rl1,Rl2,Rl3)
                    print(self.triple.child2.child1.radius, self.triple.child2.child2.radius, self.triple.child1.radius )
                time_step = time_step_tides

            if REPORT_DT or REPORT_DEBUG:
                print('min increase', time_step, maximum_time_step_factor*self.previous_dt)   
                
                                             
            time_step = min(time_step, maximum_time_step_factor*self.previous_dt)  


        elif self.has_donor() and self.triple.bin_type == 'detached' and self.triple.child2.bin_type == 'detached':
            #find beginning of RLOF carefully
            time_step = time_step_factor_find_RLOF*self.previous_dt
            if REPORT_DT:
                print('find rlof', time_step_factor_find_RLOF,self.previous_dt)

            #resetting is_donor
            self.check_RLOF()

        else:
            #during stable mass transfer   
            time_step = min(time_step_wind, min( min(time_step_stellar_code), time_step_max))
            time_step = min(time_step, self.determine_time_step_stable_mt())
            if REPORT_DT or REPORT_DEBUG:
                print('donor time:',  self.determine_time_step_stable_mt())

            if REPORT_DT or REPORT_DEBUG:
                print('min increase', time_step, maximum_time_step_factor*self.previous_dt)                                
            time_step = min(time_step, maximum_time_step_factor*self.previous_dt)  


        if self.triple.child2.bin_type == 'detached' and self.triple.child2.previous_bin_type == 'stable_mass_transfer':
            corotation_spin = corotating_spin_angular_frequency_binary(self.triple.child2.semimajor_axis, self.triple.child2.child1.mass, self.triple.child2.child2.mass)
            self.triple.child2.child1.spin_angular_frequency = corotation_spin
            self.triple.child2.child2.spin_angular_frequency = corotation_spin
            
            if REPORT_DT or REPORT_DEBUG:
                print('prev timestep', time_step, previous_time_step)
            previous_time_step = self.triple.time - self.previous_time
            time_step = min(time_step, maximum_time_step_factor_after_stable_mt*previous_time_step)  



        if self.triple.time == quantities.zero:
            #initialization (e.g. time_derivative_of_radius)
            P_out = self.orbital_period(self.triple) #period outer binary 
            # do not take 0.1*P_in -> resonance -> large error
            time_step = min(min(P_out, time_step), 1.|units.yr)


        time_step = max(time_step, minimum_time_step)  
#        time_step = min(time_step, 0.01|units.Myr)  







#        else:
            #during run-up towards mass transfer 
#            close_to_RLOF_bool, ratio_rad_rlof, time_step_close_to_mt = self.close_to_RLOF()     
#            if close_to_RLOF_bool:
#                if ratio_rad_rlof >= 1.0:
#                    print('ratio_rad_rlof >1?')
#                    exit(0)
#                else:
#                    print('close to rlof', time_step_close_to_mt, time_step, ratio_rad_rlof)
#                    exit(0)
#                    time_step = min(time_step, time_step_close_to_mt)


      # in case secular code finds RLOF 
        if self.fixed_timestep > 0.|units.yr:
            if REPORT_DT or REPORT_DEBUG:
                print('time:', self.fixed_timestep)
                print('donor time:',  self.determine_time_step_stable_mt())

            t_donor = self.determine_time_step_stable_mt()*time_step_factor_stable_mt # extra small for safety
            t_donor = max(minimum_time_step, min(time_step, t_donor))#
            # although fixed_timestep < time_step, fixed_timestep can be > time_step_stable_mt
            
            if self.fixed_timestep < t_donor:
                time_step = t_donor                
                self.secular_code.parameters.check_for_inner_RLOF = False
            else: 
                time_step = self.fixed_timestep
            self.fixed_timestep = -1|units.Myr
            return time_step # as previous timestep was effectively zero -> min increase is zero           

           


        return time_step
    #-------

    #-------
    #SN functions
    #based on SeBa (Portegies Zwart et al. 1996, Toonen et al. 2012)
    
    def anomaly_converter(self, e_a, ecc, m_a):
    #    mean_anomaly = eccentric_anomaly - ecc*sin(eccentric_anomaly)
        return e_a - ecc * np.sin(e_a) - m_a

    #Paczynski 1990, 348, 485
    def kick_velocity_paczynski(self):
        return self.random_paczynski_velocity(270|units.km/units.s)

    def kick_velocity_SeBa_std(self):
        return self.random_paczynski_velocity(600|units.km/units.s)

    #Paczynski 1990, 348, 485
    #taken from SeBa (Portegies Zwart et al. 1996, Toonen et al. 2012)
    def random_paczynski_velocity(self, v_disp):
    #Velocity distribution used by Paczynski, B., 1990, ApJ 348, 485.
    #with a dispersion of 270 km/s.
    #Phinney uses the same distribution with a dispersion of 600 km/s.
    #The distribution:
    #P(u)du = \frac{4}{\pi} \cdot \frac{du}{(1+u^2)^2},
    #u=v/v_d,
    #v_d is the dispersion velocity = 270 or 600 km/s respectively.

        max_distr = 4./np.pi
        v_max = 4
        prob = max_distr*np.random.uniform(0., 1.)
        velocity = v_max*np.random.uniform(0., 1.)
        dist_value = max_distr/(1+velocity**2.)**2.
    
        while(prob > dist_value):
            prob = max_distr*np.random.uniform(0., 1.)
            velocity = v_max*np.random.uniform(0., 1.)
            dist_value = max_distr/(1+velocity**2.)**2.
        
        vector = self.random_direction()
        return velocity*v_disp*vector

    #Hansen & Phinney 1997, 291, 569
    # V_3d = 300, V_1d = 190 km/s
    def kick_velocity_hansen(self):
        vector = self.random_direction()
        return vector*maxwell.rvs(loc=0, scale = 190)|units.km/units.s
    
    #Hobbs, Lorimer, Lyne & Kramer, 2005, 360, 974
    # V_3d = 400, V_1d = 265 km/s
    def kick_velocity_hobbs(self):
        vector = self.random_direction()
        return vector*maxwell.rvs(loc=0, scale = 265)|units.km/units.s
     
     #Arzoumanian ea 2002, 568, 289
    def kick_velocity_arzoumanian(self):
        vector = self.random_direction()
        if np.random.uniform(0,1) < 0.4:
            return vector*maxwell.rvs(loc=0, scale = 90)|units.km/units.s
        else:
            return vector*maxwell.rvs(loc=0, scale = 500)|units.km/units.s


     # Verbunt, Igoshev & Cator, 2017, 608, 57 - combination of two maxwellians
    def kick_velocity_verbunt(self):
        vector = self.random_direction()
        if np.random.uniform(0,1) < 0.42:
            return vector*maxwell.rvs(loc=0, scale = 75)|units.km/units.s
        else:
            return vector*maxwell.rvs(loc=0, scale = 316)|units.km/units.s


    
    def random_direction(self):
        theta_kick = np.arccos(np.random.uniform(-1, 1))
        phi_kick   = np.random.uniform(0, 2*np.pi)
        x_kick = np.sin(theta_kick)*np.cos(phi_kick)
        y_kick = np.sin(theta_kick)*np.sin(phi_kick)
        z_kick = np.cos(theta_kick)
        return np.array([x_kick, y_kick, z_kick])
	    

    def get_SN_kick(self, star):
        v_kick = [0.,0.,0.]|units.kms
        if star.stellar_type != star.previous_stellar_type and star.stellar_type in stellar_types_SN_remnants:
            if self.SN_kick_distr == 0: 
                v_kick =  [0.,0.,0.]|units.kms
            elif self.SN_kick_distr in [1,2]: # Hobbs, Lorimer, Lyne & Kramer, 2005, 360, 974
                v_kick = self.kick_velocity_hobbs()
            elif self.SN_kick_distr in [3,4]: #Arzoumanian ea 2002, 568, 289
                v_kick = self.kick_velocity_arzoumanian()
            elif self.SN_kick_distr in [5,6]: #Hansen & Phinney 1997, 291, 569
                v_kick = self.kick_velocity_hansen()
            elif self.SN_kick_distr in [7,8]: # Paczynski 1990, 348, 485
                v_kick = self.kick_velocity_paczynski()
            elif self.SN_kick_distr in [9,10]: # Verbunt, Igoshev & Cator, 2017, 608, 57
                v_kick = self.kick_velocity_verbunt()
            else:
                print("SN kick distribution not specified, assuming no kick")
                v_kick =  [0.,0.,0.]|units.kms

            #reduce kick of BH by conservation of momentum
            if self.SN_kick_distr in [2,4,6,8,10] and star.stellar_type == 14|units.stellar_type:
                    v_kick *= (kanonical_neutron_star_mass / star.mass)
        return v_kick

    def save_mean_anomalies_at_SN(self, inner_mean_anomaly, outer_mean_anomaly, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_SN_remnants:
                stellar_system.inner_mean_anomaly = inner_mean_anomaly
                stellar_system.outer_mean_anomaly = outer_mean_anomaly                
        else:
            self.save_mean_anomalies_at_SN(inner_mean_anomaly, outer_mean_anomaly,stellar_system.child1)        
            self.save_mean_anomalies_at_SN(inner_mean_anomaly, outer_mean_anomaly,stellar_system.child2)

    def adjust_spin_after_supernova(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_SN_remnants:
                stellar_system.spin_angular_frequency *= stellar_system.previous_moment_of_inertia_of_star/stellar_system.moment_of_inertia_of_star
        else:
            self.adjust_spin_after_supernova(stellar_system.child1)        
            self.adjust_spin_after_supernova(stellar_system.child2)
 
        
    def adjust_system_after_supernova_kick(self):
        if self.triple.is_star:
            #only velocity kick needs to be taken into account here
            #not followed currently
            return False
        elif self.secular_code.parameters.ignore_tertiary == True:
            #SN kick in binary
            #not implemented currently
            print("Supernova in binary at time = ",self.triple.time) 
            exit(1)                   
        elif not self.is_triple():
            print('SN only implemented in triple')
            exit(1)
                    
           
        #SN in triple
        if self.triple.child1.is_star:
            star = self.triple.child1
            bin = self.triple.child2
        else:
            star = self.triple.child2
            bin = self.triple.child1

        #needs BH mass    
        vel_kick1 = self.get_SN_kick(bin.child1)
        vel_kick2 = self.get_SN_kick(bin.child2)
        vel_kick3 = self.get_SN_kick(star)
    
        #convention in secular code: dm > 0
        dm1 = bin.child1.previous_mass - bin.child1.mass
        dm2 = bin.child2.previous_mass - bin.child2.mass
        dm3 = star.previous_mass - star.mass
 
        #secular code requires the mass of the stars to be the mass prior to the SN
        bin.child1.mass = bin.child1.previous_mass 
        bin.child2.mass = bin.child2.previous_mass
        star.mass = star.previous_mass 
 
        self.save_snapshot()                    
 
 
        # determine stellar anomaly
        # mean anomaly increases uniformly from 0 to 2\pi radians during each orbit
        inner_ecc = bin.eccentricity

        inner_mean_anomaly = np.random.uniform(0., 2.*np.pi)  
        inner_eccentric_anomaly = optimize.brentq(self.anomaly_converter, 0., 2.*np.pi, args=(inner_ecc, inner_mean_anomaly))
        inner_true_anomaly = 2.* np.arctan2(np.sqrt(1-inner_ecc) * np.cos(inner_eccentric_anomaly/2.), np.sqrt(1+inner_ecc) * np.sin(inner_eccentric_anomaly/2.))
        # min teken verschil met cos_true_anomaly = (np.cos(eccentric_anomaly)-ecc) / (1-ecc*np.cos(eccentric_anomaly))
        
        

        outer_ecc = self.triple.eccentricity
        outer_mean_anomaly = np.random.uniform(0., 2.*np.pi)  
        outer_eccentric_anomaly = optimize.brentq(self.anomaly_converter, 0., 2.*np.pi, args=(outer_ecc, outer_mean_anomaly))
        outer_true_anomaly = 2.* np.arctan2(np.sqrt(1-outer_ecc) * np.cos(outer_eccentric_anomaly/2.), np.sqrt(1+outer_ecc) * np.sin(outer_eccentric_anomaly/2.))


        if REPORT_SN_EVOLUTION:
            print('before SN')
            print('\n\ntime:', self.triple.time, self.previous_time, self.triple.time-self.previous_time)
            print('eccentricity:', inner_ecc, outer_ecc)
            print('semi-major axis:', bin.semimajor_axis, self.triple.semimajor_axis)
            print('inner anomaly:', inner_mean_anomaly, inner_eccentric_anomaly, inner_true_anomaly)
            print('outer anomaly:', outer_mean_anomaly, outer_eccentric_anomaly, outer_true_anomaly)
            print('vel_kicks:', vel_kick1, vel_kick2, vel_kick3  )
            print('dm:', dm1, dm2, dm3)        

        self.channel_to_secular.copy()
        V1, V2, V3 = self.secular_code.compute_effect_of_SN_on_triple(vel_kick1, vel_kick2, vel_kick3, dm1, dm2, dm3, inner_true_anomaly, outer_true_anomaly)
        self.channel_from_secular.copy() 
        
        #as secular code required the mass of the stars to be the mass prior to the SN
        #reset mass to after SN
        bin.child1.mass = bin.child1.previous_mass - dm1
        bin.child2.mass = bin.child2.previous_mass - dm2
        star.mass = star.previous_mass - dm3



        if REPORT_SN_EVOLUTION:
            print('after SN')
            print('eccentricity:', bin.eccentricity, self.triple.eccentricity)
            print('semi-major axis:', bin.semimajor_axis, self.triple.semimajor_axis)


        if bin.eccentricity >= 1.0 or bin.eccentricity < 0.0 or bin.semimajor_axis <=0.0|units.RSun or isnan(bin.semimajor_axis.value_in(units.RSun)):
            if REPORT_SN_EVOLUTION:
                print("Inner orbit dissociated by SN at time = ",self.triple.time)                 
            bin.bin_type = bin_type['disintegrated']   
            return False
        elif self.triple.eccentricity >= 1.0 or self.triple.eccentricity < 0.0 or self.triple.semimajor_axis <=0.0|units.RSun or isnan(self.triple.semimajor_axis.value_in(units.RSun)):
            if REPORT_SN_EVOLUTION:
                print("Outer orbit dissociated by SN at time = ",self.triple.time) 
            self.triple.bin_type = bin_type['disintegrated']   

            self.triple.semimajor_axis = 1e100|units.RSun
            self.triple.eccentricity = 0

            self.secular_code.parameters.ignore_tertiary = True
            self.secular_code.parameters.check_for_dynamical_stability = False
            self.secular_code.parameters.check_for_outer_collision = False
            self.secular_code.parameters.check_for_outer_RLOF = False

            # The evolution could be followed further by uncommenting the next line
            # For now we stop the simulation if we do not have a bound triple system anymore
            return False 

        self.adjust_spin_after_supernova() 


        if self.stop_at_SN:
            if REPORT_SN_EVOLUTION:
                print("Supernova at time = ",self.triple.time )
                print(self.triple.child2.child1.mass, self.triple.child2.child1.stellar_type)
                print(self.triple.child2.child2.mass, self.triple.child2.child2.stellar_type)
                print(self.triple.child1.mass, self.triple.child1.stellar_type)
            return False            

        #check for rlof -> if rlof than stop?
#        self.check_RLOF()


        return True
        
    #-------

    #-------
    #evolution
    def resolve_stellar_interaction(self, stellar_system = None):
    # the most inner binary is calculated first, and then move outwards

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if REPORT_TRIPLE_EVOLUTION:
                print('single stellar evolution')

            print('for now no single stellar evolution - exiting program')
            exit(2)
                
            return
        elif self.is_binary(stellar_system):
            if REPORT_TRIPLE_EVOLUTION:
                print('\n perform stellar interaction: binary - double star')
#            stellar_system = perform_stellar_interaction(stellar_system, self)
            stopping_condition = perform_stellar_interaction(stellar_system, self)
            return stopping_condition #stellar interaction
        else:
            if not stellar_system.child1.is_star: #child1 is a multiple
                stopping_condition = self.resolve_stellar_interaction(stellar_system.child1)  
                if not stopping_condition: #stellar interaction
                    return False                                     
            elif not stellar_system.child2.is_star: #child2 is a multiple
                stopping_condition = self.resolve_stellar_interaction(stellar_system.child2)        
                if not stopping_condition: #stellar interaction
                    return False                                     
            else:
                print('resolve_stellar_interaction: structure stellar system unknown')
                print('both children are binaries')
                exit(2)
            
            if REPORT_TRIPLE_EVOLUTION:
                print('\n perform stellar interaction: binary')
#            stellar_system = perform_stellar_interaction(stellar_system, self)            
            stopping_condition = perform_stellar_interaction(stellar_system, self)            
            if not stopping_condition: #stellar interaction
                return False                                     
            return True                  
                                                                
    def determine_mass_transfer_timescale(self, stellar_system = None):

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if REPORT_TRIPLE_EVOLUTION:
                print('single stellar evolution')
            return
        elif self.is_binary(stellar_system):
            if REPORT_TRIPLE_EVOLUTION:
                print('\n determine_mass_transfer_timescale: binary - double star')
            mass_transfer_stability(stellar_system, self)
            if REPORT_TRIPLE_EVOLUTION:
                print('mt rate double star:', stellar_system.mass_transfer_rate)
        else:
            self.determine_mass_transfer_timescale(stellar_system.child1)        
            self.determine_mass_transfer_timescale(stellar_system.child2)        

            if REPORT_TRIPLE_EVOLUTION:
                print('\n determine_mass_transfer_timescale: binary')
            mass_transfer_stability(stellar_system, self)            
            if REPORT_TRIPLE_EVOLUTION:
                print('mt rate binary:', stellar_system.mass_transfer_rate)




    def solve_for_partial_time_step_stable_mass_transfer(self):
        full_dt = self.triple.time - self.previous_time
        self.secular_code.evolve_model(self.previous_time + full_dt * self.triple.child2.part_dt_mt)

        self.channel_from_secular.copy()
        self.check_RLOF() 
        if self.has_donor():
            print('After partial timestep the system should be detached...')
            print(self.has_donor())
            print(self.triple.child1.mass, self.triple.child2.child1.mass, self.triple.child2.child2.mass)
            return False

        self.secular_code.model_time = self.triple.time # not redundant!
        self.triple.child2.part_dt_mt = 1.                    
        return True                      
                    
    
    def safety_check_time_step(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        successfull_dr = True
        
        if stellar_system.is_star:
            dm = (stellar_system.previous_mass - stellar_system.mass)/stellar_system.mass
#            if REPORT_DEBUG:
#                print('wind mass loss rate:', stellar_system.wind_mass_loss_rate, stellar_system.mass,stellar_system.core_mass,)
#                print('relative wind mass losses:', dm)
                        
            if (abs(dm) > error_dm) and not (stellar_system.stellar_type != stellar_system.previous_stellar_type ) and not (stellar_system.stellar_type in stellar_types_SN_remnants):
                successfull_dr = False 
                print('should this statement be here? does this give problems during mt?')

            if self.secular_code.parameters.include_inner_tidal_terms or self.secular_code.parameters.include_outer_tidal_terms:
                dr = (stellar_system.radius - stellar_system.previous_radius)/stellar_system.radius
    
#                if REPORT_DEBUG:    
#                    print('change in radius over time:',  stellar_system.time_derivative_of_radius,)
#                    print('relative change in radius:', dr)
                                     

                if (abs(dr) > error_dr) and not (stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_dr) and  not stellar_system.parent.bin_type == 'stable_mass_transfer':
                    #during stable mass transfer we reset the system to synchronization, 
                    #so tides and spin change due to moment of inertia change are not important, 
                    #therefore the radius change is not important

                    successfull_dr = False                    
            return successfull_dr, 1, stellar_system

        else:
            successfull_dr1, nr_dr1, star_dr1 = self.safety_check_time_step(stellar_system.child1)        
            successfull_dr2, nr_dr2, star_dr2 = self.safety_check_time_step(stellar_system.child2)
            if successfull_dr1==False and successfull_dr2==False:
                if nr_dr1 > 1:
                    return False, 3, np.array([star_dr1[0], star_dr1[1], star_dr2])
                elif nr_dr2 > 1:
                    return False, 3, np.array([star_dr1, star_dr2[0], star_dr2[1]])
                else:
                    return False, 2, np.array([star_dr1, star_dr2])
            elif successfull_dr1==False: 
                return False, nr_dr1, star_dr1
            elif successfull_dr2==False: 
                return False, nr_dr2, star_dr2
            else:
                return True, 0, 0    
            
          
    def recall_memory_one_step_stellar(self, nr_unsuccessfull, star_unsuccessfull):          
        if REPORT_DEBUG:
            print('recall_memory_one_step_stellar')

        self.update_time_derivative_of_radius()
        dt = self.triple.time-self.previous_time
        dt_new = max(minimum_time_step, min(self.determine_time_step(), 0.9*dt))
        self.stellar_code.particles.recall_memory_one_step()
        self.triple.time += (dt_new - dt)                     
        
        self.stellar_code.evolve_model(self.triple.time)
        
#                    self.channel_from_stellar.copy()
        self.channel_from_stellar.copy_attributes(["age", "mass", "core_mass", "radius", "core_radius", "convective_envelope_radius",  "convective_envelope_mass", "stellar_type", "luminosity", "wind_mass_loss_rate",  "temperature"])  #"gyration_radius_sq"                          
#        self.triple.child2.child1.radius = 3.75|units.RSun 
        self.update_stellar_parameters()          
                
        if nr_unsuccessfull > 1:
            for i_rec in range(nr_unsuccessfull):
                triple_number = self.triple.number
                if dt_new <= minimum_time_step:
                    triple_number += 0.5
                dr = (star_unsuccessfull[i_rec].radius - star_unsuccessfull[i_rec].previous_radius)/star_unsuccessfull[i_rec].radius        
        else:
            dr = (star_unsuccessfull.radius - star_unsuccessfull.previous_radius)/star_unsuccessfull.radius        
            triple_number = self.triple.number
            if dt_new <= minimum_time_step:
                triple_number += 0.5

        if dt_new <= minimum_time_step:
            return True, 0, 0    

        return self.safety_check_time_step()                        
        
    #when the stellar code finds RLOF
    def rewind_to_begin_of_rlof_stellar(self, dt):          
        if REPORT_DEBUG:
            print('rewind to begin of rlof stellar')
        dt_new = dt
        dt_min = max(minimum_time_step, self.determine_time_step_stable_mt())
        while self.has_donor() == True and dt_new > dt_min:
            dt_old = self.triple.time - self.previous_time
            dt_new = max(minimum_time_step, 0.5*dt_old)

            self.stellar_code.particles.recall_memory_one_step()
            self.triple.time += (dt_new - dt_old)   
            if REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
                print('\n\ntime:', self.triple.time, self.has_donor())
            self.stellar_code.evolve_model(self.triple.time)
            self.channel_from_stellar.copy_attributes(["age", "mass", "core_mass", "radius", "core_radius", "convective_envelope_radius",  "convective_envelope_mass", "stellar_type", "luminosity", "wind_mass_loss_rate",  "temperature"])  #"gyration_radius_sq"                          
            self.update_stellar_parameters()              
            self.check_RLOF()                    
                    

    #when the secular code finds RLOF
    def rewind_to_begin_of_rlof_secular(self):          
        if REPORT_DEBUG:
            print('rewind to begin of rlof secular')

        #find beginning of RLOF
        self.fixed_timestep = self.secular_code.model_time - self.previous_time
        self.triple.time = self.previous_time
        self.secular_code.model_time = self.previous_time                       
        if self.fixed_timestep < 0.|units.yr:
            print('fixed_timestep < 0: should not be possible', self.triple.time, self.secular_code.model_time, self.fixed_timestep)
            exit(1)                    
        #rewind system
        self.stellar_code.particles.recall_memory_one_step()
        self.refresh_memory() 
                                   
            
    def check_stopping_conditions_stellar_interaction(self):             
        if self.stop_at_merger and self.has_merger():
            if REPORT_TRIPLE_EVOLUTION:
                print('Merger at time = ',self.triple.time )                             
            return False    
        if self.stop_at_disintegrated and self.has_disintegrated():
            if REPORT_TRIPLE_EVOLUTION:
                print('Disintegration of system at time = ',self.triple.time )             
            return False
        return True            
            
            
#    def check_stopping_conditions_stellar(self):  
#        if self.stop_at_outer_mass_transfer and self.has_tertiary_donor():
#            if self.secular_code.model_time < self.triple.time:
#                self.triple.time = self.secular_code.model_time
#
#            if REPORT_TRIPLE_EVOLUTION:
#                print('Mass transfer in outer binary of triple at time = ',self.triple.time)
#            self.triple.bin_type = bin_type['rlof']
#            
#            self.determine_mass_transfer_timescale() # to set the stability #obsolete?
#            return False                                   
#
#        elif self.stop_at_mass_transfer and self.has_donor():
#            if self.secular_code.model_time < self.triple.time:
#                self.triple.time = self.secular_code.model_time
#        
#            if REPORT_TRIPLE_EVOLUTION:
#                print('Mass transfer at time = ',self.triple.time)
#            if self.is_binary(self.triple.child2):
#                self.triple.child2.bin_type = bin_type['rlof'] 
#            elif self.is_binary(self.triple.child1):
#                self.triple.child1.bin_type = bin_type['rlof']    
#            else:
#                print('currently not implemented')
#                exit(-1)    
#                
#            self.determine_mass_transfer_timescale() # to set the stability #obsolete? 
#            return False
#            
#        return True



    def check_stopping_conditions_stellar(self, stellar_system = None):  
        #before check stopping conditions_stellar, always make sure the stability labels are up to date
        #by running self.determine_mass_transfer_timescale()

        if stellar_system == None:
            stellar_system = self.triple

        if self.is_triple(stellar_system):
            if self.stop_at_outer_mass_transfer and self.has_tertiary_donor():
                if self.secular_code.model_time < self.triple.time:
                    self.triple.time = self.secular_code.model_time
    
                if REPORT_TRIPLE_EVOLUTION:
                    print('Mass transfer in outer binary of triple at time = ',self.triple.time)
                self.triple.bin_type = bin_type['rlof']
                return False                                   
            else:
                if not self.check_stopping_conditions_stellar(stellar_system.child1):
                    return False
                if not self.check_stopping_conditions_stellar(stellar_system.child2):
                    return False
                return True
        elif stellar_system.is_star:
            return True
        elif self.is_binary(stellar_system):
            if (self.has_donor() and (self.stop_at_mass_transfer or
                (self.stop_at_stable_mass_transfer and stellar_system.is_mt_stable) or
                (self.stop_at_unstable_mass_transfer and not stellar_system.is_mt_stable) or
                (self.stop_at_eccentric_stable_mass_transfer and stellar_system.is_mt_stable and stellar_system.eccentricity > minimum_eccentricity*5.) or
                (self.stop_at_eccentric_unstable_mass_transfer and not stellar_system.is_mt_stable and stellar_system.eccentricity > minimum_eccentricity*5.)   )):

                if self.secular_code.model_time < self.triple.time:
                    self.triple.time = self.secular_code.model_time
            
                if REPORT_TRIPLE_EVOLUTION:
                    print('Mass transfer in inner binary at time = ',self.triple.time)
                    print(self.stop_at_mass_transfer,self.stop_at_stable_mass_transfer, self.stop_at_unstable_mass_transfer, self.stop_at_eccentric_stable_mass_transfer, self.stop_at_eccentric_unstable_mass_transfer, stellar_system.is_mt_stable)
                if self.is_binary(self.triple.child2):
                    self.triple.child2.bin_type = bin_type['rlof'] 
                elif self.is_binary(self.triple.child1):
                    self.triple.child1.bin_type = bin_type['rlof']    
                else:
                    print('currently not implemented')
                    exit(-1)                        
                return False
            
        return True
        
            
    def check_stopping_conditions(self):
        if self.check_stopping_conditions_stellar()==False:
            return False
        
        if self.stop_at_dynamical_instability and self.triple.dynamical_instability == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time
            self.triple.dynamical_instability = True   
            if REPORT_TRIPLE_EVOLUTION:
                print("Dynamical instability at time = ",self.triple.time)
            return False
        if self.stop_at_inner_collision and self.triple.inner_collision == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time
            self.triple.child2.bin_type = bin_type['collision']
            if REPORT_TRIPLE_EVOLUTION:
                print("Inner collision at time = ",self.triple.time)
            return False
        if self.stop_at_outer_collision and self.triple.outer_collision == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time
            self.triple.bin_type = bin_type['collision']
            if REPORT_TRIPLE_EVOLUTION:
                print("Outer collision at time = ",self.triple.time)
            return False
        if self.stop_at_semisecular_regime and self.triple.semisecular_regime == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time
            self.triple.bin_type = bin_type['semisecular']
            if REPORT_TRIPLE_EVOLUTION:
                print("Semisecular regime at time = ",self.triple.time)
            return False
        if self.triple.time - self.secular_code.model_time > numerical_error|units.Myr:
            print('triple time > sec time: should not be possible', self.triple.time, self.secular_code.model_time)
            print(self.has_donor(), self.secular_code.triples[0].error_flag_secular, self.secular_code.triples[0].dynamical_instability)
            return False
            
        return True


    def check_spin_angular_frequency(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.spin_angular_frequency < 0.0|1./units.Myr: 
                if REPORT_DEBUG:
                  print("Error in spin at time = ",self.triple.time)
               

                return False
            return True                
        else:
            if self.check_spin_angular_frequency(stellar_system.child1) == False or self.check_spin_angular_frequency(stellar_system.child2)==False:
                return False
        return True

    def check_error_flag_secular(self):    
        if self.triple.error_flag_secular < 0:
            if REPORT_TRIPLE_EVOLUTION:
                print("Error in secular code at time = ",self.triple.time)
                print(self.triple.error_flag_secular)
            return False
        return True                
            
#-------------------------    

    def evolve_model(self):
        start_time = time.time()
    
        if REPORT_DEBUG:
            # for plotting data
            times_array = quantities.AdaptingVectorQuantity() 
            a_in_array = quantities.AdaptingVectorQuantity()
            e_in_array = []
            g_in_array = []
            o_in_array = []
            a_out_array = quantities.AdaptingVectorQuantity()
            e_out_array = []
            g_out_array = []
            o_out_array = []
            i_relative_array = []
            m1_array = quantities.AdaptingVectorQuantity()
            m2_array = quantities.AdaptingVectorQuantity()
            m3_array = quantities.AdaptingVectorQuantity()
            spin1_array = []
            spin2_array = []
            spin3_array = []
            r1_array = []
            r2_array = []
            r3_array = []
            moi1_array = []
            moi2_array = []
            moi3_array = []
        
            times_array.append(self.triple.time)
            e_in_array.append(self.triple.child2.eccentricity)
            a_in_array.append(self.triple.child2.semimajor_axis)
            g_in_array.append(self.triple.child2.argument_of_pericenter) 
            o_in_array.append(self.triple.child2.longitude_of_ascending_node)               
            e_out_array.append(self.triple.eccentricity)
            a_out_array.append(self.triple.semimajor_axis)
            g_out_array.append(self.triple.argument_of_pericenter) 
            o_out_array.append(self.triple.longitude_of_ascending_node)               
            i_relative_array.append(self.triple.relative_inclination)        
            m1_array.append(self.triple.child2.child1.mass)
            m2_array.append(self.triple.child2.child2.mass)
            m3_array.append(self.triple.child1.mass)
            spin1_array.append(self.triple.child2.child1.spin_angular_frequency.value_in(1./units.Myr))
            spin2_array.append(self.triple.child2.child2.spin_angular_frequency.value_in(1./units.Myr))
            spin3_array.append(self.triple.child1.spin_angular_frequency.value_in(1./units.Myr))
            r1_array.append(self.triple.child2.child1.radius.value_in(units.RSun))
            r2_array.append(self.triple.child2.child2.radius.value_in(units.RSun))
            r3_array.append(self.triple.child1.radius.value_in(units.RSun))
            moi1_array.append(self.triple.child2.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
            moi2_array.append(self.triple.child2.child2.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
            moi3_array.append(self.triple.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
        
        if REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
            print('kozai timescale:', self.kozai_timescale(), self.triple.kozai_type, self.tend)

        self.determine_mass_transfer_timescale()
        self.save_snapshot()  
        while self.triple.time<self.tend: 
		
	    # if the maximum CPU time has been exceeded for the system, stop the evolution and continue with the next system
            end_time = time.time()
            self.triple.CPU_time = end_time - start_time
            if (self.stop_at_CPU_time == True) and float(end_time - start_time) > self.max_CPU_time:
                print('stopping conditions maximum CPU time')
                print("CPU time: ", self.triple.CPU_time)
                break
		
            if REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
                print('\n\n kozai timescale:', self.kozai_timescale(), self.triple.kozai_type, self.octupole_parameter())
                        
            if self.has_stellar_type_changed() or self.has_kozai_type_changed():
                self.save_snapshot()  

            dt = self.determine_time_step()  
            if not no_stellar_evolution: 
                self.update_previous_stellar_parameters()
                self.stellar_code.particles.refresh_memory()
                
                
            self.triple.time += dt   
            self.previous_dt = dt         
            if REPORT_DEBUG or REPORT_DT:
                print('\t time:', self.triple.time, self.previous_time, self.previous_dt)
                            
            #do stellar evolution 
            if not no_stellar_evolution: 
                if REPORT_DEBUG:
                    print('Stellar evolution')

                self.stellar_code.evolve_model(self.triple.time)
                self.channel_from_stellar.copy_attributes(["age", "mass", "core_mass", "radius", "core_radius", "convective_envelope_radius",  "convective_envelope_mass", "stellar_type", "luminosity", "wind_mass_loss_rate", "temperature"]) #"gyration_radius_sq"  
                self.update_stellar_parameters()     
         
                successfull_step, nr_unsuccessfull, star_unsuccessfull = self.safety_check_time_step() 
                while successfull_step == False:
                    successfull_step, nr_unsuccessfull, star_unsuccessfull = self.recall_memory_one_step_stellar(nr_unsuccessfull, star_unsuccessfull)

                # if SN has taken place
                if self.has_stellar_type_changed_into_SN_remnant():
                    if REPORT_TRIPLE_EVOLUTION:
                        print("Supernova at time = ", self.triple.time, dt)
                    if self.adjust_system_after_supernova_kick()==False:
                        break  
                    self.instantaneous_evolution = True #skip secular evolution

                self.check_RLOF()                                       
                #find beginning of RLOF
                if self.has_donor() and self.triple.bin_type == 'detached' and self.triple.child2.bin_type == 'detached' and dt > minimum_time_step:
#                    self.rewind_to_begin_of_rlof_stellar(dt)
#                    print('RLOF:', self.triple.child2.child1.is_donor, self.triple.bin_type , self.triple.child2.bin_type )

                    self.stellar_code.particles.recall_memory_one_step()
                    self.channel_from_stellar.copy_attributes(["age", "mass", "core_mass", "radius", "core_radius", "convective_envelope_radius",  "convective_envelope_mass", "stellar_type", "luminosity", "wind_mass_loss_rate",  "temperature"])  #"gyration_radius_sq"                          
#                    self.triple.child2.child1.radius = 3.75|units.RSun 
                    self.update_stellar_parameters()                              
                    self.refresh_memory()                     
                    #note that 'previous' parameters cannot be reset 
                    #resetting is_donor in determine_time_step                                    
                    continue

                                         

            #needed for nucleair timescale 
            self.update_time_derivative_of_radius() 
          
            #do stellar interaction
            if REPORT_DEBUG:
                print('Stellar interaction')

            self.determine_mass_transfer_timescale()
            if self.check_stopping_conditions_stellar()==False:
                print('stopping conditions stellar')
                break  
            if not self.resolve_stellar_interaction():
                print('stopping conditions stellar interaction')
                break            
            self.update_time_derivative_of_radius() # includes radius change from wind and ce, not stable mt
            self.update_time_derivative_of_moment_of_inertia() # includes mass and radius change from wind and mass transfer            
            if self.check_stopping_conditions_stellar_interaction()==False: 
                print('stopping conditions stellar interaction 2')
                break
            
            
                                                           
            # do secular evolution
            if self.instantaneous_evolution == False:
                #needed for refreshing memory in case secular finds RLOF    
                previous_semimajor_axis_in = self.triple.child2.semimajor_axis
                previous_eccentricity_in = self.triple.child2.eccentricity
                previous_argument_of_pericenter_in = self.triple.child2.argument_of_pericenter
                previous_longitude_of_ascending_node_in = self.triple.child2.longitude_of_ascending_node
                previous_spin1 = self.triple.child2.child1.spin_angular_frequency
                previous_spin2 = self.triple.child2.child2.spin_angular_frequency
    
                if REPORT_DEBUG:
                    print('Secular evolution')

                self.channel_to_secular.copy()
                
                # if mass transfer should only take place for a fraction of the timestep
                # e.g. at the end of mass transfer when the envelope is thin
                if self.triple.child2.part_dt_mt < 1: # inner binary, see function determine_partial_time_step_stable_mass_transfer
                    print('skipping')
                    if self.solve_for_partial_time_step_stable_mass_transfer() == False:
                        break
                else:
                    self.secular_code.evolve_model(self.triple.time)
                    
                if REPORT_DEBUG:
                    print('Secular evolution finished')
                    
                if self.triple.time - self.secular_code.model_time < -1*numerical_error|units.Myr and self.secular_code.triples[0].error_flag_secular >= 0:
                    print('triple time < sec time: should not be possible', self.triple.time, self.secular_code.model_time)
                    print(self.has_donor(), self.secular_code.triples[0].error_flag_secular)
                    break
                    
                elif self.has_donor() and self.triple.bin_type == 'detached' and self.triple.child2.bin_type == 'detached':
                    self.determine_mass_transfer_timescale()            
                    if self.check_stopping_conditions_stellar()==False:
                        print('stopping conditions stellar 2')                    
                        break

                    self.rewind_to_begin_of_rlof_secular()
                    self.triple.child2.semimajor_axis = previous_semimajor_axis_in
                    self.triple.child2.eccentricity = previous_eccentricity_in
                    self.triple.child2.previous_argument_of_pericenter_in = previous_argument_of_pericenter_in
                    self.triple.child2.previous_longitude_of_ascending_node_in = previous_longitude_of_ascending_node_in
                    self.triple.child2.child1.spin_angular_frequency = previous_spin1
                    self.triple.child2.child2.spin_angular_frequency = previous_spin2
                    continue
                                  
                self.channel_from_secular.copy() 



                if self.check_error_flag_secular()==False:
                    break
                if self.check_spin_angular_frequency()==False:
                    break
                self.determine_mass_transfer_timescale()
                if self.check_stopping_conditions()==False:
                    print('stopping conditions')                
                    break                
                if not self.stop_at_inner_collision and self.triple.inner_collision == True:
                    perform_inner_collision(self)
                                        
            else:
                if REPORT_DEBUG:
                    print('skip secular')
                self.secular_code.model_time = self.triple.time
                self.instantaneous_evolution = False
            



            if REPORT_DEBUG:
                # for plotting data
                times_array.append(self.triple.time)
                e_in_array.append(self.triple.child2.eccentricity)
                a_in_array.append(self.triple.child2.semimajor_axis)
                g_in_array.append(self.triple.child2.argument_of_pericenter) 
                o_in_array.append(self.triple.child2.longitude_of_ascending_node)               
                e_out_array.append(self.triple.eccentricity)
                a_out_array.append(self.triple.semimajor_axis)
                g_out_array.append(self.triple.argument_of_pericenter) 
                o_out_array.append(self.triple.longitude_of_ascending_node)               
                i_relative_array.append(self.triple.relative_inclination)        
                m1_array.append(self.triple.child2.child1.mass)
                m2_array.append(self.triple.child2.child2.mass)
                m3_array.append(self.triple.child1.mass)
                spin1_array.append(self.triple.child2.child1.spin_angular_frequency.value_in(1./units.Myr))
                spin2_array.append(self.triple.child2.child2.spin_angular_frequency.value_in(1./units.Myr))
                spin3_array.append(self.triple.child1.spin_angular_frequency.value_in(1./units.Myr))
                r1_array.append(self.triple.child2.child1.radius.value_in(units.RSun))
                r2_array.append(self.triple.child2.child2.radius.value_in(units.RSun))
                r3_array.append(self.triple.child1.radius.value_in(units.RSun))
                moi1_array.append(self.triple.child2.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
                moi2_array.append(self.triple.child2.child2.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
                moi3_array.append(self.triple.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
            
                                    
        self.save_snapshot()        
            
            
        if REPORT_DEBUG:
            # for plotting data
            e_in_array = np.array(e_in_array)
            g_in_array = np.array(g_in_array)
            o_in_array = np.array(o_in_array)
            e_out_array = np.array(e_out_array)
            g_out_array = np.array(g_out_array)
            o_out_array = np.array(o_out_array)
            i_relative_array = np.array(i_relative_array)
            spin1_array = np.array(spin1_array)
            spin2_array = np.array(spin2_array)
            spin3_array = np.array(spin3_array)
            r1_array = np.array(r1_array)
            r2_array = np.array(r2_array)
            r3_array = np.array(r3_array)
            moi1_array = np.array(moi1_array)
            moi2_array = np.array(moi2_array)
            moi3_array = np.array(moi3_array)
    
            self.plot_data = plot_data_container()
            self.plot_data.times_array = times_array
            self.plot_data.a_in_array = a_in_array
            self.plot_data.e_in_array = e_in_array
            self.plot_data.g_in_array = g_in_array
            self.plot_data.o_in_array = o_in_array        
            self.plot_data.a_out_array = a_out_array
            self.plot_data.e_out_array = e_out_array
            self.plot_data.g_out_array = g_out_array
            self.plot_data.o_out_array = o_out_array        
            self.plot_data.i_relative_array = i_relative_array
            self.plot_data.m1_array = m1_array
            self.plot_data.m2_array = m2_array
            self.plot_data.m3_array = m3_array
            self.plot_data.spin1_array = spin1_array
            self.plot_data.spin2_array = spin2_array
            self.plot_data.spin3_array = spin3_array
            self.plot_data.r1_array = r1_array
            self.plot_data.r2_array = r2_array
            self.plot_data.r3_array = r3_array
            self.plot_data.moi1_array = moi1_array
            self.plot_data.moi2_array = moi2_array
            self.plot_data.moi3_array = moi3_array
        
    #-------




class plot_data_container():
    def __init__(self):
        return

def plot_function(triple, dir_plots):
    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
    t_max_Myr = max(times_array_Myr)
    a_in_array_AU = triple.plot_data.a_in_array.value_in(units.AU)
    g_in_array = triple.plot_data.g_in_array
    e_in_array = triple.plot_data.e_in_array
    i_relative_array = triple.plot_data.i_relative_array
    o_in_array = triple.plot_data.o_in_array
    a_out_array_AU = triple.plot_data.a_out_array.value_in(units.AU)
    g_out_array = triple.plot_data.g_out_array
    e_out_array = triple.plot_data.e_out_array
    o_out_array = triple.plot_data.o_out_array
    m1_array = triple.plot_data.m1_array.value_in(units.MSun)
    m2_array = triple.plot_data.m2_array.value_in(units.MSun)
    m3_array = triple.plot_data.m3_array.value_in(units.MSun)
    spin1_array = triple.plot_data.spin1_array
    spin2_array = triple.plot_data.spin2_array
    spin3_array = triple.plot_data.spin3_array
    r1_array = triple.plot_data.r1_array
    r2_array = triple.plot_data.r2_array
    r3_array = triple.plot_data.r3_array
    moi1_array = triple.plot_data.moi1_array
    moi2_array = triple.plot_data.moi2_array
    moi3_array = triple.plot_data.moi3_array
    
    f = open(triple.file_name[:-4]+'.txt','w')
    f.write('#' + str(t_max_Myr) + '\n')
    for i_p in range(len(times_array_Myr)):
        f.write(str(times_array_Myr[i_p]) + '\t')
        f.write(str(m1_array[i_p] ) + '\t')
        f.write(str(m2_array[i_p] ) + '\t')
        f.write(str(m3_array[i_p] ) + '\t')
        f.write(str(spin1_array[i_p] ) + '\t')
        f.write(str(spin2_array[i_p] ) + '\t')
        f.write(str(spin3_array[i_p] ) + '\t')
        f.write(str(r1_array[i_p] ) + '\t')
        f.write(str(r2_array[i_p] ) + '\t')
        f.write(str(r3_array[i_p] ) + '\t')
        f.write(str(moi1_array[i_p] ) + '\t')
        f.write(str(moi2_array[i_p] ) + '\t')
        f.write(str(moi3_array[i_p] ) + '\t')
        f.write(str(i_relative_array[i_p] ) + '\t')
        f.write(str(a_in_array_AU[i_p] ) + '\t')
        f.write(str(g_in_array[i_p] ) + '\t')
        f.write(str(e_in_array[i_p] ) + '\t')
        f.write(str(o_in_array[i_p] ) + '\t')
        f.write(str(a_out_array_AU[i_p] ) + '\t')
        f.write(str(g_out_array[i_p] ) + '\t')
        f.write(str(e_out_array[i_p] ) + '\t')
        f.write(str(o_out_array[i_p] ) + '\t')
        f.write('\n')
    f.close()
    
    for i_s in range(len(times_array_Myr)):
        print(times_array_Myr[i_s], end = ' ')
        print( a_in_array_AU[i_s], end = ' ')
        print( g_in_array[i_s], end = ' ')
        print( e_in_array[i_s], end = ' ')
        print( i_relative_array[i_s], end = ' ')
        print( o_in_array[i_s], end = ' ')
        print( a_out_array_AU[i_s], end = ' ')
        print( g_out_array[i_s], end = ' ')
        print( e_out_array[i_s], end = ' ')
        print( o_out_array[i_s], end = ' ')
        print( m1_array[i_s], end = ' ')
        print( m2_array[i_s], end = ' ')
        print( m3_array[i_s], end = ' ')
        print( spin1_array[i_s], end = ' ')
        print( spin2_array[i_s], end = ' ')  
        print( spin3_array[i_s], end = ' ')  
        print( r1_array[i_s], end = ' ')    
        print( r2_array[i_s], end = ' ')    
        print( r3_array[i_s], end = ' ')    
        print( moi1_array[i_s], end = ' ')   
        print( moi2_array[i_s], end = ' ')   
        print( moi3_array[i_s])   
    
    
            
    ### plots to test secular code ###
    import amuse.plot as aplt
    import matplotlib.pyplot as plt
    
#    generic_name = '_M'+str(m1_array[0]) + '_m'+str(m2_array[0]) +'_n'+str(m3_array[0]) + '_a'+str(a_in_array_AU[0]) + '_A'+str(a_out_array_AU[0]) + '_e'+str(e_in_array[0]) + '_E'+str(e_out_array[0]) + '_i'+str(i_relative_array[0]/np.pi*180.0) + '_g'+str(g_in_array[0]) + '_G'+str(g_out_array[0]) + '_o'+str(o_in_array[0]) + '_O'+str(o_out_array[0]) + '_t'+str(t_max_Myr) + '_maxdr'+str(triple.maximum_radius_change_factor)+'_edr'+str(error_dr)
    generic_name = ''

    figure = plt.figure(figsize=(10,13))
    N_subplots = 4
    
    plot_e = figure.add_subplot(N_subplots,1,1)
    plot_i_relative = figure.add_subplot(N_subplots,1,2)
    plot_a_in = figure.add_subplot(N_subplots,1,3)
    plot_a_out = figure.add_subplot(N_subplots,1,4)
    
    plot_e.plot(times_array_Myr,e_in_array, label= '$e_\mathrm{in}$')
    plot_e.plot(times_array_Myr,e_out_array, label= '$e_\mathrm{out}$')
    plot_e.set_xlim(0,t_max_Myr)
    plot_e.set_xlabel('$t/\mathrm{Myr}$')
    plot_e.set_ylabel('$e$')
    plot_e.legend(loc=0)
    
    plot_i_relative.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative.set_xlim(0,t_max_Myr)
    plot_i_relative.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')
    
    plot_a_in.plot(times_array_Myr,a_in_array_AU)
    plot_a_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_in.set_ylabel('$a_\mathrm{in}$')

    plot_a_out.plot(times_array_Myr,a_out_array_AU)
    plot_a_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_out.set_ylabel('$a_\mathrm{out}$')
    
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)
    plt.savefig(dir_plots+'TRES'+generic_name+'.pdf')
#    plt.show()
    plt.close()








    figure = plt.figure(figsize=(10,13))
    N_subplots = 4

    plot_e_in = figure.add_subplot(N_subplots,1,1)
    plot_i_relative = figure.add_subplot(N_subplots,1,2)
    plot_e_in_g_in = figure.add_subplot(N_subplots,1,3)
    plot_a_in = figure.add_subplot(N_subplots,1,4)


    plot_e_in.plot(times_array_Myr,e_in_array)
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_e_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_in.set_ylabel('$e_\mathrm{in}$')

    plot_i_relative.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative.set_xlim(0,t_max_Myr)
    plot_i_relative.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')

    plot_e_in_g_in.plot(np.cos(g_in_array),e_in_array)
    plot_e_in_g_in.set_xlabel('$\cos(g_\mathrm{in})$')
    plot_e_in_g_in.set_ylabel('$e_\mathrm{in}$')

    plot_a_in.plot(times_array_Myr,a_in_array_AU)
    plot_a_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_in.set_ylabel('$a_\mathrm{in}$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    plt.savefig(dir_plots+'TRES_inner_orbit'+generic_name+'.pdf')
#    plt.show()
    plt.close()


    dyn_inst =  2.8 / (1-e_out_array) * (1-0.3*i_relative_array/np.pi) * ((1+m3_array/(m1_array+m2_array))*(1+e_out_array)/(np.sqrt(1-e_out_array)))**0.4 / (a_out_array_AU/a_in_array_AU)
    oct = (m1_array-m2_array)/(m1_array+m2_array) * a_in_array_AU/a_out_array_AU * e_out_array/(1-e_out_array**2)
    semiseq = (5.*np.pi* m3_array/(m1_array+m2_array) * (a_in_array_AU/a_out_array_AU/(1-e_out_array))**3)/ np.sqrt(1-e_in_array) 

    plt.semilogy(times_array_Myr, dyn_inst, '.')
    plt.semilogy(times_array_Myr, oct, '.')
    plt.semilogy(times_array_Myr, semiseq, '.')
    plt.xlabel('t (Myr)')
    plt.ylabel('stability,oct and semiseq')
    plt.legend(loc=0)
    plt.savefig(dir_plots+'dyn_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    plt.plot(times_array_Myr,e_in_array)
    plt.plot(times_array_Myr,e_in_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$e_\mathrm{in}$')
    plt.savefig(dir_plots+'e_in_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


    plt.plot(times_array_Myr,e_in_array)
    plt.plot(times_array_Myr,e_in_array, '.')
    plt.plot(times_array_Myr,e_out_array)
    plt.plot(times_array_Myr,e_out_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$e_\mathrm{in}$')
    plt.savefig(dir_plots+'e_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()
    
    plt.plot(times_array_Myr,g_in_array)
    plt.plot(times_array_Myr,g_in_array, '.')
    plt.plot(times_array_Myr,g_out_array)
    plt.plot(times_array_Myr,g_out_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$g_\mathrm{in}$')
    plt.savefig(dir_plots+'g_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()
    
    
    plt.plot(e_in_array,(g_in_array%np.pi)/np.pi*180)
    plt.plot(e_in_array,(g_in_array%np.pi)/np.pi*180, '.')
#    plt.plot(e_out_array,g_out_array)
#    plt.plot(e_out_array,g_out_array, '.')
#    plt.xlim(0,1)
    plt.xlabel('$e_\mathrm{in}$')
    plt.ylabel('$g_\mathrm{in}$')
    plt.savefig(dir_plots+'g_e_inner'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    

    plt.plot(times_array_Myr,o_in_array)
    plt.plot(times_array_Myr,o_in_array, '.')
    plt.plot(times_array_Myr,o_out_array)
    plt.plot(times_array_Myr,o_out_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$o_\mathrm{in}$')
    plt.savefig(dir_plots+'o_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


#    a_in_final_theory =  a_in_array_AU[0] * (m1_array[0] * m2_array[0] / m1_array / m2_array)**2 #stable mt
    a_in_final_theory =  a_in_array_AU[0] * (m1_array[0] + m2_array[0]) / (m1_array + m2_array)#wind 
    plt.plot(times_array_Myr,a_in_array_AU)
    plt.plot(times_array_Myr,a_in_array_AU, '.')
    plt.plot(times_array_Myr,a_in_final_theory)
    plt.plot(times_array_Myr,a_in_final_theory, '.')
    
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$a_\mathrm{in}$')
    plt.savefig(dir_plots+'semi_in_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()



    constants = 6283728.92847 # constants.G #in au and solar mass per Myr
    corot_spin_inner = 1./np.sqrt(a_in_array_AU**3/(m1_array+m2_array))*constants
    corot_spin_outer = 1./np.sqrt(a_out_array_AU**3/(m1_array+m2_array+m3_array))*constants
    RSun_in_AU = 0.00464913034382
    critical_rot_1 = constants*np.sqrt(m1_array/(RSun_in_AU*r1_array)**3)
    critical_rot_2 = constants*np.sqrt(m2_array/(RSun_in_AU*r2_array)**3)
    critical_rot_3 = constants*np.sqrt(m3_array/(RSun_in_AU*r3_array)**3)

    plt.plot(times_array_Myr,spin1_array, 'b-')
    plt.plot(times_array_Myr,spin1_array, 'b.')
    plt.plot(times_array_Myr,spin2_array, 'g-')
    plt.plot(times_array_Myr,spin2_array, 'g.')
    plt.plot(times_array_Myr,spin3_array, 'r-')
    plt.plot(times_array_Myr,spin3_array, 'r.')

    plt.plot(times_array_Myr,corot_spin_inner, 'c-')
    plt.plot(times_array_Myr,corot_spin_inner, 'c,')
    plt.plot(times_array_Myr,corot_spin_outer, 'm-')
    plt.plot(times_array_Myr,corot_spin_outer, 'm,')

    plt.plot(times_array_Myr,critical_rot_1, 'y-')
    plt.plot(times_array_Myr,critical_rot_1, 'y,')
    plt.plot(times_array_Myr,critical_rot_2, 'k-')
    plt.plot(times_array_Myr,critical_rot_2, 'k,')
    plt.plot(times_array_Myr,critical_rot_3, 'k-')
    plt.plot(times_array_Myr,critical_rot_3, 'ko')

    
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$spin$')
    plt.savefig(dir_plots+'spin_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()



#        ms = 1|units.MSun
#        rs = 1|units.AU
#        J2 = ms**3 * rs * constants.G
#        J = np.sqrt(J2)
#        print(J)  #2.9071938904e+11 [RSun**2 * MSun * Myr**-1]       

#        rs = 1|units.RSun
#        print(J)  #19822565357. [RSun**2 * MSun * Myr**-1]   

    constants_Jorb = 2.9071938904e+11 #[RSun**2 * MSun * Myr**-1]
    J_orb2 = m1_array**2 * m2_array**2 / (m1_array+m2_array) * a_in_array_AU * ( 1-e_in_array**2)
    J_orb = np.sqrt(J_orb2)*constants_Jorb

    J_spin1 =  spin1_array * moi1_array
    J_spin2 =  spin2_array * moi2_array
    J_spin3 =  spin3_array * moi3_array

    plt.plot(times_array_Myr, J_orb)
    plt.plot(times_array_Myr,J_orb, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J orbit$')
    plt.savefig(dir_plots+'Jorbit_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    plt.plot(times_array_Myr,J_spin1)
    plt.plot(times_array_Myr,J_spin1, '.')
    plt.plot(times_array_Myr,J_spin2)
    plt.plot(times_array_Myr,J_spin2, '.')
    plt.plot(times_array_Myr,J_spin3)
    plt.plot(times_array_Myr,J_spin3, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J spin$')
    plt.savefig(dir_plots+'Jspin_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


    plt.plot(times_array_Myr, J_orb)
    plt.plot(times_array_Myr,J_orb, '.')
    plt.plot(times_array_Myr,J_spin1)
    plt.plot(times_array_Myr,J_spin1, '.')
    plt.plot(times_array_Myr,J_spin2)
    plt.plot(times_array_Myr,J_spin2, '.')
    plt.plot(times_array_Myr,J_spin3)
    plt.plot(times_array_Myr,J_spin3, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J spin$')
    plt.savefig(dir_plots+'Js_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    plt.semilogy(times_array_Myr,moi1_array)
    plt.semilogy(times_array_Myr,moi1_array, '.')
    plt.semilogy(times_array_Myr,moi2_array)
    plt.semilogy(times_array_Myr,moi2_array, '.')
    plt.semilogy(times_array_Myr,moi3_array)
    plt.semilogy(times_array_Myr,moi3_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$moi$')
    plt.savefig(dir_plots+'moment_of_inertia_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    

    plt.semilogy(times_array_Myr,r1_array)
    plt.semilogy(times_array_Myr,r1_array, '.')
    plt.semilogy(times_array_Myr,r2_array)
    plt.semilogy(times_array_Myr,r2_array, '.')
    plt.semilogy(times_array_Myr,r3_array)
    plt.semilogy(times_array_Myr,r3_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$radius$')
    plt.savefig(dir_plots+'radius_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()
    
    dr1_array =r1_array[1:]-r1_array[:-1]
    dr2_array =r2_array[1:]-r2_array[:-1]
    dr3_array =r3_array[1:]-r3_array[:-1]
    dt_array = times_array_Myr[1:] - times_array_Myr[:-1]
    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array)
    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array, 'b.')
    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array*-1., 'b*')
    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array)
    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array, 'g.')
    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array*-1., 'g*')
    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array)
    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array, 'r.')
    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array*-1., 'r*')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$dr/dt$')
    plt.savefig(dir_plots+'drdt_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()




    plt.semilogy(times_array_Myr[1:], dr1_array/r1_array[1:])
    plt.semilogy(times_array_Myr[1:], dr1_array/r1_array[1:], 'b.')
    plt.semilogy(times_array_Myr[1:], -1.*dr1_array/r1_array[1:], 'b*')
    plt.semilogy(times_array_Myr[1:], dr2_array/r2_array[1:])
    plt.semilogy(times_array_Myr[1:], dr2_array/r2_array[1:], 'g.')
    plt.semilogy(times_array_Myr[1:], -1.*dr2_array/r2_array[1:], 'g*')
    plt.semilogy(times_array_Myr[1:], dr3_array/r3_array[1:])
    plt.semilogy(times_array_Myr[1:], dr3_array/r3_array[1:], 'r.')
    plt.semilogy(times_array_Myr[1:], -1.*dr3_array/r3_array[1:], 'r*')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$dr/r$')
    plt.savefig(dir_plots+'dr_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()



    dm1_array =m1_array[:-1]-m1_array[1:]
    dm2_array =m2_array[:-1]-m2_array[1:]
    dm3_array =m3_array[:-1]-m3_array[1:]
    
    plt.plot(times_array_Myr[1:], m1_array[1:])
    plt.plot(times_array_Myr[1:], m1_array[1:], '.')
    plt.plot(times_array_Myr[1:], m2_array[1:])
    plt.plot(times_array_Myr[1:], m2_array[1:], '.')
    plt.plot(times_array_Myr[1:], m3_array[1:])
    plt.plot(times_array_Myr[1:], m3_array[1:], '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$mass$')
    plt.savefig(dir_plots+'mass_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


#    plt.semilogy(times_array_Myr[1:], dm1_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm1_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dm2_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm2_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dm3_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm3_array/dt_array, '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dm/dt$')
#    plt.savefig(dir_plots+'dmdt_time'+generic_name+'.pdf')
#    plt.show()
#    plt.close()

#
#    plt.semilogy(times_array_Myr[1:], dm1_array/m1_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm1_array/m1_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dm2_array/m2_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm2_array/m2_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dm3_array/m3_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm3_array/m3_array[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dm/m$')
#    plt.savefig(dir_plots+'dm_time'+generic_name+'.pdf')
#    plt.show()
#    plt.close()

 
    plt.semilogy(times_array_Myr[1:], dt_array)
    plt.semilogy(times_array_Myr[1:], dt_array, '.')
    plt.semilogy(times_array_Myr[1:], dt_array)
    plt.semilogy(times_array_Myr[1:], dt_array, '.')
    plt.semilogy(times_array_Myr[1:], dt_array)
    plt.semilogy(times_array_Myr[1:], dt_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$dt$')
    plt.savefig(dir_plots+'dt_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()





#   wind a = ai * Mti/Mt
#    Mtot = m1_array+m2_array    
#    plt.plot(times_array_Myr,a_in_array_AU)
#    plt.plot(times_array_Myr,a_in_array_AU, '.')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot, '.')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:])
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{in}$')
#    plt.savefig(dir_plots+'semi_inner_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_inner_rel_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()


#   cons mt a = ai * (m1i*m2i*/m1/m2)**2
#    plt.plot(times_array_Myr,a_in_array_AU, 'b-')
#    plt.plot(times_array_Myr,a_in_array_AU, 'b.')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2, 'g-')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2, 'g.')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*(m1_array[:-1]*m2_array[:-1]/m1_array[1:]/m2_array[1:])**2, 'r-')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*(m1_array[:-1]*m2_array[:-1]/m1_array[1:]/m2_array[1:])**2, 'r.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{in}$')
#    plt.savefig(dir_plots+'semi_inner_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_inner_rel_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()



#    dm = (m1_array[1:] - m1_array[:-1] )
#    dt = (times_array_Myr[1:] - times_array_Myr[:-1])
#    dmdt = (m1_array[1:] - m1_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])

#    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
    
    
    
    
    
#    plt.plot(times_array_Myr, g_in_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$g_\mathrm{in}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(g_in_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(g_\mathrm{in})$')
#    plt.show()
#    plt.close()
#
#
#    plt.plot(times_array_Myr, o_in_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$o_\mathrm{in}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(o_in_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(o_\mathrm{in})$')
#    plt.show()
#    plt.close()


    #outer binary
    figure = plt.figure(figsize=(10,13))
    N_subplots = 4

    plot_e_out = figure.add_subplot(N_subplots,1,1)
    plot_i_relative2 = figure.add_subplot(N_subplots,1,2)
    plot_e_out_g_out = figure.add_subplot(N_subplots,1,3)
    plot_a_out = figure.add_subplot(N_subplots,1,4)

#    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
#    t_max_Myr = max(times_array_Myr)
#    i_relative_array = triple.plot_data.i_relative_array

    plot_e_out.plot(times_array_Myr,e_out_array)
    plot_e_out.set_xlim(0,t_max_Myr)
    plot_e_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_out.set_ylabel('$e_\mathrm{out}$')

    plot_i_relative2.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative2.set_xlim(0,t_max_Myr)
    plot_i_relative2.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative2.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative2.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')

    plot_e_out_g_out.plot(np.cos(g_out_array),e_out_array)
    plot_e_out_g_out.set_xlabel('$\cos(g_\mathrm{out})$')
    plot_e_out_g_out.set_ylabel('$e_\mathrm{out}$')

    plot_a_out.plot(times_array_Myr,a_out_array_AU)
    plot_a_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_out.set_ylabel('$a_\mathrm{out}$')

    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)
    plt.savefig(dir_plots+'TRES_outer_orbit'+generic_name+'.pdf')
#    plt.show()
    plt.close()



#    plt.plot(times_array_Myr,e_out_array)
#    plt.xlim(0,t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$e_\mathrm{out}$')
#    plt.savefig(dir_plots+'e_out_time'+generic_name+'.pdf')
#    plt.show()
#    plt.close()


#    Mtott = m1_array+m2_array+m3_array    
#    plt.plot(times_array_Myr,a_out_array_AU)
#    plt.plot(times_array_Myr,a_out_array_AU, '.')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott, '.')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:])
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{out}$')
#    plt.savefig(dir_plots+'semi_outer_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_outer_rel_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#
#    m_in_array = m1_array+m2_array
#    plt.plot(times_array_Myr,a_out_array_AU, 'b-')
#    plt.plot(times_array_Myr,a_out_array_AU, 'b.')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2, 'g-')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2, 'g.')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*(m_in_array[:-1]*m3_array[:-1]/m_in_array[1:]/m3_array[1:])**2, 'r-')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*(m_in_array[:-1]*m3_array[:-1]/m_in_array[1:]/m3_array[1:])**2, 'r.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{out}$')
#    plt.savefig(dir_plots+'semi_outer_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_outer_rel_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#

#    dm = (m3_array[1:] - m3_array[:-1] )
#    dmdt = (m3_array[1:] - m3_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])

#    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#    plt.close()



#    plt.plot(np.cos(g_out_array), e_in_array)
#    plt.xlabel('$\cos(g_\mathrm{out})$')
#    plt.ylabel('$e_\mathrm{in}$')
#    plt.show()
#    plt.close()
#
#    plt.plot(np.cos(g_in_array), np.cos(g_out_array))
#    plt.xlabel('$\cos(g_\mathrm{in})$')
#    plt.ylabel('$\cos(g_\mathrm{out})$')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, g_out_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$g_\mathrm{out}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(g_out_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(g_\mathrm{out})$')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, o_out_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$o_\mathrm{out}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(o_out_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(o_\mathrm{out})$')
#    plt.show()
#    plt.close()
#
#    aplt.plot(times_array_Myr, m1_array)
#    aplt.plot(times_array_Myr, m1_array, '.')
#    aplt.plot(times_array_Myr, m2_array)
#    aplt.plot(times_array_Myr, m2_array, '.')
#    aplt.plot(times_array_Myr, m3_array)
#    aplt.plot(times_array_Myr, m3_array, '.')
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$M/\mathrm{MSun}$')
#    plt.show()
#    plt.close()
#    
    
#-----
#for running TRES.py from other routines
def main(inner_primary_mass= 1.3|units.MSun, inner_secondary_mass= 0.5|units.MSun, outer_mass= 0.5|units.MSun,
            inner_semimajor_axis= 1.0 |units.AU, outer_semimajor_axis= 100.0 |units.AU,
            inner_eccentricity= 0.1, outer_eccentricity= 0.5,
            relative_inclination= 80.0*np.pi/180.0,
            inner_argument_of_pericenter= 0.1, outer_argument_of_pericenter= 0.5,
            inner_longitude_of_ascending_node= 0.0,
            metallicity= 0.02,
            tend= 5.0 |units.Myr, number = 0, maximum_radius_change_factor = 0.005,
            stop_at_mass_transfer = True, stop_at_init_mass_transfer = True, stop_at_outer_mass_transfer = True,
            stop_at_stable_mass_transfer = True, stop_at_eccentric_stable_mass_transfer = True,
            stop_at_unstable_mass_transfer = False, stop_at_eccentric_unstable_mass_transfer = False,
            stop_at_merger = True, stop_at_disintegrated = True, stop_at_inner_collision = True, stop_at_outer_collision = True, 
            stop_at_dynamical_instability = True, stop_at_semisecular_regime = False, stop_at_SN = False,  SN_kick_distr = 2, stop_at_CPU_time = False,
            max_CPU_time = 3600.0, file_name = "triple.hdf", file_type = "hdf5", dir_plots = ""):


    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    inner_eccentricity = float(inner_eccentricity)
    outer_eccentricity = float(outer_eccentricity)
    relative_inclination = float(relative_inclination)
    inner_argument_of_pericenter = float(inner_argument_of_pericenter)
    outer_argument_of_pericenter = float(outer_argument_of_pericenter)
    inner_longitude_of_ascending_node = float(inner_longitude_of_ascending_node)

    triple_class_object = Triple_Class(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, 
            metallicity, tend, number, maximum_radius_change_factor,  
            stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer,
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime, stop_at_SN, SN_kick_distr, stop_at_CPU_time,
            max_CPU_time, file_name, file_type, dir_plots)


    if triple_class_object.triple.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given triple are incorrect.')
        return triple_class_object # no codes initialized yet
    elif stop_at_semisecular_regime == True and triple_class_object.triple.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is in the semisecular regime at initialization.')
    elif triple_class_object.triple.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is dynamically unstable at initialization.')
    elif triple_class_object.triple.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given triple at initialization.')
    else:    
        triple_class_object.evolve_model()
        if REPORT_DEBUG:
            plot_function(triple_class_object, dir_plots)
            triple_class_object.print_stellar_system()
    triple_class_object.stellar_code.stop()
    triple_class_object.secular_code.stop()


    return triple_class_object
#-----

#-----
#for running triple.py from the commandline
def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-M", unit=units.MSun, 
                      dest="inner_primary_mass", type="float", default = 1.3|units.MSun,
                      help="inner primary mass [%default]")
    parser.add_option("-m",  unit=units.MSun, 
                      dest="inner_secondary_mass", type="float", default = 0.5|units.MSun,
                      help="inner secondary mass [%default]")
    parser.add_option("-l",  unit=units.MSun, 
                      dest="outer_mass", type="float", default = 0.5|units.MSun,
                      help="outer mass [%default]")

    parser.add_option("-A", unit=units.RSun,
                      dest="inner_semimajor_axis", type="float", 
                      default = 200.0 |units.RSun,
                      help="inner semi major axis [%default]")
    parser.add_option("-a", unit=units.RSun,
                      dest="outer_semimajor_axis", type="float", 
                      default = 20000.0 |units.RSun,
                      help="outer semi major axis [%default]")
    parser.add_option("-E",
                      dest="inner_eccentricity", type="float", default = 0.1,
                      help="inner eccentricity [%default]")
    parser.add_option("-e",
                      dest="outer_eccentricity", type="float", default = 0.5,
                      help="outer eccentricity [%default]")
    parser.add_option("-i","-I",
                      dest="relative_inclination", type="float", default = 80.0*np.pi/180.0,
                      help="relative inclination [rad] [%default]")
    parser.add_option("-G",
                      dest="inner_argument_of_pericenter", type="float", default = 0.1,
                      help="inner argument of pericenter [rad] [%default]")
    parser.add_option("-g",
                      dest="outer_argument_of_pericenter", type="float", default = 0.5,
                      help="outer argument of pericenter [rad] [%default]")
    parser.add_option("-O",
                      dest="inner_longitude_of_ascending_node", type="float", default = 0.0,
                      help="inner longitude of ascending node [rad] [%default]")
##             outer longitude of ascending nodes = inner - pi               
#    parser.add_option("-o",
#                      dest="outer_longitude_of_ascending_node", type="float", default = 0.0,
#                      help="outer longitude of ascending node [rad] [%default]")


    #0  No kick 
    #1  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
    #2  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  scaled down for bh
    #3  Arzoumanian ea 2002, 568, 289
    #4  Arzoumanian ea 2002, 568, 289 scaled down for bh
    #5  Hansen & Phinney 1997, 291, 569
    #6  Hansen & Phinney 1997, 291, 569 scaled down for bh
    #7  Paczynski 1990, 348, 485
    #8  Paczynski 1990, 348, 485 scaled down for bh
    #9  Verbunt, Igoshev & Cator, 2017, 608, 57
    #10  Verbunt, Igoshev & Cator, 2017, 608, 57 scaled down for bh #default
    parser.add_option("--SN_kick_distr", dest="SN_kick_distr",  type="int", default = 10,
                      help="which supernova kick distribution [%default]")                      


    parser.add_option("-z", dest="metallicity", type="float", default = 0.02,
                      help="metallicity [%default] %unit")
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 5.0 |units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("-N", "-n", dest="number", type="int", default = 0,
                      help="number of system [%default]")
    parser.add_option("-r", dest="maximum_radius_change_factor", type="float", default = 0.01,
                      help="maximum_radius_change_factor [%default] %unit")

#    parser.add_option("--tidal", dest="tidal_terms", action="store_false", default = True, 
#                      help="tidal terms included [%default] %unit")

    parser.add_option("--no_stop_at_mass_transfer", dest="stop_at_mass_transfer", action="store_false", default = True,
                      help="stop at mass transfer [%default] %unit")
    parser.add_option("--no_stop_at_init_mass_transfer", dest="stop_at_init_mass_transfer", action="store_false", default = True,
                      help="stop if initially mass transfer[%default] %unit")
    parser.add_option("--no_stop_at_outer_mass_transfer", dest="stop_at_outer_mass_transfer", action="store_false", default = True,
                      help="stop at triple mass transfer [%default] %unit")
                      
#   if stop_at_mass_transfer is False, the following 4 stopping conditions can be used to further specify.
#   if stop_at_mass_transfer is True, the following 4 are ignored.
    parser.add_option("--stop_at_stable_mass_transfer", dest="stop_at_stable_mass_transfer", action="store_true", default = False,
                      help="stop at stable mass transfer [%default] %unit")
    parser.add_option("--stop_at_eccentric_stable_mass_transfer", dest="stop_at_eccentric_stable_mass_transfer", action="store_true",                                               
                    default = False, help="stop at eccentric stable mass transfer [%default] %unit")
    #unstable mass transfer leads to common-envelope evolution
    parser.add_option("--stop_at_unstable_mass_transfer", dest="stop_at_unstable_mass_transfer", action="store_true", 
                    default = False, help="stop at unstable mass transfer [%default] %unit")
    parser.add_option("--stop_at_eccentric_unstable_mass_transfer", dest="stop_at_eccentric_unstable_mass_transfer", 
                    action="store_true", default = False, help="stop at eccentric unstable mass transfer [%default] %unit")

    parser.add_option("--no_stop_at_merger", dest="stop_at_merger", action="store_false", default = True, 
                      help="stop at merger [%default] %unit")
    parser.add_option("--no_stop_at_disintegrated", dest="stop_at_disintegrated", action="store_false", default = True,
                      help="stop at disintegrated [%default] %unit")
    parser.add_option("--no_stop_at_inner_collision", dest="stop_at_inner_collision", action="store_false",default = True,
                      help="stop at collision in inner binary[%default] %unit")
    parser.add_option("--no_stop_at_outer_collision", dest="stop_at_outer_collision", action="store_false",default = True,
                      help="stop at collision in outer binary[%default] %unit")
    parser.add_option("--no_stop_at_dynamical_instability", dest="stop_at_dynamical_instability", action="store_false", default = True,
                      help="stop at dynamical instability [%default] %unit")
    parser.add_option("--stop_at_semisecular_regime", dest="stop_at_semisecular_regime", action="store_true", default = False,
                      help="stop at semisecular regime [%default] %unit")

    parser.add_option("--stop_at_SN", dest="stop_at_SN", action="store_true", default = False,
                      help="stop at supernova [%default] %unit")
    parser.add_option("--stop_at_CPU_time", dest="stop_at_CPU_time", action="store_true", default = False,
                      help="stop at CPU time [%default] %unit")
    parser.add_option("--max_CPU_time", dest="max_CPU_time", type="float", default = 3600.0,
                      help="max CPU time [%default] %unit")
                      
                      
    parser.add_option("-f", dest="file_name", type ="string", default = "TRES.hdf",#"TRES.txt"
                      help="file name[%default]")
    parser.add_option("-F", dest="file_type", type ="string", default = "hdf5",#"txt"
                      help="file type[%default]")
    parser.add_option("--dir_plots", dest="dir_plots", type ="string", default = "",#"txt"
                      help="directory for plots for debugging mode [%default]")



                                           
    options, args = parser.parse_args()
    return options.__dict__


if __name__ == '__main__':
    options = parse_arguments()

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")


    triple_class_object = Triple_Class(**options)  

    if triple_class_object.triple.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given triple are incorrect.' )   
        # no codes initialized yet
        exit(0)
    elif options['stop_at_semisecular_regime'] == True and triple_class_object.triple.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is in the semisecular regime at initialization.')
    elif triple_class_object.triple.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is dynamically unstable at initialization.')
    elif triple_class_object.triple.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given triple at initialization.')
    else:    
        triple_class_object.evolve_model()
        if REPORT_DEBUG:
            plot_function(triple_class_object, options['dir_plots'])
            triple_class_object.print_stellar_system()



        if REPORT_TRIPLE_EVOLUTION:
            print('Simulation has finished succesfully')
            
    print('\nYou have used the TRES triple evolution code. Literature reference:')
    print('** Toonen, Hamers & Portegies Zwart 2016, ComAC, 3, 6T:')
    print('... "The evolution of hierarchical triple star-systems" ')
            
    triple_class_object.stellar_code.stop()
    triple_class_object.secular_code.stop()
