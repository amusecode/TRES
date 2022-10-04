"""
Interface for SecularTriple
Designed to work together with TPS

Adrian Hamers 2015
"""

from amuse.community import *
from amuse.units import units,constants
import sys

### units used internally in the ODE solver ###
unit_l = units.AU
unit_m = units.MSun
unit_t = 1.0e6*units.yr
unit_h = unit_m*unit_l**2/unit_t ### angular momentum

### CVODE flags ###
CV_SUCCESS=0
CV_ROOT_RETURN=2
CV_WARNING=99
CV_TOO_MUCH_WORK=-1
CV_TOO_MUCH_ACC=-2
CV_ERR_FAILURE=-3
CV_CONV_FAILURE=-4
CV_LINIT_FAIL=-5
CV_LSETUP_FAIL=-6
CV_LSOLVE_FAIL=-7
CV_RHSFUNC_FAIL=-8
CV_FIRST_RHSFUNC_ERR=-9
CV_REPTD_RHSFUNC_ERR=-10
CV_UNREC_RHSFUNC_ERR=-11
CV_RTFUNC_FAIL=-12
CV_MEM_FAIL=-20
CV_MEM_NULL=-21
CV_ILL_INPUT=-22
CV_NO_MALLOC=-23
CV_BAD_K=-24
CV_BAD_T=-25
CV_BAD_DKY=-26
CV_TOO_CLOSE=-27


class SecularTripleInterface(CodeInterface):
    include_headers = ['src/main_code.h','src/ODE_system.h']

    def __init__(self, **options):
         CodeInterface.__init__(self, **options)
        
    @legacy_function
    def evolve():
        function = LegacyFunctionSpecification()
        function.addParameter('stellar_type1', dtype='int32', direction=function.IN)
        function.addParameter('stellar_type2', dtype='int32', direction=function.IN)
        function.addParameter('stellar_type3', dtype='int32', direction=function.IN)
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('m1_convective_envelope', dtype='float64', direction=function.IN)
        function.addParameter('m2_convective_envelope', dtype='float64', direction=function.IN)
        function.addParameter('m3_convective_envelope', dtype='float64', direction=function.IN)
        function.addParameter('R1', dtype='float64', direction=function.IN)
        function.addParameter('R2', dtype='float64', direction=function.IN)
        function.addParameter('R3', dtype='float64', direction=function.IN)   
        function.addParameter('R1_convective_envelope', dtype='float64', direction=function.IN)
        function.addParameter('R2_convective_envelope', dtype='float64', direction=function.IN)
        function.addParameter('R3_convective_envelope', dtype='float64', direction=function.IN)   
        function.addParameter('luminosity_star1', dtype='float64', direction=function.IN)
        function.addParameter('luminosity_star2', dtype='float64', direction=function.IN)
        function.addParameter('luminosity_star3', dtype='float64', direction=function.IN)   
        function.addParameter('spin_angular_frequency1', dtype='float64', direction=function.IN)   
        function.addParameter('spin_angular_frequency2', dtype='float64', direction=function.IN)   
        function.addParameter('spin_angular_frequency3', dtype='float64', direction=function.IN)                   
        function.addParameter('AMC_star1', dtype='float64', direction=function.IN)   
        function.addParameter('AMC_star2', dtype='float64', direction=function.IN)   
        function.addParameter('AMC_star3', dtype='float64', direction=function.IN)           
        function.addParameter('gyration_radius_star1', dtype='float64', direction=function.IN)   
        function.addParameter('gyration_radius_star2', dtype='float64', direction=function.IN)   
        function.addParameter('gyration_radius_star3', dtype='float64', direction=function.IN)   
        function.addParameter('moment_of_inertia_star1', dtype='float64', direction=function.IN)
        function.addParameter('moment_of_inertia_star2', dtype='float64', direction=function.IN)
        function.addParameter('moment_of_inertia_star3', dtype='float64', direction=function.IN)
        function.addParameter('moment_of_inertia_dot_star1', dtype='float64', direction=function.IN)
        function.addParameter('moment_of_inertia_dot_star2', dtype='float64', direction=function.IN)
        function.addParameter('moment_of_inertia_dot_star3', dtype='float64', direction=function.IN)
#        function.addParameter('k_div_T_tides_star1', dtype='float64', direction=function.IN)   
#        function.addParameter('k_div_T_tides_star2', dtype='float64', direction=function.IN)   
#        function.addParameter('k_div_T_tides_star3', dtype='float64', direction=function.IN)   
        function.addParameter('a_in', dtype='float64', direction=function.IN)
        function.addParameter('a_out', dtype='float64', direction=function.IN)
        function.addParameter('e_in', dtype='float64', direction=function.IN)
        function.addParameter('e_out', dtype='float64', direction=function.IN)
        function.addParameter('INCL_in', dtype='float64', direction=function.IN)
        function.addParameter('INCL_out', dtype='float64', direction=function.IN)
        function.addParameter('AP_in', dtype='float64', direction=function.IN)
        function.addParameter('AP_out', dtype='float64', direction=function.IN)
        function.addParameter('LAN_in', dtype='float64', direction=function.IN)
        function.addParameter('LAN_out', dtype='float64', direction=function.IN)
        function.addParameter('star1_is_donor', dtype='bool', direction=function.IN)
        function.addParameter('star2_is_donor', dtype='bool', direction=function.IN)
        function.addParameter('star3_is_donor', dtype='bool', direction=function.IN)
        function.addParameter('wind_mass_loss_rate_star1', dtype='float64', direction=function.IN)
        function.addParameter('wind_mass_loss_rate_star2', dtype='float64', direction=function.IN)
        function.addParameter('wind_mass_loss_rate_star3', dtype='float64', direction=function.IN)
        function.addParameter('time_derivative_of_radius_star1', dtype='float64', direction=function.IN)
        function.addParameter('time_derivative_of_radius_star2', dtype='float64', direction=function.IN)
        function.addParameter('time_derivative_of_radius_star3', dtype='float64', direction=function.IN)        
        function.addParameter('inner_mass_transfer_rate', dtype='float64', direction=function.IN)
        function.addParameter('outer_mass_transfer_rate', dtype='float64', direction=function.IN)
        function.addParameter('inner_accretion_efficiency_wind_child1_to_child2', dtype='float64', direction=function.IN)
        function.addParameter('inner_accretion_efficiency_wind_child2_to_child1', dtype='float64', direction=function.IN)
        function.addParameter('outer_accretion_efficiency_wind_child1_to_child2', dtype='float64', direction=function.IN)
        function.addParameter('outer_accretion_efficiency_wind_child2_to_child1', dtype='float64', direction=function.IN)
        function.addParameter('inner_accretion_efficiency_mass_transfer', dtype='float64', direction=function.IN)        
        function.addParameter('outer_accretion_efficiency_mass_transfer', dtype='float64', direction=function.IN)                      
        function.addParameter('inner_specific_AM_loss_mass_transfer', dtype='float64', direction=function.IN)                        
        function.addParameter('outer_specific_AM_loss_mass_transfer', dtype='float64', direction=function.IN)
        function.addParameter('inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2', dtype='float64', direction=function.IN)                        
        function.addParameter('inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1', dtype='float64', direction=function.IN)                        
        function.addParameter('t', dtype='float64', direction=function.IN)        
        function.addParameter('global_time_step', dtype='float64', direction=function.IN)        
        function.addParameter('m1_output', dtype='float64', direction=function.OUT) 
        function.addParameter('m2_output', dtype='float64', direction=function.OUT) 
        function.addParameter('m3_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R1_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R2_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R3_output', dtype='float64', direction=function.OUT)
        function.addParameter('spin_angular_frequency1_output', dtype='float64', direction=function.OUT)
        function.addParameter('spin_angular_frequency2_output', dtype='float64', direction=function.OUT)
        function.addParameter('spin_angular_frequency3_output', dtype='float64', direction=function.OUT)
        function.addParameter('a_in_output', dtype='float64', direction=function.OUT) 
        function.addParameter('a_out_output', dtype='float64', direction=function.OUT) 
        function.addParameter('e_in_output', dtype='float64', direction=function.OUT) 
        function.addParameter('e_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('delta_e_in_output', dtype='float64', direction=function.OUT)
        function.addParameter('INCL_in_output', dtype='float64', direction=function.OUT)
        function.addParameter('INCL_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('INCL_in_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('AP_in_output', dtype='float64', direction=function.OUT)
        function.addParameter('AP_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('LAN_in_output', dtype='float64', direction=function.OUT)
        function.addParameter('LAN_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('t_output', dtype='float64', direction=function.OUT)
        function.addParameter('CVODE_flag', dtype='int32', direction=function.OUT)
        function.addParameter('root_finding_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def compute_orbital_vectors_from_orbital_elements():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('a1', dtype='float64', direction=function.IN)
        function.addParameter('a2', dtype='float64', direction=function.IN)
        function.addParameter('e1', dtype='float64', direction=function.IN)
        function.addParameter('e2', dtype='float64', direction=function.IN)
        function.addParameter('INCL_rel', dtype='float64', direction=function.IN)        
        function.addParameter('AP1', dtype='float64', direction=function.IN)        
        function.addParameter('AP2', dtype='float64', direction=function.IN)        
        function.addParameter('LAN1', dtype='float64', direction=function.IN)        
        function.addParameter('LAN2', dtype='float64', direction=function.IN)        
        function.addParameter('e1_vec_x', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e1_vec_y', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e1_vec_z', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e2_vec_x', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e2_vec_y', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e2_vec_z', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h1_vec_x', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h1_vec_y', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h1_vec_z', dtype='float64',direction=function.OUT,description = "")        
        function.addParameter('h2_vec_x', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h2_vec_y', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h2_vec_z', dtype='float64',direction=function.OUT,description = "")        
        function.result_type = 'float64'
        return function
    
    @legacy_function
    def compute_orbital_elements_from_orbital_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('e1_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e1_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e1_vec_z', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e2_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e2_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e2_vec_z', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h1_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h1_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h1_vec_z', dtype='float64',direction=function.IN,description = "")        
        function.addParameter('h2_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h2_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h2_vec_z', dtype='float64',direction=function.IN,description = "")  
        function.addParameter('a1', dtype='float64', direction=function.OUT)
        function.addParameter('a2', dtype='float64', direction=function.OUT)
        function.addParameter('e1', dtype='float64', direction=function.OUT)
        function.addParameter('e2', dtype='float64', direction=function.OUT)
        function.addParameter('INCL_rel', dtype='float64', direction=function.OUT)        
        function.addParameter('AP1', dtype='float64', direction=function.OUT)        
        function.addParameter('AP2', dtype='float64', direction=function.OUT)        
        function.addParameter('LAN1', dtype='float64', direction=function.OUT)        
        function.addParameter('LAN2', dtype='float64', direction=function.OUT)        
        function.result_type = 'float64'
        return function

    @legacy_function
    def compute_effect_of_SN_on_orbital_vectors():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('e1_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e1_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e1_vec_z', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e2_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e2_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('e2_vec_z', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h1_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h1_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h1_vec_z', dtype='float64',direction=function.IN,description = "")        
        function.addParameter('h2_vec_x', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h2_vec_y', dtype='float64',direction=function.IN,description = "")
        function.addParameter('h2_vec_z', dtype='float64',direction=function.IN,description = "") 
        function.addParameter('Vkick1_vec_x', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick1_vec_y', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick1_vec_z', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick2_vec_x', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick2_vec_y', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick2_vec_z', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick3_vec_x', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick3_vec_y', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('Vkick3_vec_z', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('delta_m1', dtype='float64', direction=function.IN)
        function.addParameter('delta_m2', dtype='float64', direction=function.IN)
        function.addParameter('delta_m3', dtype='float64', direction=function.IN)
        function.addParameter('f1', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('f2', dtype='float64',direction=function.IN,description = "")         
        function.addParameter('V1_prime_x', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V1_prime_y', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V1_prime_z', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V2_prime_x', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V2_prime_y', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V2_prime_z', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V3_prime_x', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V3_prime_y', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('V3_prime_z', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e1_vec_x_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e1_vec_y_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e1_vec_z_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e2_vec_x_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e2_vec_y_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('e2_vec_z_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h1_vec_x_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h1_vec_y_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h1_vec_z_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h2_vec_x_p', dtype='float64',direction=function.OUT,description = "")
        function.addParameter('h2_vec_y_p', dtype='float64',direction=function.OUT,description = "") 
        function.addParameter('h2_vec_z_p', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('cos_phi1', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('cos_phi2', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('R', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('v_sys', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('r1dotv1', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('r1dotvk', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('v1dotvk', dtype='float64',direction=function.OUT,description = "")   
        function.addParameter('r1dotvsys', dtype='float64',direction=function.OUT,description = "")           
        function.addParameter('v1dotvsys', dtype='float64',direction=function.OUT,description = "")           
        function.addParameter('r2dot2', dtype='float64',direction=function.OUT,description = "")           
        function.addParameter('r2dotvsys', dtype='float64',direction=function.OUT,description = "")           
        function.addParameter('v2dotvsys', dtype='float64',direction=function.OUT,description = "")           
        function.addParameter('r1dotr2', dtype='float64',direction=function.OUT,description = "")           
        function.addParameter('r1dotv2', dtype='float64',direction=function.OUT,description = "")           
        function.result_type = 'float64'
        return function
    
    @legacy_function
    def a_out_div_a_in_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('a_in', dtype='float64', direction=function.IN)
        function.addParameter('a_out', dtype='float64', direction=function.IN)
        function.addParameter('e_in', dtype='float64', direction=function.IN)
        function.addParameter('e_out', dtype='float64', direction=function.IN)
        function.addParameter('itot', dtype='float64', direction=function.IN)        
        function.addParameter('stability_limit_specification', dtype='int32', direction=function.IN)
        function.result_type = 'float64'
        return function

    @legacy_function
    def a_out_div_a_in_semisecular_regime():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('e_in', dtype='float64', direction=function.IN)
        function.addParameter('e_out', dtype='float64', direction=function.IN)
        function.addParameter('check_for_semisecular_regime_parameter', dtype='float64', direction=function.IN)
        function.result_type = 'float64'
        return function

    @legacy_function
    def roche_radius():
        function = LegacyFunctionSpecification()
        function.addParameter('rp', dtype='float64', direction=function.IN)
        function.addParameter('q', dtype='float64', direction=function.IN)
        function.addParameter('e', dtype='float64', direction=function.IN)
        function.addParameter('f', dtype='float64', direction=function.IN)
        function.addParameter('roche_radius_specification', dtype='int32', direction=function.IN)
        function.result_type = 'float64'
        return function


    @legacy_function
    def get_relative_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "Relative tolerance, default 1e-10")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_relative_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "Relative tolerance, default 1e-10")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "",unit=1.0/unit_t)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "",unit=1.0/unit_t)
        function.result_type = 'int32'
        return function

#    @legacy_function
#    def get_dynamical_stability_model():
#        function = LegacyFunctionSpecification()
#        function.addParameter('value', dtype='int32',direction=function.OUT,description = "")
#        function.result_type = 'int32'
#        return function
#
#    @legacy_function
#    def set_dynamical_stability_model():
#        function = LegacyFunctionSpecification()
#        function.addParameter('value', dtype='int32',direction=function.IN,description = "")
#        function.result_type = 'int32'
#        return function
        
    @legacy_function
    def get_input_precision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_input_precision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_linear_solver():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.OUT,description = "linear solver")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_linear_solver():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.IN,description = "linear solver")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_equations_of_motion_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_equations_of_motion_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_roche_radius_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_roche_radius_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_stability_limit_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_stability_limit_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_check_for_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_check_for_inner_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_inner_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_outer_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_outer_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_inner_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_inner_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_outer_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_outer_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_include_quadrupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include quadrupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_quadrupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include quadrupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function         

    @legacy_function
    def set_include_octupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include octupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_octupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include octupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function   

    @legacy_function
    def get_include_1PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 1PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 1PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_1PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 1PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 1PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_1PN_inner_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 1PN 'interaction terms' in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_inner_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 1PN 'interaction terms' in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_25PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 2.5PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_25PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 2.5PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_25PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 2.5PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_25PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 2.5PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_inner_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the inner binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_inner_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the inner binary system")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_outer_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_outer_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ignore_tertiary():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_ignore_tertiary():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_include_inner_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_inner_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_outer_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_outer_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_magnetic_braking_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_magnetic_braking_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_spin_radius_mass_coupling_terms_star1():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_spin_radius_mass_coupling_terms_star1():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_spin_radius_mass_coupling_terms_star2():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_spin_radius_mass_coupling_terms_star2():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_spin_radius_mass_coupling_terms_star3():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_spin_radius_mass_coupling_terms_star3():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_include_inner_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_inner_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_outer_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_outer_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_linear_mass_change():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_linear_mass_change():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_linear_radius_change():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_linear_radius_change():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_check_for_dynamical_stability_at_initialisation():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    

    @legacy_function
    def get_check_for_dynamical_stability_at_initialisation():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_semisecular_regime_at_initialisation():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    

    @legacy_function
    def get_check_for_semisecular_regime_at_initialisation():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_verbose():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    

    @legacy_function
    def get_verbose():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function




    ### addition Jan 2017 ###
    @legacy_function
    def get_check_for_semisecular_regime():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_semisecular_regime():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_semisecular_regime_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_semisecular_regime_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "")
        function.result_type = 'int32'
        return function




class SecularTriple(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  SecularTripleInterface(**options), **options)
        self.model_time = 0.0 | units.Myr
        self.evolve_further_after_root_was_found = False ### in some cases, it is desirable to integrate until the given end time, despite the finding of a root
        self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe = True ### this may seem like a strange `feature', but exists for consistency with the stellar/binary evolution code

    def define_parameters(self, object):
        
        object.add_method_parameter(
            "get_model_time",
            "set_model_time",
            "model_time",
            "model_time",
            default_value = 0.0 | units.Myr
        )
        object.add_method_parameter(
            "get_relative_tolerance",
            "set_relative_tolerance",
            "relative_tolerance",
            "Relative tolerance, default 1e-10",
            default_value = 1.0e-10
        )
        object.add_method_parameter(
            "get_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero",
            "set_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero",
            "threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero",
            "",
            default_value = 1.0e-12
        )
        object.add_method_parameter(
            "get_threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero",
            "set_threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero",
            "threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero",
            "",
            default_value = 1.0e-8
        )
#        object.add_method_parameter(
#            "get_dynamical_stability_model",
#            "set_dynamical_stability_model",
#            "dynamical_stability_model",
#            "",
#            default_value = 0
#        )
        object.add_method_parameter(
            "get_input_precision",
            "set_input_precision",
            "input_precision",
            "Input_precision, default 1e-5",
            default_value = 1.0e-5
        )        
        object.add_method_parameter(
            "get_linear_solver",
            "set_linear_solver",
            "linear_solver",
            "linear_solver, default 0",
            default_value = 0
        )             
        object.add_method_parameter(
            "get_equations_of_motion_specification",
            "set_equations_of_motion_specification",
            "equations_of_motion_specification",
            "equations_of_motion_specification",
            default_value = 0
        )
        object.add_method_parameter(
            "get_roche_radius_specification",
            "set_roche_radius_specification",
            "roche_radius_specification",
            "roche_radius_specification",
            default_value = 0
        )
        object.add_method_parameter(
            "get_stability_limit_specification",
            "set_stability_limit_specification",
            "stability_limit_specification",
            "stability_limit_specification",
            default_value = 0
        )
        object.add_method_parameter(
            "get_check_for_dynamical_stability",
            "set_check_for_dynamical_stability",
            "check_for_dynamical_stability",
            "check_for_dynamical_stability",
            default_value = True
        )       
        object.add_method_parameter(
            "get_check_for_inner_collision",
            "set_check_for_inner_collision",
            "check_for_inner_collision",
            "check_for_inner_collision",
            default_value = True
        )          
        object.add_method_parameter(
            "get_check_for_outer_collision",
            "set_check_for_outer_collision",
            "check_for_outer_collision",
            "check_for_outer_collision",
            default_value = False
        )  
        object.add_method_parameter(
            "get_check_for_inner_RLOF",
            "set_check_for_inner_RLOF",
            "check_for_inner_RLOF",
            "check_for_inner_RLOF",
            default_value = False
        )
        object.add_method_parameter(
            "get_check_for_outer_RLOF",
            "set_check_for_outer_RLOF",
            "check_for_outer_RLOF",
            "check_for_outer_RLOF",
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_quadrupole_terms",
            "set_include_quadrupole_terms",
            "include_quadrupole_terms",
            "Whether or not to include quadrupole terms in the secular equations of motion", 
            default_value = True
        )        
        object.add_method_parameter(
            "get_include_octupole_terms",
            "set_include_octupole_terms",
            "include_octupole_terms",
            "Whether or not to include octupole terms in the secular equations of motion", 
            default_value = True
        )
        object.add_method_parameter(
            "get_include_1PN_inner_terms",
            "set_include_1PN_inner_terms",
            "include_1PN_inner_terms",
            "Include 1PN inner binary terms in the equations of motion", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_1PN_outer_terms",
            "set_include_1PN_outer_terms",
            "include_1PN_outer_terms",
            "Include 1PN outer binary terms in the equations of motion", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_1PN_inner_outer_terms",
            "set_include_1PN_inner_outer_terms",
            "include_1PN_inner_outer_terms",
            "Include 1PN interaction terms in the equations of motion", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_25PN_inner_terms",
            "set_include_25PN_inner_terms",
            "include_25PN_inner_terms",
            "Include 25PN inner binary terms in the equations of motion", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_25PN_outer_terms",
            "set_include_25PN_outer_terms",
            "include_25PN_outer_terms",
            "Include 25PN outer binary terms in the equations of motion", 
            default_value = False
        )
        object.add_method_parameter(
            "get_ignore_tertiary",
            "set_ignore_tertiary",
            "ignore_tertiary",
            "", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_inner_tidal_terms",
            "set_include_inner_tidal_terms",
            "include_inner_tidal_terms",
            "Whether or not to include the effects of tides in the inner binary system", 
            default_value = False
        )                           
        object.add_method_parameter(
            "get_include_outer_tidal_terms",
            "set_include_outer_tidal_terms",
            "include_outer_tidal_terms",
            "Whether or not to include the effects of tides in the outer binary system", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_inner_wind_terms",
            "set_include_inner_wind_terms",
            "include_inner_wind_terms",
            "..", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_outer_wind_terms",
            "set_include_outer_wind_terms",
            "include_outer_wind_terms",
            "..", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_magnetic_braking_terms",
            "set_include_magnetic_braking_terms",
            "include_magnetic_braking_terms",
            "..", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_spin_radius_mass_coupling_terms_star1",
            "set_include_spin_radius_mass_coupling_terms_star1",
            "include_spin_radius_mass_coupling_terms_star1",
            "..", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_spin_radius_mass_coupling_terms_star2",
            "set_include_spin_radius_mass_coupling_terms_star2",
            "include_spin_radius_mass_coupling_terms_star2",
            "..", 
            default_value = False
        )  
        object.add_method_parameter(
            "get_include_spin_radius_mass_coupling_terms_star3",
            "set_include_spin_radius_mass_coupling_terms_star3",
            "include_spin_radius_mass_coupling_terms_star3",
            "..", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_inner_RLOF_terms",
            "set_include_inner_RLOF_terms",
            "include_inner_RLOF_terms",
            "..", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_outer_RLOF_terms",
            "set_include_outer_RLOF_terms",
            "include_outer_RLOF_terms",
            "..", 
            default_value = False
        )    
        object.add_method_parameter(
            "get_include_linear_mass_change",
            "set_include_linear_mass_change",
            "include_linear_mass_change",
            "..", 
            default_value = False
        )            
        object.add_method_parameter(
            "get_include_linear_radius_change",
            "set_include_linear_radius_change",
            "include_linear_radius_change",
            "..", 
            default_value = False
        )   
        object.add_method_parameter(
            "get_check_for_dynamical_stability_at_initialisation",
            "set_check_for_dynamical_stability_at_initialisation",
            "check_for_dynamical_stability_at_initialisation",
            "..", 
            default_value = True
        )     
        object.add_method_parameter(
            "get_check_for_semisecular_regime_at_initialisation",
            "set_check_for_semisecular_regime_at_initialisation",
            "check_for_semisecular_regime_at_initialisation",
            "..", 
            default_value = True
        )     
        object.add_method_parameter(
            "get_verbose",
            "set_verbose",
            "verbose",
            "..", 
            default_value = True
        )    
        
        
        
        ### addition Jan 2017 ###
        object.add_method_parameter(
            "get_check_for_semisecular_regime_parameter",
            "set_check_for_semisecular_regime_parameter",
            "check_for_semisecular_regime_parameter",
            "..", 
            default_value = 1
        )   
        object.add_method_parameter(
            "get_check_for_semisecular_regime",
            "set_check_for_semisecular_regime",
            "check_for_semisecular_regime",
            "..", 
            default_value = False
        )    
    def define_methods(self, object):
        unit_lum = unit_m*unit_l**2/(unit_t**3)
        object.add_method(
            "evolve",
            (
                units.stellar_type,         ### stellar_type1
                units.stellar_type,         ### stellar_type2
                units.stellar_type,         ### stellar_type3
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                unit_m,                     ### m1_convective_envelope
                unit_m,                     ### m2_convective_envelope
                unit_m,                     ### m3_convective_envelope
                unit_l,                     ### R1
                unit_l,                     ### R2
                unit_l,                     ### R3
                unit_l,                     ### R1_convective_envelope
                unit_l,                     ### R2_convective_envelope
                unit_l,                     ### R3_convective_envelope
                unit_lum,                   ### luminosity_star1
                unit_lum,                   ### luminosity_star2
                unit_lum,                   ### luminosity_star3
                1.0/unit_t,                 ### spin_angular_frequency1
                1.0/unit_t,                 ### spin_angular_frequency2
                1.0/unit_t,                 ### spin_angular_frequency3                                
                object.NO_UNIT,             ### AMC_star1
                object.NO_UNIT,             ### AMC_star2
                object.NO_UNIT,             ### AMC_star3
                object.NO_UNIT,             ### gyration_radius_star1
                object.NO_UNIT,             ### gyration_radius_star2
                object.NO_UNIT,             ### gyration_radius_star3
                unit_m*unit_l**2,           ### moment_of_inertia_star1
                unit_m*unit_l**2,           ### moment_of_inertia_star2
                unit_m*unit_l**2,           ### moment_of_inertia_star3
                unit_m*unit_l**2/unit_t,    ### moment_of_inertia_dot_star1
                unit_m*unit_l**2/unit_t,    ### moment_of_inertia_dot_star2
                unit_m*unit_l**2/unit_t,    ### moment_of_inertia_dot_star3
#                1.0/units.s,                ### k_div_T_tides_star1
#                1.0/units.s,                ### k_div_T_tides_star2
#                1.0/units.s,                ### k_div_T_tides_star3                                
                unit_l,                     ### a_in
                unit_l,                     ### a_out
                object.NO_UNIT,             ### e_in
                object.NO_UNIT,             ### e_out
                object.NO_UNIT,             ### INCL_in
                object.NO_UNIT,             ### INCL_out                
                object.NO_UNIT,             ### AP_in
                object.NO_UNIT,             ### AP_out
                object.NO_UNIT,             ### LAN_in
                object.NO_UNIT,             ### LAN_out
                object.NO_UNIT,             ### star1_is_donor
                object.NO_UNIT,             ### star2_is_donor
                object.NO_UNIT,             ### star3_is_donor
                unit_m/unit_t,              ### wind_mass_loss_rate_star1
                unit_m/unit_t,              ### wind_mass_loss_rate_star2
                unit_m/unit_t,              ### wind_mass_loss_rate_star3
                unit_l/unit_t,              ### time_derivative_of_radius_star1
                unit_l/unit_t,              ### time_derivative_of_radius_star2
                unit_l/unit_t,              ### time_derivative_of_radius_star3
                unit_m/unit_t,              ### inner_mass_transfer_rate                                
                unit_m/unit_t,              ### outer_mass_transfer_rate                                
                object.NO_UNIT,             ### inner_accretion_efficiency_wind_child1_to_child2
                object.NO_UNIT,             ### inner_accretion_efficiency_wind_child2_to_child1
                object.NO_UNIT,             ### outer_accretion_efficiency_wind_child1_to_child2
                object.NO_UNIT,             ### outer_accretion_efficiency_wind_child2_to_child1
                object.NO_UNIT,             ### inner_accretion_efficiency_mass_transfer
                object.NO_UNIT,             ### outer_accretion_efficiency_mass_transfer
                object.NO_UNIT,             ### inner_specific_AM_loss_mass_transfer
                object.NO_UNIT,             ### outer_specific_AM_loss_mass_transfer
                object.NO_UNIT,             ### inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2
                object.NO_UNIT,             ### inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1
                unit_t,                     ### t
                unit_t,                     ### dt
            ),
            (
                unit_m,                     ### m1_output
                unit_m,                     ### m2_output
                unit_m,                     ### m3_output
                unit_l,                     ### R1_output
                unit_l,                     ### R2_output
                unit_l,                     ### R3_output
                1.0/unit_t,                 ### spin_angular_frequency1_output
                1.0/unit_t,                 ### spin_angular_frequency2_output
                1.0/unit_t,                 ### spin_angular_frequency3_output                         
                unit_l,                     ### a_in_output
                unit_l,                     ### a_out_output
                object.NO_UNIT,             ### e_in_output
                object.NO_UNIT,             ### e_out_output
                object.NO_UNIT,             ### delta_e_in_output
                object.NO_UNIT,             ### INCL_in_output
                object.NO_UNIT,             ### INCL_out_output
                object.NO_UNIT,             ### INCL_in_out_output
                object.NO_UNIT,             ### AP_in_output
                object.NO_UNIT,             ### AP_out_output                                
                object.NO_UNIT,             ### LAN_in_output
                object.NO_UNIT,             ### LAN_out_output                                
                unit_t,                     ### t_output
                object.NO_UNIT,             ### CVODE_flag
                object.NO_UNIT,             ### root_finding_flag
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "compute_orbital_vectors_from_orbital_elements",
            (
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                unit_l,                     ### a1
                unit_l,                     ### a2
                object.NO_UNIT,             ### e1
                object.NO_UNIT,             ### e2
                object.NO_UNIT,             ### INCL_rel
                object.NO_UNIT,             ### AP1
                object.NO_UNIT,             ### AP2
                object.NO_UNIT,             ### LAN1
                object.NO_UNIT,             ### LAN2
            ),
            (
                object.NO_UNIT,             ### e1_vec_x
                object.NO_UNIT,             ### e1_vec_y
                object.NO_UNIT,             ### e1_vec_z
                object.NO_UNIT,             ### e2_vec_x
                object.NO_UNIT,             ### e2_vec_y
                object.NO_UNIT,             ### e2_vec_z
                unit_h,                     ### h1_vec_x
                unit_h,                     ### h1_vec_y
                unit_h,                     ### h1_vec_z
                unit_h,                     ### h2_vec_x
                unit_h,                     ### h2_vec_y
                unit_h,                     ### h2_vec_z
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "compute_orbital_elements_from_orbital_vectors",
            (
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                object.NO_UNIT,             ### e1_vec_x
                object.NO_UNIT,             ### e1_vec_y
                object.NO_UNIT,             ### e1_vec_z
                object.NO_UNIT,             ### e2_vec_x
                object.NO_UNIT,             ### e2_vec_y
                object.NO_UNIT,             ### e2_vec_z
                unit_h,                     ### h1_vec_x
                unit_h,                     ### h1_vec_y
                unit_h,                     ### h1_vec_z
                unit_h,                     ### h2_vec_x
                unit_h,                     ### h2_vec_y
                unit_h,                     ### h2_vec_z
            ),
            (
                unit_l,                     ### a1
                unit_l,                     ### a2
                object.NO_UNIT,             ### e1
                object.NO_UNIT,             ### e2
                object.NO_UNIT,             ### INCL_rel
                object.NO_UNIT,             ### AP1
                object.NO_UNIT,             ### AP2
                object.NO_UNIT,             ### LAN1
                object.NO_UNIT,             ### LAN2
                object.ERROR_CODE,
            )
        )

        object.add_method(
            "compute_effect_of_SN_on_orbital_vectors",
            (
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                object.NO_UNIT,             ### e1_vec_x
                object.NO_UNIT,             ### e1_vec_y
                object.NO_UNIT,             ### e1_vec_z
                object.NO_UNIT,             ### e2_vec_x
                object.NO_UNIT,             ### e2_vec_y
                object.NO_UNIT,             ### e2_vec_z
                unit_h,                     ### h1_vec_x
                unit_h,                     ### h1_vec_y
                unit_h,                     ### h1_vec_z
                unit_h,                     ### h2_vec_x
                unit_h,                     ### h2_vec_y
                unit_h,                     ### h2_vec_z
                unit_l/unit_t,              ### Vkick1_vec_x
                unit_l/unit_t,              ### Vkick1_vec_y
                unit_l/unit_t,              ### Vkick1_vec_z
                unit_l/unit_t,              ### Vkick2_vec_x
                unit_l/unit_t,              ### Vkick2_vec_y
                unit_l/unit_t,              ### Vkick2_vec_z
                unit_l/unit_t,              ### Vkick3_vec_x
                unit_l/unit_t,              ### Vkick3_vec_y
                unit_l/unit_t,              ### Vkick3_vec_z
                unit_m,                     ### delta_m1
                unit_m,                     ### delta_m2
                unit_m,                     ### delta_m3
                object.NO_UNIT,             ### f1
                object.NO_UNIT,             ### f2
            ),
            (
                unit_l/unit_t,              ### V1_prime_x
                unit_l/unit_t,              ### V1_prime_y
                unit_l/unit_t,              ### V1_prime_z
                unit_l/unit_t,              ### V2_prime_x
                unit_l/unit_t,              ### V2_prime_y
                unit_l/unit_t,              ### V2_prime_z
                unit_l/unit_t,              ### V3_prime_x
                unit_l/unit_t,              ### V3_prime_y
                unit_l/unit_t,              ### V3_prime_z
                object.NO_UNIT,             ### e1_vec_x_p
                object.NO_UNIT,             ### e1_vec_y_p
                object.NO_UNIT,             ### e1_vec_z_p
                object.NO_UNIT,             ### e2_vec_x_p
                object.NO_UNIT,             ### e2_vec_y_p
                object.NO_UNIT,             ### e2_vec_z_p
                unit_h,                     ### h1_vec_x_p
                unit_h,                     ### h1_vec_y_p
                unit_h,                     ### h1_vec_z_p
                unit_h,                     ### h2_vec_x_p
                unit_h,                     ### h2_vec_y_p
                unit_h,                     ### h2_vec_z_p
                object.NO_UNIT,             ### cos_phi1
                object.NO_UNIT,             ### cos_phi2
                unit_l,                     ### R
                unit_l/unit_t,              ### v_sys
                unit_l*unit_l/unit_t,       ### r1dotvk
                unit_l*unit_l/unit_t,       ### r1dotv1
                unit_l**2/(unit_t**2),      ### v1dotvk
                unit_l*unit_l/unit_t,       ### r1dotvsys
                unit_l**2/(unit_t**2),      ### v1dotvsys
                unit_l*unit_l/unit_t,       ### r2dotv2
                unit_l*unit_l/unit_t,       ### r2dotvsys
                unit_l**2/(unit_t**2),      ### v2dotvsys
                unit_l**2,                  ### r1dotr2
                unit_l**2/unit_t,           ### r1dotv2
                object.ERROR_CODE,
            )
        )
        
        object.add_method(
            "a_out_div_a_in_dynamical_stability",
            (
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                unit_l,                     ### a_in
                unit_l,                     ### a_out
                object.NO_UNIT,             ### e_in
                object.NO_UNIT,             ### e_out
                object.NO_UNIT,             ### itot -- NOTE: should be in radians
                object.NO_UNIT,             ### stability_limit_specification
                
            ),
            (
                object.NO_UNIT,             ### a_out/a_in for dynamical stability
            )
        )

        object.add_method(
            "a_out_div_a_in_semisecular_regime",
            (
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                object.NO_UNIT,             ### e_in
                object.NO_UNIT,             ### e_out
                object.NO_UNIT,             ### check_for_semisecular_regime_parameter
            ),
            (
                object.NO_UNIT,             ### a_out/a_in for semisecular regime
            )
        )


        object.add_method(
            "roche_radius",
            (
                unit_l,                     ### rp
                object.NO_UNIT,             ### q
                object.NO_UNIT,             ### e
                object.NO_UNIT,             ### f
                object.NO_UNIT,             ### roche_radius_specification
            ),
            (
                unit_l,                     ### Roche radius
            )
        )

        """
        object.add_method(
            "get_time",
            (),
            (units.Myr, object.ERROR_CODE,)
        )
        object.add_method(
            "set_time",
            (units.Myr, ),
            (object.ERROR_CODE,)
        )       
        object.add_method(
            "get_relative_tolerance",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_relative_tolerance",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_equations_of_motion_specification",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_equations_of_motion_specification",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_roche_radius_specification",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_roche_radius_specification",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_stability_limit_specification",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_stability_limit_specification",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_include_quadrupole_terms",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_include_quadrupole_terms",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_include_octupole_terms",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_include_octupole_terms",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )        
        object.add_method(
            "get_f_mass_tranfer",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_mass_tranfer",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )   
        object.add_method(
            "get_f_1PN_in",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_1PN_in",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )        
        object.add_method(
            "get_f_1PN_out",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_1PN_out",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_1PN_in_out",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_1PN_in_out",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_25PN_in",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_25PN_in",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_25PN_out",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_25PN_out",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        """
    def before_get_parameter(self):
        """
        Called everytime just before a parameter is retrieved in using::
            instance.parameter.name
        """
        pass
        
    def before_set_parameter(self):
        """
        Called everytime just before a parameter is updated in using::
            instance.parameter.name = newvalue
        """
        pass

    def define_particle_sets(self,object):
        object.define_inmemory_set('triples')

    def check_for_dynamical_stability(self):

        for index_triple, triple in enumerate(self.triples):
            inner_binary,outer_binary,star1,star2,star3 = give_binaries_and_stars(self,triple)
            m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out = give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary,triple)
            
            a_out_div_a_in_dynamical_stability = self.a_out_div_a_in_dynamical_stability(m1,m2,m3,a_in,a_out,e_in, e_out,triple.relative_inclination, self.parameters.stability_limit_specification)
            if a_out/a_in <= a_out_div_a_in_dynamical_stability:
                if self.parameters.verbose == True:
                    print('SecularTriple -- triple system is dynamically unstable: a_out/a_in = ',a_out/a_in,', whereas for dynamical stability, a_out/a_in should be > ',a_out_div_a_in_dynamical_stability)
                triple.dynamical_instability = True
            else:
                triple.dynamical_instability = False

    def check_for_semisecular_regime(self):

        for index_triple, triple in enumerate(self.triples):
            inner_binary,outer_binary,star1,star2,star3 = give_binaries_and_stars(self,triple)
            m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out = give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary,triple)
            
            check_for_semisecular_regime_parameter = self.parameters.check_for_semisecular_regime_parameter
            a_out_div_a_in_semisecular_regime = self.a_out_div_a_in_semisecular_regime(m1,m2,m3,e_in,e_out,check_for_semisecular_regime_parameter)
            if a_out/a_in <= a_out_div_a_in_semisecular_regime:
                if self.parameters.verbose == True:
                    print('SecularTriple -- triple system is within the semisecular regime: a_out/a_in = ',a_out/a_in,'; a_out/a_in for semisecular regime: ',a_out_div_a_in_semisecular_regime)
                triple.semisecular_regime = True
            else:
                triple.semisecular_regime = False
                


    def evolve_model(self,end_time):
        if end_time is None:
            if self.parameters.verbose == True:
                print('SecularTriple -- please specify end time!')
            return

        parameters = self.parameters
        triples = self.triples
        
        for index_triple, triple in enumerate(triples):

            ### extract data from triple ###
            self.time_step = end_time - self.model_time            

            inner_binary,outer_binary,star1,star2,star3 = give_binaries_and_stars(self,triple)
            args,skip_integration = extract_data_and_give_args(self,triple,inner_binary,outer_binary,star1,star2,star3)
            if skip_integration == True: 
                continue

            ### solve system of ODEs ###
            m1,m2,m3,R1,R2,R3, \
            spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3, \
            a_in,a_out,e_in,e_out, delta_e_in, \
            INCL_in,INCL_out,INCL_in_out, \
            AP_in,AP_out,LAN_in,LAN_out, \
            end_time_cvode,CVODE_flag,root_finding_flag = self.evolve(*args)
            if self.parameters.verbose == True:
                print('SecularTriple -- done; t/Myr = ', end_time.value_in(units.Myr),' a_in/AU=',a_in.value_in(units.AU),' a_out/AU=',a_out.value_in(units.AU),'; e_in=',e_in,'e_out=',e_out,' rel_INCL = ',INCL_in,' spin_freq1/day^{-1} = ',spin_angular_frequency1.value_in(1.0/units.day),'LAN_in',LAN_in,'LAN_out',LAN_out,'AP_in',AP_in,'AP_out',AP_out)

            print_CVODE_output(self,CVODE_flag)
            triple.error_flag_secular = CVODE_flag

            ####################
            ### root finding ###
            ####################
            triple.dynamical_instability = False
            triple.inner_collision = False
            triple.outer_collision = False
            triple.semisecular_regime = False
            
            ### The secular code will stop the integration if a root is found.
            ### In some cases, it is desirable to, nevertheless, continue integrating
            ### until the end time originally specified from triple.py.
            ### In those cases, evolve_further_until_original_end_time is to be set to `True'.

            if CVODE_flag==CV_ROOT_RETURN:
                
                ### dynamical instability ###
                if root_finding_flag==1:
                    if self.parameters.verbose == True:
                        print('SecularTriple -- triple dynamical instability')
                    triple.dynamical_instability = True
    
                ### inner collision ###
                if root_finding_flag==2: 
                    triple.inner_collision = True
                    if self.parameters.verbose == True:
                        print( 'SecularTriple --inner collision')
                ### outer collision ###
                if root_finding_flag==3: 
                    triple.outer_collision = True
                    if self.parameters.verbose == True:
                        print( 'SecularTriple --outer collision')
    
                ### RLOF star1 ###
                if root_finding_flag==4: 
                    star1.is_donor = True
                    if self.parameters.verbose == True:
                        print( 'SecularTriple -- star 1 has filled its Roche Lobe during secular integration')
                if root_finding_flag==-4: 
                    star1.is_donor = False
                    if self.parameters.verbose == True:
                        print( 'SecularTriple -- star 1 no longer fills its Roche Lobe during secular integration')
                    self.evolve_further_after_root_was_found = True
    
                ### RLOF star2 ###
                if root_finding_flag==5: 
                    star2.is_donor = True
                    if self.parameters.verbose == True:
                        print( 'SecularTriple -- star 2 has filled its Roche Lobe during secular integration')
                if root_finding_flag==-5: 
                    star2.is_donor = False
                    if self.parameters.verbose == True:
                        print( 'SecularTriple -- star 2 no longer fills its Roche Lobe during secular integration')
                    self.evolve_further_after_root_was_found = True
    
                ### RLOF star3 ###
                if root_finding_flag==6: 
                    star3.is_donor = True
                    if self.parameters.verbose == True:
                        print( 'SecularTriple -- star 3 has filled its Roche Lobe during secular integration')
                if root_finding_flag==-6: 
                    star3.is_donor = False
                    if self.parameters.verbose == True:
                        print( 'SecularTriple -- star 3 no longer fills its Roche Lobe during secular integration')
                    self.evolve_further_after_root_was_found = True


                ### semisecular regime ###
                if root_finding_flag==7:
                    if self.parameters.verbose == True:
                        print( 'SecularTriple -- semisecular regime')
                    triple.semisecular_regime = True

            ### update model time ###
            self.model_time += end_time_cvode

            ### update triple particle ###
            if parameters.include_inner_tidal_terms == True:
                star1.spin_angular_frequency = spin_angular_frequency1
                star2.spin_angular_frequency = spin_angular_frequency2
            if parameters.include_outer_tidal_terms == True:
                star3.spin_angular_frequency = spin_angular_frequency3

            if parameters.include_magnetic_braking_terms == True:
                star1.spin_angular_frequency = spin_angular_frequency1
                star2.spin_angular_frequency = spin_angular_frequency2
                star3.spin_angular_frequency = spin_angular_frequency3

            if parameters.include_spin_radius_mass_coupling_terms_star1 == True:
                star1.spin_angular_frequency = spin_angular_frequency1
            if parameters.include_spin_radius_mass_coupling_terms_star2 == True:
                star2.spin_angular_frequency = spin_angular_frequency2
            if parameters.include_spin_radius_mass_coupling_terms_star3 == True:
                star3.spin_angular_frequency = spin_angular_frequency3
                
            inner_binary.semimajor_axis = a_in
            inner_binary.eccentricity = e_in
            outer_binary.semimajor_axis = a_out
            outer_binary.eccentricity = e_out
            triple.relative_inclination = INCL_in
            inner_binary.argument_of_pericenter = AP_in
            outer_binary.argument_of_pericenter = AP_out
            inner_binary.longitude_of_ascending_node = LAN_in
            outer_binary.longitude_of_ascending_node = LAN_out
            inner_binary.delta_e_in = delta_e_in

            if self.evolve_further_after_root_was_found == True:
                original_time_step = self.time_step ### time-step given from triple.py
                old_secular_time_step = end_time_cvode ### time-step made by secular code until root was found
                new_secular_time_step = original_time_step - old_secular_time_step ### the remaining time that the secular code should integrate for
                new_end_time = self.model_time + new_secular_time_step

                if self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe == False:
                    string = "NOT"
                else:
                    string = ""
                if self.parameters.verbose == True:
                    print( 'SecularTriple -- root was found at t = ',self.model_time,'; integrating further until end time t = ',new_end_time,', specified from triple.py; RLOF is ',string,' taken into account.')
                
                ### In the remaining time, do not check for RLOF.
                ### However, for consistency with the stellar/binary evolution code, if the latter assumed RLOF from the beginning of the time-step, do take into account effect of RLOF on the orbit.
                if root_finding_flag==-4:
                    parameters.check_for_inner_RLOF = False 
                    if self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe == True:
                        star1.is_donor = True 
                if root_finding_flag==-5:
                    parameters.check_for_inner_RLOF = False
                    if self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe == True:
                        star2.is_donor = True
                if root_finding_flag==-6:
                    parameters.check_for_outer_RLOF = False
                    if self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe == True:
                        star3.is_donor = True
                self.evolve_further_after_root_was_found = False
                self.evolve_model(new_end_time)

                ### Check for RLOF again in the next call by triple.py -- this parameter must originally have been True, otherwise a root could never have been found in the first place.
                ### Set is_donor to False for next call by triple.py.
                if root_finding_flag==-4:
                    parameters.check_for_inner_RLOF = True 
                    star1.is_donor = False
                if root_finding_flag==-5:
                    parameters.check_for_inner_RLOF = True
                    star2.is_donor = False
                if root_finding_flag==-6:
                    parameters.check_for_outer_RLOF = True
                    star3.is_donor = False

    def give_roche_radii(self,triple):
        if triple is None:
            if self.parameters.verbose == True:
                print( 'SecularTriple -- please give triple particle')
            return

        inner_binary,outer_binary,star1,star2,star3 = give_binaries_and_stars(self,triple)
        m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out = give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary,triple)

        spin_angular_frequency1 = star1.spin_angular_frequency
        spin_angular_frequency2 = star2.spin_angular_frequency
        spin_angular_frequency3 = star3.spin_angular_frequency

        rp_in = a_in*(1.0-e_in)
        rp_out = a_out*(1.0-e_out)

        spin_angular_frequency_inner_orbit_periapse = numpy.sqrt( constants.G*(m1+m2)*(1.0+e_in)/(rp_in**3) ) 
        spin_angular_frequency_outer_orbit_periapse = numpy.sqrt( constants.G*(m1+m2+m3)*(1.0+e_out)/(rp_out**3) ) 

        f1 = spin_angular_frequency1/spin_angular_frequency_inner_orbit_periapse
        f2 = spin_angular_frequency2/spin_angular_frequency_inner_orbit_periapse
        f3 = spin_angular_frequency3/spin_angular_frequency_outer_orbit_periapse      

        R_L_star1 = self.roche_radius(rp_in,m1/m2,e_in,f1, self.parameters.roche_radius_specification)
        R_L_star2 = self.roche_radius(rp_in,m2/m1,e_in,f2, self.parameters.roche_radius_specification)
        R_L_star3 = self.roche_radius(rp_out,m3/(m1+m2),e_out,f3, self.parameters.roche_radius_specification)
                        
        return R_L_star1,R_L_star2,R_L_star3

    def compute_effect_of_SN_on_triple(self,Vkick_1,Vkick_2,Vkick_3,delta_m_1,delta_m_2,delta_m_3,inner_true_anomaly,outer_true_anomaly):
        triples = self.triples
        
        
        #################################################
        ### convention: new mass = old mass - delta_m ###
        ### i.e., delta_m>0 for mass loss             ###
        #################################################
        
        for index_triple, triple in enumerate(triples):
            inner_binary,outer_binary,star1,star2,star3 = give_binaries_and_stars(self,triple)
            m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out = give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary,triple)
            
            if self.parameters.verbose == True:
                print( 'PRE orb',a_in.value_in(units.AU),a_out.value_in(units.AU),e_in,e_out,INCL_in,AP_in,AP_out,LAN_in,LAN_out)
            
            ### pre-SN -- from orbital elements to orbital vectors ###
            e_in_vec_x,e_in_vec_y,e_in_vec_z,e_out_vec_x,e_out_vec_y,e_out_vec_z,h_in_vec_x,h_in_vec_y,h_in_vec_z,h_out_vec_x,h_out_vec_y,h_out_vec_z = self.compute_orbital_vectors_from_orbital_elements(m1,m2,m3,a_in,a_out,e_in,e_out,INCL_in,AP_in,AP_out,LAN_in,LAN_out)
            
            if self.parameters.verbose == True:
                print( 'PRE vec',e_in_vec_x,e_in_vec_y,e_in_vec_z,e_out_vec_x,e_out_vec_y,e_out_vec_z,h_in_vec_x,h_in_vec_y,h_in_vec_z,h_out_vec_x,h_out_vec_y,h_out_vec_z)
            
            ### effect of SN on orbital vectors ###
            V1_prime_x,V1_prime_y,V1_prime_z, \
                V2_prime_x,V2_prime_y,V2_prime_z, \
                V3_prime_x,V3_prime_y,V3_prime_z, \
                e_in_vec_x_prime,e_in_vec_y_prime,e_in_vec_z_prime,e_out_vec_x_prime,e_out_vec_y_prime,e_out_vec_z_prime, \
                h_in_vec_x_prime,h_in_vec_y_prime,h_in_vec_z_prime,h_out_vec_x_prime,h_out_vec_y_prime,h_out_vec_z_prime, \
                cos_phi1,cos_phi2,R,v_sys, \
                r1dotv1, r1dotvk, v1dotvk, \
                r1dotvsys, v1dotvsys, r2dotv2, r2dotvsys, v2dotvsys, r1dotr2, r1dotv2 = \
                self.compute_effect_of_SN_on_orbital_vectors(m1,m2,m3, \
                    e_in_vec_x,e_in_vec_y,e_in_vec_z,e_out_vec_x,e_out_vec_y,e_out_vec_z, \
                    h_in_vec_x,h_in_vec_y,h_in_vec_z,h_out_vec_x,h_out_vec_y,h_out_vec_z, \
                    Vkick_1[0],Vkick_1[1],Vkick_1[2],Vkick_2[0],Vkick_2[1],Vkick_2[2],Vkick_3[0],Vkick_3[1],Vkick_3[2],
                    delta_m_1,delta_m_2,delta_m_3,inner_true_anomaly,outer_true_anomaly)
            
            unit_v = unit_l/unit_t
            V1_prime = [V1_prime_x.value_in(unit_v),V1_prime_y.value_in(unit_v),V1_prime_z.value_in(unit_v)] | unit_v
            V2_prime = [V2_prime_x.value_in(unit_v),V2_prime_y.value_in(unit_v),V2_prime_z.value_in(unit_v)] | unit_v
            V3_prime = [V3_prime_x.value_in(unit_v),V3_prime_y.value_in(unit_v),V3_prime_z.value_in(unit_v)] | unit_v
            
            if self.parameters.verbose == True:
                print( 'POST vec',e_in_vec_x_prime,e_in_vec_y_prime,e_in_vec_z_prime,e_out_vec_x_prime,e_out_vec_y_prime,e_out_vec_z_prime, \
                h_in_vec_x_prime,h_in_vec_y_prime,h_in_vec_z_prime,h_out_vec_x_prime,h_out_vec_y_prime,h_out_vec_z_prime)
            
            ### post-SN -- from orbital vectors to orbital elements (with new reference frame) ###
            a_in_prime,a_out_prime,e_in_prime,e_out_prime,INCL_rel_prime,AP_in_prime,AP_out_prime,LAN_in_prime,LAN_out_prime = \
                self.compute_orbital_elements_from_orbital_vectors(m1-delta_m_1, m2-delta_m_2, m3-delta_m_3, \
                    e_in_vec_x_prime,e_in_vec_y_prime,e_in_vec_z_prime, \
                    e_out_vec_x_prime,e_out_vec_y_prime,e_out_vec_z_prime, \
                    h_in_vec_x_prime,h_in_vec_y_prime,h_in_vec_z_prime, \
                    h_out_vec_x_prime,h_out_vec_y_prime,h_out_vec_z_prime)
            
            if self.parameters.verbose == True:
                print( 'POST orb',a_in_prime.value_in(units.AU),a_out_prime.value_in(units.AU),e_in_prime,e_out_prime,INCL_rel_prime,AP_in_prime,AP_out_prime,LAN_in_prime,LAN_out_prime)
    
            inner_binary.semimajor_axis = a_in_prime
            inner_binary.eccentricity = e_in_prime
            outer_binary.semimajor_axis = a_out_prime
            outer_binary.eccentricity = e_out_prime
            triple.relative_inclination = INCL_rel_prime
            inner_binary.argument_of_pericenter = AP_in_prime
            outer_binary.argument_of_pericenter = AP_out_prime
            inner_binary.longitude_of_ascending_node = LAN_in_prime
            outer_binary.longitude_of_ascending_node = LAN_out_prime

            #return V1_prime,V2_prime,V3_prime,cos_phi1,cos_phi2,R,v_sys,r1dotv1, r1dotvk, v1dotvk, r1dotvsys, v1dotvsys, r2dotv2, r2dotvsys, v2dotvsys, r1dotr2, r1dotv2
            return V1_prime,V2_prime,V3_prime

def print_CVODE_output(self,CVODE_flag):
    if self.parameters.verbose == False: return
    
    if CVODE_flag==CV_SUCCESS:
        print( "SecularTriple -- ODE integration proceeded successfully.")
    elif CVODE_flag==CV_ROOT_RETURN:
        print( "SecularTriple -- root was found during ODE integration.")
    elif CVODE_flag==CV_WARNING:
        print( "SecularTriple -- (recoverable) warnings occurred during ODE integration.")
    if CVODE_flag<0:
        if CVODE_flag==CV_TOO_MUCH_WORK:
            message = "too many steps were taken during ODE integration."
        elif CVODE_flag==CV_TOO_MUCH_ACC:
            message = "required accuracy could not be obtained."
        elif CVODE_flag==CV_ERR_FAILURE:
            message = "error test failures occurred too many times during one internal time-step or minimum step size was reached."
        elif CVODE_flag==CV_CONV_FAILURE:
            message = "convergence test failures occurred too many times during one internal time-step or minimum step size was reached."
        elif CVODE_flag==CV_LINIT_FAIL:    
            message = "the linear solver's initialization function failed."
        elif CVODE_flag==CV_LSETUP_FAIL:    
            message = "the linear solver's setup function failed in an unrecoverable manner."
        elif CVODE_flag==CV_LSOLVE_FAIL:    
            message = "the linear solver's solve function failed in an unrecoverable manner."
        elif CVODE_flag==CV_RHSFUNC_FAIL:
            message = "the right-hand side function failed in an unrecoverable manner."
        elif CVODE_flag==CV_FIRST_RHSFUNC_ERR:    
            message = "the right-hand side function failed at the first call."
        elif CVODE_flag==CV_REPTD_RHSFUNC_ERR:    
            message = "the right-hand side function had repetead recoverable errors."
        elif CVODE_flag==CV_UNREC_RHSFUNC_ERR:    
            message = "the right-hand side function had a recoverable error, but no recovery is possible."
        elif CVODE_flag==CV_RTFUNC_FAIL:    
            message = "the rootfinding function failed in an unrecoverable manner."
        elif CVODE_flag==CV_MEM_FAIL:    
            message = "a memory allocation failed."
        elif CVODE_flag==CV_MEM_NULL:    
            message = "the cvode mem argument was NULL."
        elif CVODE_flag==CV_ILL_INPUT:    
            message = "one of the function inputs is illegal."
        elif CVODE_flag==CV_NO_MALLOC:    
            message = "the cvode memory block was not allocated by a call to CVodeMalloc."
        elif CVODE_flag==CV_BAD_K:    
            message = "the derivative order k is larger than the order used."
        elif CVODE_flag==CV_BAD_T:    
            message = "the time t is outside the last step taken."
        elif CVODE_flag==CV_BAD_DKY:    
            message = "the output derivative vector is NULL."
        elif CVODE_flag==CV_TOO_CLOSE:    
            message = "the output and initial times are too close to each other."
        else:
            message = "unknown error."
        print( "SecularTriple -- unrecoverable error occurred during secular integration (CVODE_flag ",str(CVODE_flag),"): ",message)


def give_binaries_and_stars(self,triple):
    ### the 'old' method ###
#    inner_binary = triple.child1
#    outer_binary = triple.child2

    ### the 'new' method -- removes the superflous 'triple' layer ###
    inner_binary = triple.child2
    outer_binary = triple

    star1 = inner_binary.child1
    star2 = inner_binary.child2
    star3 = outer_binary.child1    

    return inner_binary,outer_binary,star1,star2,star3

def give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary,triple=None):
    
    ### stellar parameters ###
    m1 = star1.mass
    m2 = star2.mass
    m3 = star3.mass

    R1 = star1.radius
    R2 = star2.radius
    R3 = star3.radius   
    
    ### binary parameters
    a_in = inner_binary.semimajor_axis
    a_out = outer_binary.semimajor_axis
    e_in = inner_binary.eccentricity
    e_out = outer_binary.eccentricity

    if triple==None:
        INCL_in = 0.0
        if self.parameters.verbose == True:
            print( 'SecularTriple -- warning: relative inclination is being set zero in method give_stellar_masses_radii_and_binary_parameters')
    else:
        INCL_in = triple.relative_inclination ### in radians
    INCL_out = 0.0
    AP_in = inner_binary.argument_of_pericenter
    AP_out = outer_binary.argument_of_pericenter
    LAN_in = inner_binary.longitude_of_ascending_node
    LAN_out = outer_binary.longitude_of_ascending_node
    
    return m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out

def extract_data_and_give_args(self,triple,inner_binary,outer_binary,star1,star2,star3):
    parameters = self.parameters
    skip_integration = False ### if True, no integration will be done in evolve function
        
    m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out = give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary,triple)
    
    ### By default, RLOF is taken into account for the remaining integration time if
    ### a root is found corresponding to the Roche lobe no longer being filled.
    ### However, in the latter case, if initially R>RL and is_donor = False, then RLOF is not taken into account during the remaining integration.
    self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe = True
    RL1,RL2,RL3 = self.give_roche_radii(triple)
    if (R1>=RL1 and star1.is_donor == False):
        if self.parameters.verbose == True:
            print( 'SecularTriple -- warning: R1>=RL1 at initialisation while star1.is_donor = False; effects of mass transfer on the orbit will not be taken into account.')
        self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe = False
    if (R2>=RL2 and star2.is_donor == False):
        if self.parameters.verbose == True:
            print( 'SecularTriple -- warning: R2>=RL2 at initialisation while star2.is_donor = False; effects of mass transfer on the orbit will not be taken into account.')
        self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe = False
    if (R3>=RL3 and star3.is_donor == False):
        if self.parameters.verbose == True:
            print( 'SecularTriple -- warning: R3>=RL3 at initialisation while star3.is_donor = False; effects of mass transfer on the orbit will not be taken into account.')
        self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe = False        

    ### if enabled, check for dynamical stability at initialisation ###
    if parameters.check_for_dynamical_stability_at_initialisation == True:
        a_out_div_a_in_dynamical_stability = self.a_out_div_a_in_dynamical_stability(m1,m2,m3,a_in,a_out,e_in,e_out,triple.relative_inclination, self.parameters.stability_limit_specification)
        if a_out/a_in <= a_out_div_a_in_dynamical_stability:
            if self.parameters.verbose == True:
                print( 'SecularTriple -- code parameter "check_for_dynamical_stability_at_initialisation" = True')
                print( 'SecularTriple -- given system is initially dynamically unstable: a_out/a_in = ',a_out/a_in,', whereas for dynamical stability, a_out/a_in should be > ',a_out_div_a_in_dynamical_stability)
                print( 'SecularTriple -- no integration will be carried out')
            skip_integration = True
            triple.dynamical_instability = True

    ### if enabled, check for semisecular regime at initialisation ###
    if parameters.check_for_semisecular_regime_at_initialisation == True:
        check_for_semisecular_regime_parameter = self.parameters.check_for_semisecular_regime_parameter
        a_out_div_a_in_semisecular_regime = self.a_out_div_a_in_semisecular_regime(m1,m2,m3,e_in,e_out,check_for_semisecular_regime_parameter)
        if a_out/a_in <= a_out_div_a_in_semisecular_regime:
            if self.parameters.verbose == True:
                print( 'SecularTriple -- code parameter "check_for_semisecular_regime_at_initialisation" = True')
                print( 'SecularTriple -- given system is initially in semisecular regime: a_out/a_in = ',a_out/a_in,'; a_out/a_in for semisecular regime: ',a_out_div_a_in_semisecular_regime)
                print( 'SecularTriple -- no integration will be carried out')
            skip_integration = True
            triple.semisecular_regime = True

    wind_mass_loss_rate_star1 = wind_mass_loss_rate_star2 = wind_mass_loss_rate_star3 = 0.0 | units.MSun/units.yr
    if parameters.include_linear_mass_change == True:
        try:
            wind_mass_loss_rate_star1 = star1.wind_mass_loss_rate
            wind_mass_loss_rate_star2 = star2.wind_mass_loss_rate
            wind_mass_loss_rate_star3 = star3.wind_mass_loss_rate
        except AttributeError:
            sys.exit('SecularTriple -- mass time_derivative_of_mass is needed for all three stars if include_linear_mass_change==True! exiting')
    time_derivative_of_radius_star1=time_derivative_of_radius_star2=time_derivative_of_radius_star3 = 0.0 | units.RSun/units.s
    if parameters.include_linear_radius_change == True:
        try:
            time_derivative_of_radius_star1 = star1.time_derivative_of_radius
            time_derivative_of_radius_star2 = star2.time_derivative_of_radius
            time_derivative_of_radius_star3 = star3.time_derivative_of_radius
        except AttributeError:
            sys.exit('SecularTriple -- time_derivative_of_radius is needed for all three stars if include_linear_radius_change==True! exiting')

    ### mass variation parameters ###
    star1_is_donor = star2_is_donor = star3_is_donor = False
    inner_mass_transfer_rate = outer_mass_transfer_rate = 0.0 | units.MSun/units.yr
    inner_accretion_efficiency_wind_child1_to_child2 = inner_accretion_efficiency_wind_child2_to_child1 = 0.0
    outer_accretion_efficiency_wind_child1_to_child2 = outer_accretion_efficiency_wind_child2_to_child1 = 0.0
    inner_accretion_efficiency_mass_transfer = outer_accretion_efficiency_mass_transfer = 0.0
    inner_specific_AM_loss_mass_transfer = outer_specific_AM_loss_mass_transfer = 0.0

    if parameters.include_inner_wind_terms == True:
        try:
            wind_mass_loss_rate_star1 = star1.wind_mass_loss_rate
            wind_mass_loss_rate_star2 = star2.wind_mass_loss_rate
            inner_accretion_efficiency_wind_child1_to_child2 = inner_binary.accretion_efficiency_wind_child1_to_child2
            inner_accretion_efficiency_wind_child2_to_child1 = inner_binary.accretion_efficiency_wind_child2_to_child1
        except AttributeError:
            sys.exit("SecularTriple -- more attributes required for inner wind! exiting")
    if parameters.include_outer_wind_terms == True:
        try:
            wind_mass_loss_rate_star3 = star3.wind_mass_loss_rate
            inner_mass_transfer_rate = inner_binary.mass_transfer_rate
            outer_accretion_efficiency_wind_child1_to_child2 = outer_binary.accretion_efficiency_wind_child1_to_child2
            outer_accretion_efficiency_wind_child2_to_child1 = outer_binary.accretion_efficiency_wind_child2_to_child1
            
        except AttributeError:
            sys.exit("SecularTriple -- more attributes required for outer wind! exiting")
    if parameters.include_inner_RLOF_terms == True:
        try:
            star1_is_donor = star1.is_donor
            star2_is_donor = star2.is_donor
            inner_mass_transfer_rate = inner_binary.mass_transfer_rate
            inner_accretion_efficiency_mass_transfer = inner_binary.accretion_efficiency_mass_transfer
            inner_specific_AM_loss_mass_transfer = inner_binary.specific_AM_loss_mass_transfer
        except AttributeError:
            sys.exit("SecularTriple -- more attributes required for inner RLOF! exiting")
    if parameters.include_outer_RLOF_terms == True:
        try:
            star3_is_donor = star3.is_donor            
            outer_mass_transfer_rate = outer_binary.mass_transfer_rate
            outer_accretion_efficiency_mass_transfer = outer_binary.accretion_efficiency_mass_transfer
            outer_specific_AM_loss_mass_transfer = outer_binary.specific_AM_loss_mass_transfer
        except AttributeError:
            sys.exit("SecularTriple -- more attributes required for outer RLOF! exiting")

    ### RLOF checks ###
    spin_angular_frequency1 = spin_angular_frequency2 = spin_angular_frequency3 = 0.0 | 1.0/units.s
    gyration_radius_star1 = gyration_radius_star2 = gyration_radius_star3 = 0.0
    if parameters.check_for_inner_RLOF == True:
        try:
            spin_angular_frequency1 = star1.spin_angular_frequency
            spin_angular_frequency2 = star2.spin_angular_frequency
        except AttributeError:
            sys.exit( "SecularTriple -- more attributes required for inner RLOF check! exiting")
    if parameters.check_for_outer_RLOF == True:
        try:
            spin_angular_frequency3 = star3.spin_angular_frequency
        except AttributeError:
            sys.exit("SecularTriple -- more attributes required for outer RLOF check! exiting")

    ### magnetic braking ###
    if parameters.include_magnetic_braking_terms == True:
        try:
            spin_angular_frequency1 = star1.spin_angular_frequency
            spin_angular_frequency2 = star2.spin_angular_frequency
            spin_angular_frequency3 = star3.spin_angular_frequency
            wind_mass_loss_rate_star1 = star1.wind_mass_loss_rate
            wind_mass_loss_rate_star2 = star2.wind_mass_loss_rate
            wind_mass_loss_rate_star3 = star3.wind_mass_loss_rate
            gyration_radius_star1 = star1.gyration_radius
            gyration_radius_star2 = star2.gyration_radius                    
            gyration_radius_star3 = star3.gyration_radius                    
        except AttributeError:
            sys.exit("SecularTriple -- more attributes required for magnetic braking terms! exiting")

    ### spin-radius-mass coupling ###
    moment_of_inertia_star1 = moment_of_inertia_star2 = moment_of_inertia_star3 = 0.0 | units.MSun*units.RSun**2
    moment_of_inertia_dot_star1 = moment_of_inertia_dot_star2 = moment_of_inertia_dot_star3 = 0.0 | units.MSun*units.RSun**2/units.yr
#    moment_of_inertia_star1 = moment_of_inertia_star2 = moment_of_inertia_star3 = 1.0e8 | units.MSun*units.RSun**2
#    moment_of_inertia_dot_star1 = moment_of_inertia_dot_star2 = moment_of_inertia_dot_star3 = 0.0 | units.MSun*units.RSun**2/units.yr
    inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2 = inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1 = 0.0
    if parameters.include_spin_radius_mass_coupling_terms_star1 == True or parameters.include_spin_radius_mass_coupling_terms_star2 == True or parameters.include_spin_radius_mass_coupling_terms_star3 == True:
        try:
            spin_angular_frequency1 = star1.spin_angular_frequency
            spin_angular_frequency2 = star2.spin_angular_frequency
            spin_angular_frequency3 = star3.spin_angular_frequency
            wind_mass_loss_rate_star1 = star1.wind_mass_loss_rate
            wind_mass_loss_rate_star2 = star2.wind_mass_loss_rate
            wind_mass_loss_rate_star3 = star3.wind_mass_loss_rate
            time_derivative_of_radius_star1 = star1.time_derivative_of_radius
            time_derivative_of_radius_star2 = star2.time_derivative_of_radius
            time_derivative_of_radius_star3 = star3.time_derivative_of_radius
            moment_of_inertia_star1 = star1.moment_of_inertia_of_star
            moment_of_inertia_star2 = star2.moment_of_inertia_of_star
            moment_of_inertia_star3 = star3.moment_of_inertia_of_star
            moment_of_inertia_dot_star1 = star1.time_derivative_of_moment_of_inertia
            moment_of_inertia_dot_star2 = star2.time_derivative_of_moment_of_inertia
            moment_of_inertia_dot_star3 = star3.time_derivative_of_moment_of_inertia
#            inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2 = inner_binary.spin_angular_momentum_wind_accretion_efficiency_child1_to_child2
#            inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1 = inner_binary.spin_angular_momentum_wind_accretion_efficiency_child2_to_child1
        except AttributeError:
            sys.exit( "SecularTriple -- more attributes required for spin-radius-mass coupling terms! exiting")


    ### tides ###
    stellar_type1 = stellar_type2 = stellar_type3 = 0 | units.stellar_type
    m1_convective_envelope = m2_convective_envelope = m3_convective_envelope = 0.0 | units.MSun
    R1_convective_envelope = R2_convective_envelope = R3_convective_envelope = 0.0 | units.RSun
    luminosity_star1 = luminosity_star2 = luminosity_star3 = 0.0 | units.LSun
    AMC_star1 = AMC_star2 = AMC_star3 = 0.0
#    k_div_T_tides_star1 = k_div_T_tides_star2 = k_div_T_tides_star3 = 0.0 | 1.0/units.s
    
    if parameters.include_inner_tidal_terms == True:
        try:
            stellar_type1 = star1.stellar_type
            stellar_type2 = star2.stellar_type
            stellar_type3 = star3.stellar_type
            m1_convective_envelope = star1.convective_envelope_mass
            m2_convective_envelope = star2.convective_envelope_mass
            m3_convective_envelope = star3.convective_envelope_mass
            R1_convective_envelope = star1.convective_envelope_radius
            R2_convective_envelope = star2.convective_envelope_radius
            R3_convective_envelope = star3.convective_envelope_radius                
            luminosity_star1 = star1.luminosity
            luminosity_star2 = star2.luminosity
            luminosity_star3 = star3.luminosity
            AMC_star1 = star1.apsidal_motion_constant
            AMC_star2 = star2.apsidal_motion_constant
            AMC_star3 = star3.apsidal_motion_constant
            gyration_radius_star1 = star1.gyration_radius
            gyration_radius_star2 = star2.gyration_radius
            gyration_radius_star3 = star3.gyration_radius
            spin_angular_frequency1 = star1.spin_angular_frequency
            spin_angular_frequency2 = star2.spin_angular_frequency
            spin_angular_frequency3 = star3.spin_angular_frequency

            ### TODO: implement checks for unphysical input values ###
            
            
        except AttributeError:
            sys.exit("SecularTriple -- more attributes required for tides! exiting")

        if (m1_convective_envelope.value_in(unit_m)<0):
            sys.exit('SecularTriple -- m1_convective_envelope must be positive! exiting')
        if (m2_convective_envelope.value_in(unit_m)<0):
            sys.exit('SecularTriple -- m2_convective_envelope must be positive! exiting')
        if (m3_convective_envelope.value_in(unit_m)<0):
            sys.exit('SecularTriple -- m3_convective_envelope must be positive! exiting')

        if (R1_convective_envelope.value_in(unit_l)<0):
            sys.exit('SecularTriple -- R1_convective_envelope must be positive! exiting')
        if (R2_convective_envelope.value_in(unit_l)<0):
            sys.exit('SecularTriple -- R2_convective_envelope must be positive! exiting')
        if (R3_convective_envelope.value_in(unit_l)<0):
            sys.exit('SecularTriple -- R3_convective_envelope must be positive! exiting')

    if self.parameters.verbose == True:
        print( 'SecularTriple -- initialization; a_in/AU=',a_in.value_in(units.AU),'; e_in=',e_in,'e_out=',e_out,' rel_INCL = ',INCL_in)

    if ((star1.is_donor == False) and (star2.is_donor == False)):
        inner_mass_transfer_rate = 0.0 | units.MSun/units.yr

    if (star3.is_donor == False):
        outer_mass_transfer_rate = 0.0 | units.MSun/units.yr

    args = [stellar_type1,stellar_type2,stellar_type3,
        m1,m2,m3,
        m1_convective_envelope,m2_convective_envelope,m3_convective_envelope,
        R1,R2,R3,
        R1_convective_envelope,R2_convective_envelope,R3_convective_envelope,
        luminosity_star1,luminosity_star2,luminosity_star3,
        spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3,
        AMC_star1,AMC_star2,AMC_star3,
        gyration_radius_star1,gyration_radius_star2,gyration_radius_star3,
        moment_of_inertia_star1,moment_of_inertia_star2,moment_of_inertia_star3,
        moment_of_inertia_dot_star1,moment_of_inertia_dot_star2,moment_of_inertia_dot_star3,
#        k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3,
        a_in,a_out,
        e_in,e_out,
        INCL_in,INCL_out,
        AP_in,AP_out,LAN_in,LAN_out,
        star1_is_donor,star2_is_donor,star3_is_donor,
        wind_mass_loss_rate_star1,wind_mass_loss_rate_star2,wind_mass_loss_rate_star3,
        time_derivative_of_radius_star1,time_derivative_of_radius_star2,time_derivative_of_radius_star3,
        inner_mass_transfer_rate,outer_mass_transfer_rate,
        inner_accretion_efficiency_wind_child1_to_child2,inner_accretion_efficiency_wind_child2_to_child1,
        outer_accretion_efficiency_wind_child1_to_child2,outer_accretion_efficiency_wind_child2_to_child1,
        inner_accretion_efficiency_mass_transfer,outer_accretion_efficiency_mass_transfer,
        inner_specific_AM_loss_mass_transfer,outer_specific_AM_loss_mass_transfer,
        inner_spin_angular_momentum_wind_accretion_efficiency_child1_to_child2,inner_spin_angular_momentum_wind_accretion_efficiency_child2_to_child1,
        self.model_time,self.time_step] ### NOTE: only time_step is used at the moment
        
    return args,skip_integration
