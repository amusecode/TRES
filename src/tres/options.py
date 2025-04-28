import numpy as np
from amuse.units import units

#--------------------------------------------------------------------------------------------------------------------
# TRES general settings
REPORT_USER_WARNINGS = True
REPORT_DEBUG = False
REPORT_DT = False 
REPORT_SN_EVOLUTION = False
REPORT_TRIPLE_EVOLUTION = False 
MAKE_PLOTS = False

REPORT_BINARY_EVOLUTION = False
REPORT_FUNCTION_NAMES = False
REPORT_MASS_TRANSFER_STABILITY = False

no_stellar_evolution = False

#reset in case of MESA by function options_mesa below
GET_GYRATION_RADIUS_FROM_STELLAR_CODE = False
GET_AMC_FROM_STELLAR_CODE = False


#--------------------------------------------------------------------------------------------------------------------
#TRES constants
time_step_factor_stable_mt = 0.01 #1% mass loss during mass transfer
# lowering this to 0.005 makes the code twice as slow
time_step_factor_ecc = 0.01
#Rl_fraction = 0.8
# 0.01 -> error in the semi-major axis of about 0.5%
maximum_wind_mass_loss_factor = 0.01 
error_dm = 0.05
#maximum_radius_change_factor = 0.005
error_dr = 0.05 #0.01
#reset in case of MESA by function options_mesa below
minimum_time_step = 1.e-9 |units.Myr


max_mass = 100 |units.MSun
min_mass = 0.08 |units.MSun # for primary stars

maximum_time_step_factor = 100.
maximum_time_step_factor_after_stable_mt = 5. 
time_step_factor_find_RLOF = 0.5
#Rl_fraction = 0.9#1.0-10.*error_dr # ratio or star radius over Roche lobe at which time step is decreased
                              # radius grows maximally by error_dr
time_step_factor_kozai = 0.025 # 0.2*0.1, 0.2-> for error in kozai timescale, 0.1-> 10 steps per cycle
kozai_type_factor = 10.
maximum_time_step = np.inf|units.Myr

kanonical_neutron_star_mass = 1.4|units.MSun
fall_back_mass = 41 |units.MSun

#--------------------------------------------------------------------------------------------------------------------
#TPS general settings
REPORT_TPS = False
REPORT_USER_WARNINGS_TPS = False
EXCLUDE_SSO = True #in order to not simulate systems with exoplanet or brown dwarf secondaries and tertiaries

#--------------------------------------------------------------------------------------------------------------------
#TPS constants
precision = 1.e-10 
absolute_max_mass = 100 |units.MSun
#  for secondaries and tertiaries
if EXCLUDE_SSO:
    absolute_min_mass = 0.0075|units.MSun 
else:
    absolute_min_mass = 0.2|units.MJupiter  

#--------------------------------------------------------------------------------------------------------------------

def options_mesa(stellar_code):   

    #suggested
#    GET_GYRATION_RADIUS_FROM_STELLAR_CODE = True
#    GET_AMC_FROM_STELLAR_CODE = True
#    minimum_time_step = 1.e-3 |units.Myr
        
    for i in range(len(stellar_code.particles)):
          
        # Opacity
        stellar_code.particles[i].set_kap('use_Type2_opacities',True)
        stellar_code.particles[i].set_kap('Zbase', 0.02)

            
        # End phase of simulation (here end of central C burning)
        stellar_code.particles[i].set_control('xa_central_lower_limit_species(1)','c12')
        stellar_code.particles[i].set_control('xa_central_lower_limit(1)',1e-12)
            
        # Winds
        # The 'Dutch' prescription combines Vink, de Jager and Nugis & Lamers
        stellar_code.particles[i].set_control('hot_wind_scheme','Dutch')
        
        # RGB winds
        stellar_code.particles[i].set_control('cool_wind_rgb_scheme','Dutch')
      
        # AGB winds    
        stellar_code.particles[i].set_control('cool_wind_agb_scheme','Dutch')
        
        # Default factor in SeBa
        stellar_code.particles[i].set_control('Dutch_scaling_factor',0.333333)

            
            
        # Overshoot  
        # Here a step overshoot is chosen. This is not the only option
        # 'exponential' overshoot is also available in MESA

        # Overshoot during central H burning 
        stellar_code.particles[i].set_control("overshoot_scheme(1)", 'step')
        stellar_code.particles[i].set_control("overshoot_zone_type(1)", 'burn_H')
        stellar_code.particles[i].set_control("overshoot_zone_loc(1)", 'core')
        stellar_code.particles[i].set_control("overshoot_bdy_loc(1)", 'top')


        # alpha_ov calibrated so that maximum radius in the MS matches that of SeBa 
        # for the 50Msun star, which corresponds to the upper limit of Pols+98 grid.
        # This value of a_ov also allows a pretty good match for all the models with M <= 50Msun.        
        stellar_code.particles[i].set_control("overshoot_f(1)",0.31)# alpha_ov/f_ov for step/exponential overshot
        stellar_code.particles[i].set_control("overshoot_f0(1)", 0.03)# insert a value strictly smaller than 10% of alpha_ov 
        
        # Overshoot during central He burning
        stellar_code.particles[i].set_control("overshoot_scheme(2)", 'step')
        stellar_code.particles[i].set_control("overshoot_zone_type(2)", 'burn_He')
        stellar_code.particles[i].set_control("overshoot_zone_loc(2)", 'core')
        stellar_code.particles[i].set_control("overshoot_bdy_loc(2)", 'top')

        stellar_code.particles[i].set_control("overshoot_f(2)",0.31)# alpha_ov/f_ov for step/exponential overshot
        stellar_code.particles[i].set_control("overshoot_f0(2)", 0.03)# insert a value strictly smaller than 10% of alpha_ov 

            
        # Semiconvection
        stellar_code.particles[i].set_control('alpha_semiconvection',0)
            
        # Semiconvection only works with Ledoux (Schwarzschild by default)
        stellar_code.particles[i].set_control('use_Ledoux_criterion', False)
            
            
        # Mixing length
        stellar_code.particles[i].set_control('mixing_length_alpha', 2.0)
            
        # Nuclear network
        stellar_code.particles[i].set_control('default_net_name','mesa_49.net')
        	
		#Time and space resolution
        stellar_code.particles[i].set_control('mesh_delta_coeff', 0.5)
        stellar_code.particles[i].set_control('time_delta_coeff', 0.4)
        stellar_code.particles[i].set_control('dH_hard_limit', 0.001)
        stellar_code.particles[i].set_control('varcontrol_target', 0.01)
            
        
        # MLT ++ ('okay_to_reduce_gradT_excess' = True)
        # We recommend to use MLT++ for massive stars (>= 10 Msun)
        # For very massive stars (>= 50 Msun), the stars usually will contract a bit at the beginning of the simulation
        # during model convergence when MLT++ is used, which may lead to TRES imposing a timestep reduction (because of too high radius_change).
        # A possible solution may be to prevent the timestep reduction in this specific case.
        stellar_code.particles[i].set_control('okay_to_reduce_gradT_excess', True)
        
        # Convergence
        stellar_code.particles[i].set_control('Pextra_factor', 2.0)
