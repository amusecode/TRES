##              Binary evolution
##              computes the evolution of a given binary
##              given any initial conditions (M, m, a, e, g, T, z).
 
import sys
import numpy as np
# from interactions import *
# from tidal_friction_constant import *
from amuse.units import units
from amuse.support.console import set_printing_strategy

try:
    from amuse.community.seba.interface import SeBa 
except ImportError:
    SeBa = None
try:
    from amuse.community.sse.interface import SSE 
except ImportError:
    SSE = None
try:
    from amuse.community.mesa.interface import Mesa 
except ImportError:
    Mesa = None
from seculartriple_TPS.interface import SecularTriple


from stellarsystem_class import StellarSystem_Class
from TRES_plotting import plot_data_container, plot_function
from TRES_setup import make_particle_sets, setup_stellar_code, make_dic_args
from TRES_options import REPORT_DEBUG, \
                         REPORT_EVOLUTION, \
                         MAKE_PLOTS, \
                         REPORT_USER_WARNINGS
from interactions import corotating_spin_angular_frequency_binary, \
                        lang_spin_angular_frequency, \
                        break_up_angular_frequency, \
                        criticial_angular_frequency_CHE



#-----
#for running BIN.py from other routines
def main(primary_mass = 1.3|units.MSun, secondary_mass = 0.5|units.MSun, 
            semimajor_axis = 1.0 |units.AU, eccentricity = 0.1, argument_of_pericenter = 0.5,
            longitude_of_ascending_node = 0.0,
            metallicity = 0.02, tend = 5.0 |units.Myr, tinit = 0.0|units.Myr, 
            number = 0, maximum_radius_change_factor = 0.005,
            stop_at_mass_transfer = True, stop_at_init_mass_transfer = True, 
            stop_at_stable_mass_transfer = True, stop_at_eccentric_stable_mass_transfer = True,
            stop_at_unstable_mass_transfer = False, stop_at_eccentric_unstable_mass_transfer = False, which_common_envelope = 2,
            stop_at_no_CHE = False, include_CHE = False, 
            stop_at_merger = True, stop_at_disintegrated = True, stop_at_collision = True, 
            stop_at_SN = False, SN_kick_distr = 2, impulse_kick_for_black_holes = True, fallback_kick_for_black_holes = True, 
            stop_at_CPU_time = False, max_CPU_time = 3600.0, file_name = "BIN.hdf", file_type = "hdf5", dir_plots = "", seed = -1,
            stellar_code=None, secular_code=None):

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    #set seed if specified, otherwise random
    if args["seed"]>=0:
        np.random.seed(args["seed"]) 

    eccentricity = float(eccentricity)
    argument_of_pericenter = float(argument_of_pericenter)
    longitude_of_ascending_node = float(longitude_of_ascending_node)
    #redundant parameters for binaries
    relative_inclination = 0
    stop_at_outer_mass_transfer = True
    stop_at_outer_collision = True
    stop_at_dynamical_instability = True
    stop_at_semisecular_regime = True
    stop_at_inner_collision = stop_at_collision
    args = make_dic_args(relative_inclination, tend, tinit, number, maximum_radius_change_factor,  
        stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
        stop_at_no_CHE, include_CHE,
        stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
        stop_at_dynamical_instability, stop_at_semisecular_regime, 
        stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes, fallback_kick_for_black_holes,
        stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots)

    stars, bin, correct_params = make_particle_sets(
        [primary_mass, secondary_mass], 
        [semimajor_axis], [eccentricity], [relative_inclination], 
        [argument_of_pericenter], [longitude_of_ascending_node])

    clean_up_stellar_code = False
    clean_up_secular_code = False

    if stellar_code is None or stellar_code.__module__.split(".")[-2]=="seba":
        stellar_code = SeBa()
    #    stellar_code = SeBa(redirection='none')
    #    stellar_code = SeBa(redirection='file', redirect_file='output_SeBa_BIN.txt')
    elif stellar_code.__module__.split(".")[-2]=="sse":                               
        stellar_code = SSE()
    elif stellar_code.__module__.split(".")[-2]=="mesa_r15140":                
        stellar_code = Mesa()
    else:
        print('No valid stellar evolution code selected. Options are SeBa (default), SSE or MESA')
        return bin_class_object # no codes initialized yet #silvia object??
    clean_up_stellar_code = True

    stellar_code.parameters.metallicity = metallicity
    if secular_code is None:
        secular_code = SecularTriple()
    #    secular_code = SecularTriple(redirection='none')
    #    secular_code = SecularTriple(redirection='file', redirect_file='output_SecularTriple_BIN.txt')
        clean_up_secular_code = True
        

    bin_class_object = StellarSystem_Class(stars, bin, correct_params, stellar_code, secular_code, args)            
    bin_class_object.secular_code.parameters.ignore_tertiary == True
    bin_class_object.secular_code.parameters.check_for_dynamical_stability = False
    bin_class_object.secular_code.parameters.check_for_outer_collision = False
    bin_class_object.secular_code.parameters.check_for_outer_RLOF = False

    if bin_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given binary are incorrect.')
        return bin_class_object # no codes initialized yet
    elif bin_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given binary at initialization.')
    elif stop_at_no_CHE == True and bin_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. No chemically homogeneous evolution at initialization')
    else:    
        bin_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(bin_class_object, dir_plots)
            bin_class_object.print_stellar_system()
            

    stellar_code.particles.remove_particles(stars)
    bin_set = bin_class_object.triple.as_set()    
    secular_code.triples.remove_particles(bin_set)         
    del stars, bin, bin_set

    if clean_up_stellar_code:
        bin_class_object.stellar_code.stop()
        if REPORT_USER_WARNINGS:
            print('Cleaning stellar evolution code')
    if clean_up_secular_code:
        bin_class_object.secular_code.stop()
        if REPORT_USER_WARNINGS:
            print('Cleaning secular code')
    
    return bin_class_object


def main_developer(stars, bin, correct_params, stellar_code, secular_code,
            metallicity = 0.02, tend = 5.0 |units.Myr, tinit = 0.0 |units.Myr, 
            number = 0, maximum_radius_change_factor = 0.005,
            stop_at_mass_transfer = True, 
            stop_at_stable_mass_transfer = True, stop_at_eccentric_stable_mass_transfer = True,
            stop_at_unstable_mass_transfer = False, stop_at_eccentric_unstable_mass_transfer = False, which_common_envelope = 2, 
            stop_at_no_CHE = False, include_CHE = False, 
            stop_at_merger = True, stop_at_disintegrated = True, stop_at_collision = True,
            stop_at_SN = False, SN_kick_distr = 2, impulse_kick_for_black_holes = True, fallback_kick_for_black_holes = True,
            stop_at_CPU_time = False, max_CPU_time = 3600.0, file_name = "BIN.hdf", file_type = "hdf5", dir_plots = "", seed = -1):

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    #set seed if specified, otherwise random
    if args["seed"]>=0:
        np.random.seed(args["seed"]) 

    bin.eccentricity = float(bin.eccentricity)
    bin.argument_of_pericenter = float(bin.argument_of_pericenter)
    bin.longitude_of_ascending_node = float(bin.longitude_of_ascending_node)   

    #redundant parameters for binaries
    relative_inclination = 0
    stop_at_outer_mass_transfer = True
    stop_at_outer_collision = True
    stop_at_dynamical_instability = True
    stop_at_semisecular_regime = True
    stop_at_inner_collision = stop_at_collision
    args = make_dic_args(relative_inclination, tend, tinit, number, maximum_radius_change_factor,  
        stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
        stop_at_no_CHE, include_CHE,
        stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
        stop_at_dynamical_instability, stop_at_semisecular_regime, 
        stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes, fallback_kick_for_black_holes,
        stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots)

    bin_class_object = StellarSystem_Class(stars, bins, correct_params, stellar_code, secular_code, args)
    bin_class_object.secular_code.parameters.ignore_tertiary == True    
    bin_class_object.secular_code.parameters.check_for_dynamical_stability = False
    bin_class_object.secular_code.parameters.check_for_outer_collision = False
    bin_class_object.secular_code.parameters.check_for_outer_RLOF = False

    if bin_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given binary are incorrect.')
        return bin_class_object # no codes initialized yet
    elif bin_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given binary at initialization.')
    elif stop_at_no_CHE == True and bin_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. No chemically homogeneous evolution at initialization')
    else:    
        bin_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(bin_class_object, dir_plots)
            bin_class_object.print_stellar_system()
   
    return bin_class_object



#-----

#-----
#for running bin.py from the commandline
def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-M", "--M1", unit=units.MSun, 
                      dest="primary_mass", type="float", default = 1.3|units.MSun,
                      help="primary mass [%default]")
    parser.add_option("-m", "--M2", unit=units.MSun, 
                      dest="secondary_mass", type="float", default = 0.5|units.MSun,
                      help="secondary mass [%default]")

    parser.add_option("-A", "-a",  unit=units.RSun,
                      dest="semimajor_axis", type="float", 
                      default = 200.0 |units.RSun,
                      help="semi major axis [%default]")
    parser.add_option("-E", "-e",
                      dest="eccentricity", type="float", default = 0.1,
                      help="eccentricity [%default]")
    parser.add_option("-G", "-g",
                      dest="argument_of_pericenter", type="float", default = 0.1,
                      help="argument of pericenter [rad] [%default]")
    parser.add_option("-O", "-o",
                      dest="longitude_of_ascending_node", type="float", default = 0.0,
                      help="longitude of ascending node [rad] [%default]")

    parser.add_option("-z", "-Z", dest="metallicity", type="float", default = 0.02,
                      help="metallicity [%default] %unit")
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 5.0 |units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("--initial_time", unit=units.Myr, 
                      dest="tinit", type="float", default = 0.0 |units.Myr,
                      help="initial time [%default] %unit")
    parser.add_option("-N", dest="number", type="int", default = 0,
                      help="number ID of system [%default]")
    parser.add_option("-s", dest="seed", type=int", default = -1,
                      help="seed (int) [%default]")
    parser.add_option("-r", dest="maximum_radius_change_factor", type="float", default = 0.01,
                      help="maximum_radius_change_factor [%default] %unit")

#    parser.add_option("--tidal", dest="tidal_terms", action="store_false", default = True, 
#                      help="tidal terms included [%default] %unit")

    parser.add_option("--no_stop_at_mass_transfer", dest="stop_at_mass_transfer", action="store_false", default = True,
                      help="stop at mass transfer [%default] %unit")
    parser.add_option("--no_stop_at_init_mass_transfer", dest="stop_at_init_mass_transfer", action="store_false", default = True,
                      help="stop if initially mass transfer[%default] %unit")
                      
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
    #0  alpha-ce + alpha-dce
    #1  gamma-ce + alpha-dce
    #2  seba style; combination of gamma-ce, alpha-ce & alpha-dce
    parser.add_option("--CE", dest="which_common_envelope",  type="int", default = 2,
                      help="which common envelope modeling [%default]")                      

    parser.add_option("--stop_at_no_CHE", dest="stop_at_no_CHE", 
                    action="store_true", default = False, help="stop if no chemically homogeneous evolution [%default] %unit")
    parser.add_option("--include_CHE", dest="include_CHE", 
                    action="store_true", default = False, help="include chemically homogeneous evolution in the stellar evolution [%default] %unit")

    parser.add_option("--no_stop_at_merger", dest="stop_at_merger", action="store_false", default = True, 
                      help="stop at merger [%default] %unit")
    parser.add_option("--no_stop_at_disintegrated", dest="stop_at_disintegrated", action="store_false", default = True,
                      help="stop at disintegrated [%default] %unit")
    parser.add_option("--no_stop_at_collision", dest="stop_at_inner_collision", action="store_false",default = True,
                      help="stop at collision in binary[%default] %unit")

    parser.add_option("--stop_at_SN", dest="stop_at_SN", action="store_true", default = False,
                      help="stop at supernova [%default] %unit")
    #0  No kick 
    #1  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
    #2  Arzoumanian ea 2002, 568, 289
    #3  Hansen & Phinney 1997, 291, 569
    #4  Paczynski 1990, 348, 485
    #5  Verbunt, Igoshev & Cator, 2017, 608, 57
    parser.add_option("--SN_kick_distr", dest="SN_kick_distr",  type="int", default = 5,
                      help="which supernova kick distribution [%default]")                      
    parser.add_option("--no_impulse_kick_for_black_holes", dest="impulse_kick_for_black_holes",  action="store_false", default = True,
                      help="do not rescale the BH SN kick by mass -> impulse kick [%default]")                      
    parser.add_option("--no_fallback_kick_for_black_holes", dest="fallback_kick_for_black_holes",  action="store_false", default = True,
                      help="do not rescale the BH SN kick with fallback  [%default]")                      
                      
                      
    parser.add_option("--stop_at_CPU_time", dest="stop_at_CPU_time", action="store_true", default = False,
                      help="stop at CPU time [%default] %unit")
    parser.add_option("--max_CPU_time", dest="max_CPU_time", type="float", default = 3600.0,
                      help="max CPU time [%default] %unit")

    parser.add_option("--stellar_evolution_code", dest="SE_code",  type="int", default = 0,
                      help="which stellar evolution [%default]")                                            
                      
    parser.add_option("-f", dest="file_name", type ="string", default = "BIN.hdf",#"BIN.txt"
                      help="file name[%default]")
    parser.add_option("-F", dest="file_type", type ="string", default = "hdf5",#"txt"
                      help="file type[%default]")
    parser.add_option("--dir_plots", dest="dir_plots", type ="string", default = "",#"txt"
                      help="directory for plots for debugging mode [%default]")
                                           
    options, args = parser.parse_args()
    return options.__dict__


if __name__ == '__main__':
    args = parse_arguments()

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    #set seed if specified, otherwise random
    if args["seed"]>=0:
        np.random.seed(args["seed"]) 

    #redundant parameters for binaries
    args["relative_inclination"] = 0.
    args["stop_at_outer_mass_transfer"] = True
    args["stop_at_outer_collision"] = True
    args["stop_at_dynamical_instability"] = True
    args["stop_at_semisecular_regime"] = True

    stars, bins, correct_params = make_particle_sets(
        [args["primary_mass"], args["secondary_mass"]],
        [args["semimajor_axis"]],
        [args["eccentricity"]],
        [args["relative_inclination"]],
        [args["argument_of_pericenter"]],
        [args["longitude_of_ascending_node"]])


    if args["SE_code"] == 1:
        stellar_code = SSE()
    elif args["SE_code"] == 2:
        stellar_code = Mesa()
    else:
        stellar_code = SeBa()    
#        stellar_code = SeBa(redirection='none')
#        stellar_code = SeBa(redirection='file', redirect_file='output_SeBa_BIN.txt')
    stellar_code.parameters.metallicity = args["metallicity"]

    secular_code = SecularTriple()
#    secular_code = SecularTriple(redirection='none')
#    secular_code = SecularTriple(redirection='file', redirect_file='output_SecularTriple_BIN.txt')

    bin_class_object = StellarSystem_Class(stars, bins, correct_params, stellar_code, secular_code, args)
    bin_class_object.secular_code.parameters.ignore_tertiary == True
    bin_class_object.secular_code.parameters.check_for_dynamical_stability = False
    bin_class_object.secular_code.parameters.check_for_outer_collision = False
    bin_class_object.secular_code.parameters.check_for_outer_RLOF = False
    
    if bin_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given binary are incorrect.' )   
        # no codes initialized yet
        sys.exit('Choose a different system. The parameters of the given binary are incorrect.')
    elif bin_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given bin at initialization.')
    elif args["stop_at_no_CHE"] == True and bin_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. No chemically homogeneous evolution at initialization')
    else:    
        bin_class_object.evolve_model(args["tend"])
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(bin_class_object, args['dir_plots'])
            bin_class_object.print_stellar_system()



        if REPORT_EVOLUTION:
            print('Simulation has finished succesfully')
            
    print('\nYou have used the TRES/BIN evolution code. Literature reference:')
    print('** Toonen, Hamers & Portegies Zwart 2016, ComAC, 3, 6T:')
    print('... "The evolution of hierarchical star-systems" ')
            
    bin_class_object.stellar_code.stop()
    bin_class_object.secular_code.stop()

