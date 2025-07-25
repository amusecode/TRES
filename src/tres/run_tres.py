# to do
# minus sign in mean anomaly

## Triple:      Triple evolution
##              computes the evolution of a given triple
##              given any initial conditions (M, m, l, A, a, E, e, i, G, g, O, o, T, z).

import sys
import argparse
import numpy as np

# from interactions import *
# from tidal_friction_constant import *
from amuse.units import units
from amuse.support.console import set_printing_strategy

try:
    from amuse.community.seba import Seba
except ImportError:
    Seba = None
try:
    from amuse.community.sse import Sse
except ImportError:
    Sse = None
try:
    from amuse.community.mesa import Mesa
except ImportError:
    Mesa = None
from tres.seculartriple import Seculartriple


from tres.triple_class import Triple_Class
from tres.plotting import plot_data_container, plot_function
from tres.setup import make_particle_sets, setup_stellar_code
from tres.options import REPORT_DEBUG, REPORT_TRIPLE_EVOLUTION, MAKE_PLOTS, REPORT_USER_WARNINGS
from tres.interactions import (
    corotating_spin_angular_frequency_binary,
    lang_spin_angular_frequency,
    break_up_angular_frequency,
    criticial_angular_frequency_CHE,
)


def initialize_triple_class(
    stars,
    bins,
    correct_params,
    stellar_code,
    secular_code,
    relative_inclination=80.0 * np.pi / 180.0,
    metallicity=0.02,
    tend=5.0 | units.Myr,
    tinit=0.0 | units.Myr,
    number=0,
    maximum_radius_change_factor=0.005,
    stop_at_mass_transfer=True,
    stop_at_init_mass_transfer=True,
    stop_at_outer_mass_transfer=True,
    stop_at_stable_mass_transfer=True,
    stop_at_eccentric_stable_mass_transfer=True,
    stop_at_unstable_mass_transfer=False,
    stop_at_eccentric_unstable_mass_transfer=False,
    which_common_envelope=2,
    stop_at_no_CHE=False,
    include_CHE=False,
    stop_at_merger=True,
    stop_at_disintegrated=True,
    stop_at_inner_collision=True,
    stop_at_outer_collision=True,
    stop_at_dynamical_instability=True,
    stop_at_semisecular_regime=False,
    stop_at_SN=False,
    SN_kick_distr=2,
    impulse_kick_for_black_holes=True,
    fallback_kick_for_black_holes=True,
    stop_at_CPU_time=False,
    max_CPU_time=3600.0,
    file_name="TRES.amuse",
    file_type="amuse",
    dir_plots="",
):

    triple = Triple_Class(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        relative_inclination,
        tend,
        tinit,
        number,
        maximum_radius_change_factor,
        stop_at_mass_transfer,
        stop_at_init_mass_transfer,
        stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer,
        stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer,
        stop_at_eccentric_unstable_mass_transfer,
        which_common_envelope,
        stop_at_no_CHE,
        include_CHE,
        stop_at_merger,
        stop_at_disintegrated,
        stop_at_inner_collision,
        stop_at_outer_collision,
        stop_at_dynamical_instability,
        stop_at_semisecular_regime,
        stop_at_SN,
        SN_kick_distr,
        impulse_kick_for_black_holes,
        fallback_kick_for_black_holes,
        stop_at_CPU_time,
        max_CPU_time,
        file_name,
        file_type,
        dir_plots,
    )
    triple.stellar_code.parameters.metallicity = metallicity

    return triple


# -----
# for running TRES.py from other routines
def run_tres(
    inner_primary_mass=1.3 | units.MSun,
    inner_secondary_mass=0.5 | units.MSun,
    outer_mass=0.5 | units.MSun,
    inner_semimajor_axis=1.0 | units.AU,
    outer_semimajor_axis=100.0 | units.AU,
    inner_eccentricity=0.1,
    outer_eccentricity=0.5,
    relative_inclination=80.0 * np.pi / 180.0,
    inner_argument_of_pericenter=0.1,
    outer_argument_of_pericenter=0.5,
    inner_longitude_of_ascending_node=0.0,
    metallicity=0.02,
    tend=5.0 | units.Myr,
    tinit=0.0 | units.Myr,
    number=0,
    maximum_radius_change_factor=0.005,
    stop_at_mass_transfer=True,
    stop_at_init_mass_transfer=True,
    stop_at_outer_mass_transfer=True,
    stop_at_stable_mass_transfer=True,
    stop_at_eccentric_stable_mass_transfer=True,
    stop_at_unstable_mass_transfer=False,
    stop_at_eccentric_unstable_mass_transfer=False,
    which_common_envelope=2,
    stop_at_no_CHE=False,
    include_CHE=False,
    stop_at_merger=True,
    stop_at_disintegrated=True,
    stop_at_inner_collision=True,
    stop_at_outer_collision=True,
    stop_at_dynamical_instability=True,
    stop_at_semisecular_regime=False,
    stop_at_SN=False,
    SN_kick_distr=2,
    impulse_kick_for_black_holes=True,
    fallback_kick_for_black_holes=True,
    stop_at_CPU_time=False,
    max_CPU_time=3600.0,
    file_name="TRES.amuse",
    file_type="amuse",
    dir_plots="",
    stellar_code=None,
    secular_code=None,
):

    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=11,
        prefix="",
        separator=" [",
        suffix="]",
    )

    inner_eccentricity = float(inner_eccentricity)
    outer_eccentricity = float(outer_eccentricity)
    relative_inclination = float(relative_inclination)
    inner_argument_of_pericenter = float(inner_argument_of_pericenter)
    outer_argument_of_pericenter = float(outer_argument_of_pericenter)
    inner_longitude_of_ascending_node = float(inner_longitude_of_ascending_node)

    stars, bins, correct_params = make_particle_sets(
        inner_primary_mass,
        inner_secondary_mass,
        outer_mass,
        inner_semimajor_axis,
        outer_semimajor_axis,
        inner_eccentricity,
        outer_eccentricity,
        relative_inclination,
        inner_argument_of_pericenter,
        outer_argument_of_pericenter,
        inner_longitude_of_ascending_node,
    )

    clean_up_stellar_code = False
    clean_up_secular_code = False

    if stellar_code is None or stellar_code.__module__.split(".")[-2].lower() == "seba":
        stellar_code = Seba()
    #    stellar_code = Seba(redirection='none')
    #    stellar_code = Seba(redirection='file', redirect_file='output_SeBa_TRES.txt')
    elif stellar_code.__module__.split(".")[-2].lower() == "sse":
        stellar_code = Sse()
    elif stellar_code.__module__.split(".")[-2].lower() == "mesa_r15140":
        stellar_code = Mesa()
    else:
        print("No valid stellar evolution code selected. Options are SeBa (default), SSE or MESA")
        return triple_class_object  # no codes initialized yet
    clean_up_stellar_code = True

    stellar_code.parameters.metallicity = metallicity
    if secular_code is None:
        secular_code = Seculartriple()
        #    secular_code = Seculartriple(redirection='none')
        #    secular_code = Seculartriple(redirection='file', redirect_file='output_SecularTriple_TRES.txt')
        clean_up_secular_code = True

    triple_class_object = Triple_Class(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        relative_inclination,
        tend,
        tinit,
        number,
        maximum_radius_change_factor,
        stop_at_mass_transfer,
        stop_at_init_mass_transfer,
        stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer,
        stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer,
        stop_at_eccentric_unstable_mass_transfer,
        which_common_envelope,
        stop_at_no_CHE,
        include_CHE,
        stop_at_merger,
        stop_at_disintegrated,
        stop_at_inner_collision,
        stop_at_outer_collision,
        stop_at_dynamical_instability,
        stop_at_semisecular_regime,
        stop_at_SN,
        SN_kick_distr,
        impulse_kick_for_black_holes,
        fallback_kick_for_black_holes,
        stop_at_CPU_time,
        max_CPU_time,
        file_name,
        file_type,
        dir_plots,
    )

    if triple_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The parameters of the given triple are incorrect.")
        return triple_class_object  # no codes initialized yet
    elif stop_at_semisecular_regime == True and triple_class_object.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The given triple is in the semisecular regime at initialization.")
    elif triple_class_object.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The given triple is dynamically unstable at initialization.")
    elif triple_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. There is mass transfer in the given triple at initialization.")
    elif stop_at_no_CHE == True and triple_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. No chemically homogeneous evolution at initialization")
    else:
        triple_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, dir_plots)
            triple_class_object.print_stellar_system()

    stellar_code.particles.remove_particles(stars)
    triple_set = triple_class_object.triple.as_set()
    secular_code.triples.remove_particles(triple_set)
    del stars, bins, triple_set

    if clean_up_stellar_code:
        triple_class_object.stellar_code.stop()
        if REPORT_USER_WARNINGS:
            print("Cleaning stellar evolution code")
    if clean_up_secular_code:
        triple_class_object.secular_code.stop()
        if REPORT_USER_WARNINGS:
            print("Cleaning secular code")

    return triple_class_object


def run_tres_developer(
    stars,
    bins,
    correct_params,
    stellar_code,
    secular_code,
    relative_inclination=80.0 * np.pi / 180.0,
    metallicity=0.02,
    tend=5.0 | units.Myr,
    tinit=0.0 | units.Myr,
    number=0,
    maximum_radius_change_factor=0.005,
    stop_at_mass_transfer=True,
    stop_at_init_mass_transfer=True,
    stop_at_outer_mass_transfer=True,
    stop_at_stable_mass_transfer=True,
    stop_at_eccentric_stable_mass_transfer=True,
    stop_at_unstable_mass_transfer=False,
    stop_at_eccentric_unstable_mass_transfer=False,
    which_common_envelope=2,
    stop_at_no_CHE=False,
    include_CHE=False,
    stop_at_merger=True,
    stop_at_disintegrated=True,
    stop_at_inner_collision=True,
    stop_at_outer_collision=True,
    stop_at_dynamical_instability=True,
    stop_at_semisecular_regime=False,
    stop_at_SN=False,
    SN_kick_distr=2,
    impulse_kick_for_black_holes=True,
    fallback_kick_for_black_holes=True,
    stop_at_CPU_time=False,
    max_CPU_time=3600.0,
    file_name="TRES.amuse",
    file_type="amuse",
    dir_plots="",
):

    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=11,
        prefix="",
        separator=" [",
        suffix="]",
    )

    bins.eccentricity[0] = float(bins.eccentricity[0])
    bins.eccentricity[1] = float(bins.eccentricity[1])
    bins.argument_of_pericenter[0] = float(bins.argument_of_pericenter[0])
    bins.argument_of_pericenter[1] = float(bins.argument_of_pericenter[1])
    bins.longitude_of_ascending_node[0] = float(bins.longitude_of_ascending_node[0])
    bins.longitude_of_ascending_node[1] = float(bins.longitude_of_ascending_node[1])
    relative_inclination = float(relative_inclination)

    triple_class_object = Triple_Class(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        relative_inclination,
        tend,
        tinit,
        number,
        maximum_radius_change_factor,
        stop_at_mass_transfer,
        stop_at_init_mass_transfer,
        stop_at_outer_mass_transfer,
        stop_at_stable_mass_transfer,
        stop_at_eccentric_stable_mass_transfer,
        stop_at_unstable_mass_transfer,
        stop_at_eccentric_unstable_mass_transfer,
        which_common_envelope,
        stop_at_no_CHE,
        include_CHE,
        stop_at_merger,
        stop_at_disintegrated,
        stop_at_inner_collision,
        stop_at_outer_collision,
        stop_at_dynamical_instability,
        stop_at_semisecular_regime,
        stop_at_SN,
        SN_kick_distr,
        impulse_kick_for_black_holes,
        fallback_kick_for_black_holes,
        stop_at_CPU_time,
        max_CPU_time,
        file_name,
        file_type,
        dir_plots,
    )

    if triple_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The parameters of the given triple are incorrect.")
        return triple_class_object  # no codes initialized yet
    elif stop_at_semisecular_regime == True and triple_class_object.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The given triple is in the semisecular regime at initialization.")
    elif triple_class_object.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The given triple is dynamically unstable at initialization.")
    elif triple_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. There is mass transfer in the given triple at initialization.")
    elif stop_at_no_CHE == True and triple_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. No chemically homogeneous evolution at initialization")
    else:
        triple_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, dir_plots)
            triple_class_object.print_stellar_system()

    return triple_class_object


# -----


# -----
# for running triple.py from the commandline
def parse_arguments():
    """
    Parse command line arguments and show default values.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-M",
        "--M1",
        type=units.MSun,
        dest="inner_primary_mass",
        default=1.3 | units.MSun,
        help="inner primary mass",
    )
    parser.add_argument(
        "-m",
        "--M2",
        type=units.MSun,
        dest="inner_secondary_mass",
        default=0.5 | units.MSun,
        help="inner secondary mass",
    )
    parser.add_argument(
        "-l",
        "--M3",
        type=units.MSun,
        dest="outer_mass",
        default=0.5 | units.MSun,
        help="outer mass",
    )

    parser.add_argument(
        "-A",
        "--Ain",
        type=units.RSun,
        dest="inner_semimajor_axis",
        default=200.0 | units.RSun,
        help="inner semi major axis",
    )
    parser.add_argument(
        "-a",
        "--Aout",
        type=units.RSun,
        dest="outer_semimajor_axis",
        default=20000.0 | units.RSun,
        help="outer semi major axis",
    )
    parser.add_argument(
        "-E", "--Ein", dest="inner_eccentricity", type=float, default=0.1, help="inner eccentricity"
    )
    parser.add_argument(
        "-e", "--Eout", dest="outer_eccentricity", type=float, default=0.5, help="outer eccentricity"
    )
    parser.add_argument(
        "-i",
        "-I",
        dest="relative_inclination",
        type=float,
        default=80.0 * np.pi / 180.0,
        help="relative inclination [rad]",
    )
    parser.add_argument(
        "-G",
        "--Gin",
        dest="inner_argument_of_pericenter",
        type=float,
        default=0.1,
        help="inner argument of pericenter [rad]",
    )
    parser.add_argument(
        "-g",
        "--Gout",
        dest="outer_argument_of_pericenter",
        type=float,
        default=0.5,
        help="outer argument of pericenter [rad]",
    )
    parser.add_argument(
        "-O",
        "--Oin",
        dest="inner_longitude_of_ascending_node",
        type=float,
        default=0.0,
        help="inner longitude of ascending node [rad]",
    )
    ##             outer longitude of ascending nodes = inner - pi
    #    parser.add_argument("-o",
    #                      dest="outer_longitude_of_ascending_node", type=float, default = 0.0,
    #                      help="outer longitude of ascending node [rad]")

    parser.add_argument("-z", "-Z", dest="metallicity", type=float, default=0.02, help="metallicity")
    parser.add_argument(
        "-t", "-T", type=units.Myr, dest="tend", default=5.0 | units.Myr, help="end time"
    )
    parser.add_argument(
        "--initial_time",
        type=units.Myr,
        dest="tinit",
        default=0.0 | units.Myr,
        help="initial time",
    )
    parser.add_argument("-N", dest="number", type=int, default=0, help="number ID of system")
    parser.add_argument(
        "-r",
        dest="maximum_radius_change_factor",
        type=float,
        default=0.01,
        help="maximum_radius_change_factor",
    )

    #    parser.add_argument("--tidal", dest="tidal_terms", action="store_false", default = True,
    #                      help="tidal terms included")

    parser.add_argument(
        "--no_stop_at_mass_transfer",
        dest="stop_at_mass_transfer",
        action="store_false",
        default=True,
        help="stop at mass transfer",
    )
    parser.add_argument(
        "--no_stop_at_init_mass_transfer",
        dest="stop_at_init_mass_transfer",
        action="store_false",
        default=True,
        help="stop if initially mass transfer",
    )
    parser.add_argument(
        "--no_stop_at_outer_mass_transfer",
        dest="stop_at_outer_mass_transfer",
        action="store_false",
        default=True,
        help="stop at triple mass transfer",
    )

    #   if stop_at_mass_transfer is False, the following 4 stopping conditions can be used to further specify.
    #   if stop_at_mass_transfer is True, the following 4 are ignored.
    parser.add_argument(
        "--stop_at_stable_mass_transfer",
        dest="stop_at_stable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at stable mass transfer",
    )
    parser.add_argument(
        "--stop_at_eccentric_stable_mass_transfer",
        dest="stop_at_eccentric_stable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at eccentric stable mass transfer",
    )
    # unstable mass transfer leads to common-envelope evolution
    parser.add_argument(
        "--stop_at_unstable_mass_transfer",
        dest="stop_at_unstable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at unstable mass transfer",
    )
    parser.add_argument(
        "--stop_at_eccentric_unstable_mass_transfer",
        dest="stop_at_eccentric_unstable_mass_transfer",
        action="store_true",
        default=False,
        help="stop at eccentric unstable mass transfer",
    )
    # 0  alpha-ce + alpha-dce
    # 1  gamma-ce + alpha-dce
    # 2  seba style; combination of gamma-ce, alpha-ce & alpha-dce
    parser.add_argument(
        "--CE", dest="which_common_envelope", type=int, default=2, help="which common envelope modeling"
    )

    parser.add_argument(
        "--stop_at_no_CHE",
        dest="stop_at_no_CHE",
        action="store_true",
        default=False,
        help="stop if no chemically homogeneous evolution",
    )
    parser.add_argument(
        "--include_CHE",
        dest="include_CHE",
        action="store_true",
        default=False,
        help="include chemically homogeneous evolution in the stellar evolution",
    )

    parser.add_argument(
        "--no_stop_at_merger",
        dest="stop_at_merger",
        action="store_false",
        default=True,
        help="stop at merger",
    )
    parser.add_argument(
        "--no_stop_at_disintegrated",
        dest="stop_at_disintegrated",
        action="store_false",
        default=True,
        help="stop at disintegrated",
    )
    parser.add_argument(
        "--no_stop_at_inner_collision",
        dest="stop_at_inner_collision",
        action="store_false",
        default=True,
        help="stop at collision in inner binary",
    )
    parser.add_argument(
        "--no_stop_at_outer_collision",
        dest="stop_at_outer_collision",
        action="store_false",
        default=True,
        help="stop at collision in outer binary",
    )
    parser.add_argument(
        "--no_stop_at_dynamical_instability",
        dest="stop_at_dynamical_instability",
        action="store_false",
        default=True,
        help="stop at dynamical instability",
    )
    parser.add_argument(
        "--stop_at_semisecular_regime",
        dest="stop_at_semisecular_regime",
        action="store_true",
        default=False,
        help="stop at semisecular regime",
    )

    parser.add_argument(
        "--stop_at_SN", dest="stop_at_SN", action="store_true", default=False, help="stop at supernova"
    )
    # 0  No kick
    # 1  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974
    # 2  Arzoumanian ea 2002, 568, 289
    # 3  Hansen & Phinney 1997, 291, 569
    # 4  Paczynski 1990, 348, 485
    # 5  Verbunt, Igoshev & Cator, 2017, 608, 57
    parser.add_argument(
        "--SN_kick_distr",
        dest="SN_kick_distr",
        type=int,
        default=5,
        help="which supernova kick distribution",
    )
    parser.add_argument(
        "--no_impulse_kick_for_black_holes",
        dest="impulse_kick_for_black_holes",
        action="store_false",
        default=True,
        help="do not rescale the BH SN kick by mass -> impulse kick",
    )
    parser.add_argument(
        "--no_fallback_kick_for_black_holes",
        dest="fallback_kick_for_black_holes",
        action="store_false",
        default=True,
        help="do not rescale the BH SN kick with fallback ",
    )

    parser.add_argument(
        "--stop_at_CPU_time",
        dest="stop_at_CPU_time",
        action="store_true",
        default=False,
        help="stop at CPU time",
    )
    parser.add_argument(
        "--max_CPU_time", dest="max_CPU_time", type=float, default=3600.0, help="max CPU time"
    )

    parser.add_argument(
        "--stellar_evolution_code", dest="SE_code", type=int, default=0, help="which stellar evolution"
    )

    parser.add_argument(
        "-f", dest="file_name", type=str, default="TRES.amuse", help="file name"  # "TRES.txt"
    )
    parser.add_argument("-F", dest="file_type", type=str, default="amuse", help="file type")  # "txt"
    parser.add_argument(
        "--dir_plots",
        dest="dir_plots",
        type=str,
        default="",  # "txt"
        help="directory for plots for debugging mode",
    )

    arguments = parser.parse_args()
    return arguments.__dict__


if __name__ == "__main__":
    args = parse_arguments()

    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=11,
        prefix="",
        separator=" [",
        suffix="]",
    )

    stars, bins, correct_params = make_particle_sets(
        args["inner_primary_mass"],
        args["inner_secondary_mass"],
        args["outer_mass"],
        args["inner_semimajor_axis"],
        args["outer_semimajor_axis"],
        args["inner_eccentricity"],
        args["outer_eccentricity"],
        args["relative_inclination"],
        args["inner_argument_of_pericenter"],
        args["outer_argument_of_pericenter"],
        args["inner_longitude_of_ascending_node"],
    )

    if args["SE_code"] == 1:
        stellar_code = SSE()
    elif args["SE_code"] == 2:
        stellar_code = Mesa()
    else:
        stellar_code = SeBa()
    #        stellar_code = SeBa(redirection='none')
    #        stellar_code = SeBa(redirection='file', redirect_file='output_SeBa_TRES.txt')
    stellar_code.parameters.metallicity = args["metallicity"]

    secular_code = Seculartriple()
    #    secular_code = Seculartriple(redirection='none')
    #    secular_code = Seculartriple(redirection='file', redirect_file='output_SecularTriple_TRES.txt')

    triple_class_object = Triple_Class(
        stars,
        bins,
        correct_params,
        stellar_code,
        secular_code,
        args["relative_inclination"],
        args["tend"],
        args["tinit"],
        args["number"],
        args["maximum_radius_change_factor"],
        args["stop_at_mass_transfer"],
        args["stop_at_init_mass_transfer"],
        args["stop_at_outer_mass_transfer"],
        args["stop_at_stable_mass_transfer"],
        args["stop_at_eccentric_stable_mass_transfer"],
        args["stop_at_unstable_mass_transfer"],
        args["stop_at_eccentric_unstable_mass_transfer"],
        args["which_common_envelope"],
        args["stop_at_no_CHE"],
        args["include_CHE"],
        args["stop_at_merger"],
        args["stop_at_disintegrated"],
        args["stop_at_inner_collision"],
        args["stop_at_outer_collision"],
        args["stop_at_dynamical_instability"],
        args["stop_at_semisecular_regime"],
        args["stop_at_SN"],
        args["SN_kick_distr"],
        args["impulse_kick_for_black_holes"],
        args["fallback_kick_for_black_holes"],
        args["stop_at_CPU_time"],
        args["max_CPU_time"],
        args["file_name"],
        args["file_type"],
        args["dir_plots"],
    )

    if triple_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The parameters of the given triple are incorrect.")
        # no codes initialized yet
        sys.exit("Choose a different system. The parameters of the given triple are incorrect.")
    elif args["stop_at_semisecular_regime"] == True and triple_class_object.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The given triple is in the semisecular regime at initialization.")
    elif triple_class_object.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. The given triple is dynamically unstable at initialization.")
    elif triple_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. There is mass transfer in the given triple at initialization.")
    elif args["stop_at_no_CHE"] == True and triple_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print("Choose a different system. No chemically homogeneous evolution at initialization")
    else:
        triple_class_object.evolve_model(args["tend"])
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, args["dir_plots"])
            triple_class_object.print_stellar_system()

        if REPORT_TRIPLE_EVOLUTION:
            print("Simulation has finished succesfully")

    print("\nYou have used the TRES triple evolution code. Literature reference:")
    print("** Toonen, Hamers & Portegies Zwart 2016, ComAC, 3, 6T:")
    print('... "The evolution of hierarchical triple star-systems" ')

    triple_class_object.stellar_code.stop()
    triple_class_object.secular_code.stop()
