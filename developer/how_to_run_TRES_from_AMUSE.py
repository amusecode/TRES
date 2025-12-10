# the default method to run TRES is described in README.md
# in some cases you may want to run TRES from within AMUSE
# aka to run TRES as a sort of community code or to allow for more flexibility
# this document provides a few examples of how to do this 

import numpy as np
import matplotlib.pyplot as plt
from amuse.datamodel import Particles
from amuse.units import units
from amuse.community.seba.interface import SeBa

import sys, os
sys.path.append(os.path.dirname(os.getcwd()))
import TRES as TRES
from seculartriple_TPS.interface import SecularTriple
import BIN as BIN

#simplest way of running TRES
def example_1():
    print('TRES example 1')
    tr = TRES.main()
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    # the stellar and secular codes are stopped inside main(), 
    # so no need to do it here:
    # tr.stellar_code.stop()
    # tr.secular_code.stop()


#simple way of running TRES with adjusting input parameters
#possible parameters: 
#inner_primary_mass,  inner_secondary_mass,  outer_mass, 
#inner_semimajor_axis, outer_semimajor_axis, inner_eccentricity, outer_eccentricity, 
#relative_inclination, inner_argument_of_pericenter, outer_argument_of_pericenter, inner_longitude_of_ascending_node,                      
#metallicity, 
#tend, tinit, 
#number #id number of (first) triple
#
#possible stopping conditions:
#stop_at_mass_transfer, stop_at_init_mass_transfer,stop_at_outer_mass_transfer, 
#stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer, 
#stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, 
#stop_at_no_CHE, include_CHE,
#stop_at_merger, stop_at_disintegrated,
#stop_at_inner_collision, stop_at_outer_collision,
#stop_at_dynamical_instability, 
#stop_at_semisecular_regime,  
#stop_at_SN, SN_kick_distr, 
#
#possible settings:
#impulse_kick_for_black_holes, 
#fallback_kick_for_black_holes,
#which_common_envelope,
#stop_at_CPU_time,
#max_CPU_time, 
#file_name, file_type, dir_plots

def example_2():
    print('TRES example 2')

    M1 = 1.5|units.MSun
    M2 = 0.6|units.MSun
    M3 = 0.6|units.MSun
    ain = 150|units.RSun
    aout = 15000|units.RSun
    incl = 81.0*np.pi/180.0
    
    tr = TRES.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3,
                        inner_semimajor_axis = ain, outer_semimajor_axis = aout, relative_inclination = incl)
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    # tr.stellar_code.stop()
    # tr.secular_code.stop()





#simple way of running TRES with adjusting input parameters
# evolves just the stars to 625Myr (without changing the triple), and continues the triple simulation until 630Myr
# useful in case of starting the simulation with evolved stars 
def example_3():
    print('TRES example 3')

    M1 = 2.5|units.MSun
    tinit = 625|units.Myr
    tend = 630|units.Myr
    
    tr = TRES.main(inner_primary_mass = M1, tinit = tinit, tend = tend)
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    # tr.stellar_code.stop()
    # tr.secular_code.stop()
           


# example of Kozai-Lidov cycles in a triple with plotting
# advanced level
# uses TRES.main_developer() in stead of TRES.main()
# evolve the triple to multiple timestamps
def example_4():    
    print('TRES example 4')

    M1 = 1.3|units.MSun
    M2 = 0.5|units.MSun
    M3 = 0.5|units.MSun
    Ain = 200|units.RSun
    Aout = 20000|units.RSun
    Ein = 0.1
    Eout = 0.5
    incl = 80.0*np.pi/180.0
    Gin = 0.1
    Gout = 0.5
    Oin = 0.
    metallicity = 0.02
        
    stars, bins, correct_params = TRES.make_particle_sets([M1,M2,M3], [Ain, Aout], [Ein, Eout], [incl], [Gin, Gout], [Oin])
    
    stellar_code = SeBa()
    stellar_code.parameters.metallicity = metallicity
    secular_code = SecularTriple()
    
    inner_eccentricity_array = []
    outer_eccentricity_array = []
    n_steps = 150
    time_array = (np.arange(n_steps)+1)*5|units.Myr/n_steps
    
    #make triple object and evolve for small timestep
    #needs to be bigger then 1e-4|units.Myr for secular code 
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, secular_code,
                             incl, tend=1e-4|units.Myr)
    
    for i in range(len(time_array)):
        tr.evolve_model(time_array[i])
        inner_eccentricity_array.append(tr.triple.child2.eccentricity)
        outer_eccentricity_array.append(tr.triple.eccentricity)

        if tr.check_stopping_conditions()==False or tr.check_stopping_conditions_stellar()==False or tr.check_stopping_conditions_stellar_interaction()==False:
            print('stopping conditions reached')
            time_array = time_array[:len(inner_eccentricity_array)]
            time_array[-1] = tr.triple.time
            break 

    plt.plot(time_array.value_in(units.Myr), inner_eccentricity_array, label='inner eccentricity')
    plt.plot(time_array.value_in(units.Myr), outer_eccentricity_array, label='outer eccentricity')
    plt.plot(time_array.value_in(units.Myr), inner_eccentricity_array, 'k.')
    plt.xlabel('time (Myr)')
    plt.ylabel('eccentricity')
    plt.legend(loc=0)
    plt.show()
    
    tr.stellar_code.stop()
    tr.secular_code.stop()


# example of triple with wind mass loss & calculating through mass transfer, with plotting 
# advanced level
# uses TRES.main_developer() in stead of TRES.main()
def example_5():    
    print('TRES example 5')

    M1 = 3.|units.MSun
    M2 = 0.5|units.MSun
    M3 = 0.5|units.MSun
    Ain = 825|units.RSun
    Aout = 80000|units.RSun
    Ein = 0.0
    Eout = 0.5
    incl = 0.0
    Gin = 0.1
    Gout = 0.5
    Oin = 0.
    metallicity = 0.02

    stars, bins, correct_params = TRES.make_particle_sets([M1,M2,M3], [Ain, Aout], [Ein, Eout], [incl], [Gin, Gout], [Oin])
    
    stellar_code = SeBa()
    stellar_code.parameters.metallicity = metallicity
    secular_code = SecularTriple()

    inner_semimajor_axis_array = np.array([])
    outer_semimajor_axis_array = np.array([])
    radius_primary_array = np.array([])
    stellar_type_primary_array = np.array([])
    
    time_array = (np.arange(25)+1)/25.*370
    time_array = np.append(time_array, (np.arange(100)+1)/100.*100 + 370)
    time_array = np.append(time_array, (np.arange(100)+1)/100.*10 + 470)
    time_array = time_array|units.Myr
    
    #make triple object and evolve for small timestep
    #needs to be bigger then 1e-4|units.Myr for secular code 
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, secular_code, incl, tend=1e-4|units.Myr, stop_at_mass_transfer=False)
    
    for i in range(len(time_array)):
        tr.evolve_model(time_array[i])
#        print(time_array[i], tr.triple.child2.bin_type, tr.instantaneous_evolution,tr.triple.child2.child1.stellar_type)
        inner_semimajor_axis_array = np.append(inner_semimajor_axis_array, tr.triple.child2.semimajor_axis.value_in(units.RSun))
        outer_semimajor_axis_array = np.append(outer_semimajor_axis_array, tr.triple.semimajor_axis.value_in(units.RSun))
        radius_primary_array = np.append(radius_primary_array, tr.triple.child2.child1.radius.value_in(units.RSun))
        stellar_type_primary_array = np.append(stellar_type_primary_array, tr.triple.child2.child1.stellar_type.value_in(units.stellar_type))
        
        if tr.check_stopping_conditions()==False or tr.check_stopping_conditions_stellar()==False or tr.check_stopping_conditions_stellar_interaction()==False:
            print('stopping conditions reached')
            time_array = time_array[:len(inner_semimajor_axis_array)]
            time_array[-1] = tr.triple.time
            break 
        

    plt.semilogy(time_array.value_in(units.Myr), inner_semimajor_axis_array, label='inner semimajor axis')
    plt.semilogy(time_array.value_in(units.Myr), outer_semimajor_axis_array, label='outer semimajor axis')
    plt.xlabel('time (Myr)')
    plt.ylabel('semimajor axis (RSun)')
    plt.legend(loc=0)
    plt.show()
    
    w_ms = np.arange(len(stellar_type_primary_array))[stellar_type_primary_array<=1]
    w_hg = np.arange(len(stellar_type_primary_array))[stellar_type_primary_array==2]
    w_gb = np.arange(len(stellar_type_primary_array))[stellar_type_primary_array==3]
    w_cheb = np.arange(len(stellar_type_primary_array))[stellar_type_primary_array==4]
    w_agb = np.arange(len(stellar_type_primary_array))[stellar_type_primary_array==5]
    w_wd = np.arange(len(stellar_type_primary_array))[(stellar_type_primary_array>=10) & (stellar_type_primary_array<=12)]

    plt.semilogy(time_array.value_in(units.Myr), radius_primary_array, 'k')
    plt.semilogy(time_array.value_in(units.Myr)[w_ms], radius_primary_array[w_ms], '.', label='MS')
    plt.semilogy(time_array.value_in(units.Myr)[w_hg], radius_primary_array[w_hg], '.', label='HG')
    plt.semilogy(time_array.value_in(units.Myr)[w_gb], radius_primary_array[w_gb], '.', label='GB')
    plt.semilogy(time_array.value_in(units.Myr)[w_cheb], radius_primary_array[w_cheb], '.', label='CHeB')
    plt.semilogy(time_array.value_in(units.Myr)[w_agb], radius_primary_array[w_agb], '.', label='AGB')
    plt.semilogy(time_array.value_in(units.Myr)[w_wd], radius_primary_array[w_wd], '.', label='WD')
    plt.xlabel('time (Myr)')
    plt.ylabel('primary radius (RSun)')
    plt.legend(loc=0)
    plt.show()
    

    tr.stellar_code.stop()
    tr.secular_code.stop()

       

    

                   
# advanced level
# uses TRES.main_developer() in stead of TRES.main()
# evolve the triple to 2Myr, 3Myr, 5Myr, then evolves just the stars to 8Myr (without changing the triple), and continues the triple simulation until 9Myr
def example_6():    
    print('TRES example 6')

    M1 = 1.3|units.MSun
    M2 = 0.5|units.MSun
    M3 = 0.5|units.MSun
    Ain = 200|units.RSun
    Aout = 20000|units.RSun
    Ein = 0.1
    Eout = 0.5
    incl = 80.0*np.pi/180.0
    Gin = 0.1
    Gout = 0.5
    Oin = 0.
    metallicity = 0.02
        
    stars, bins, correct_params = TRES.make_particle_sets([M1,M2,M3], [Ain, Aout], [Ein, Eout], [incl], [Gin, Gout], [Oin])
    
    stellar_code = SeBa()
    stellar_code.parameters.metallicity = metallicity
    secular_code = SecularTriple()
    
    #make triple object and evolve unil 2Myr
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, secular_code, incl, tend=2|units.Myr)    
    #continue evolution until 3 myr
    tr.evolve_model(3|units.Myr)
    #continue evolution until 5 myr
    tr.evolve_model(5|units.Myr)
    
    #only evolve the stars until 8Myr 
    stellar_code.evolve_model(8|units.Myr)
    channel_from_stellar = stellar_code.particles.new_channel_to(stars)
    channel_from_stellar.copy()
    #pickup triple evolution at 8Myr
    tr.triple.time = 8|units.Myr
    tr.secular_code.model_time = 8|units.Myr #not redundant!
        
#    continue triple evolution until 9Myr
    tr.evolve_model(9|units.Myr)
    
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)
    
    tr.stellar_code.stop()
    tr.secular_code.stop()


# advanced level
# uses TRES.main_developer() in stead of TRES.main()
# evolve the triple to 2Myr, 3Myr, 5Myr, then evolves just the stars to 8Myr (without changing the triple), and continues the triple simulation until 9Myr
#at 9Myr, some mass is removed from one star (without changing the triple), then the triple is evolved until 10 Myr
def example_7():    
    print('TRES example 7')

    M1 = 1.3|units.MSun
    M2 = 0.5|units.MSun
    M3 = 0.5|units.MSun
    Ain = 200|units.RSun
    Aout = 20000|units.RSun
    Ein = 0.1
    Eout = 0.5
    incl = 80.0*np.pi/180.0
    Gin = 0.1
    Gout = 0.5
    Oin = 0.
    metallicity = 0.02
        
    stars, bins, correct_params = TRES.make_particle_sets([M1,M2,M3], [Ain, Aout], [Ein, Eout], [incl], [Gin, Gout], [Oin])
    
    stellar_code = SeBa()
    stellar_code.parameters.metallicity = metallicity
    secular_code = SecularTriple()
    
    #make triple object and evolve unil 2Myr
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, secular_code, incl, tend=2|units.Myr)
    #continue evolution until 3 myr
    tr.evolve_model(3|units.Myr)
    #continue evolution until 5 myr
    tr.evolve_model(5|units.Myr)
    
    #only evolve the stars until 8Myr 
    stellar_code.evolve_model(8|units.Myr)
    channel_from_stellar = stellar_code.particles.new_channel_to(stars)
    channel_from_stellar.copy()
    #pickup triple evolution at 8Myr
    tr.triple.time = 8|units.Myr
    tr.secular_code.model_time = 8|units.Myr #not redundant!
        
#    continue triple evolution until 9Myr
    tr.evolve_model(9|units.Myr)
    
#    make modifications to stellar object
#    in this case remove some mass of the envelope (without effect on triple)
    donor_in_stellar_code = stars[0].as_set().get_intersecting_subset_in(stellar_code.particles)[0]
    donor_in_stellar_code.change_mass(-0.2|units.MSun, 0.|units.yr)    
    minimum_time_step = 1.e-9 |units.Myr
    stellar_code.evolve_model(tr.triple.time+minimum_time_step) #to get updated radii    
    channel_from_stellar.copy()
    
#    continue triple evolution until 10Myr
    tr.evolve_model(10|units.Myr)
    
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    tr.stellar_code.stop()
    tr.secular_code.stop()


# most advanced level
# uses TRES.main_developer() in stead of TRES.main()
# make triple, but don't evolve it.  evolves just the stars to 8Myr (without changing the triple), and continues the triple simulation until 9Myr
# useful in case of starting the simulation with evolved stars, stripped stars, or compact objects 
# note that the same can be achieved very simply with  tr = TRES.main(tinit = 8|units.Myr, tend = 9|units.Myr) - see example_3. 
# below option provides more flexibility
def example_8():    
    print('TRES example 8')

    M1 = 1.3|units.MSun
    M2 = 0.5|units.MSun
    M3 = 0.5|units.MSun
    Ain = 200|units.RSun
    Aout = 20000|units.RSun
    Ein = 0.1
    Eout = 0.5
    incl = 80.0*np.pi/180.0
    Gin = 0.1
    Gout = 0.5
    Oin = 0.
    metallicity = 0.02
        
    stars, bins, correct_params = TRES.make_particle_sets([M1,M2,M3], [Ain, Aout], [Ein, Eout], [incl], [Gin, Gout], [Oin])
    
    stellar_code = SeBa()
    stellar_code.parameters.metallicity = metallicity
    secular_code = SecularTriple()

    #make triple object (evolve for small timestep)
    #needs to be bigger then 1e-4|units.Myr for secular code 
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, secular_code, incl, tend=1e-4|units.Myr)
    
    #only evolve the stars until 8Myr 
    stellar_code.evolve_model(8|units.Myr)
    channel_from_stellar = stellar_code.particles.new_channel_to(stars)
    channel_from_stellar.copy()
    #pickup triple evolution at 8Myr
    tr.triple.time = 8|units.Myr
    tr.secular_code.model_time = 8|units.Myr #not redundant!
        
#    continue triple evolution until 9Myr
    tr.evolve_model(9|units.Myr)
    
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)
    
    tr.stellar_code.stop()
    tr.secular_code.stop()
    
    
# most advanced level
# uses TRES.main_developer() in stead of TRES.main()
# make triple, but don't evolve it.  evolves just the stars to 8Myr (without changing the triple), and continues the triple simulation until 9Myr
# useful in case of starting the simulation with evolved stars, stripped stars, or compact objects 
#at 9Myr, some mass is removed from one star (without changing the triple), then the triple is evolved until 10 Myr
def example_9():    
    print('TRES example 9')

    M1 = 1.3|units.MSun
    M2 = 0.5|units.MSun
    M3 = 0.5|units.MSun
    Ain = 200|units.RSun
    Aout = 20000|units.RSun
    Ein = 0.1
    Eout = 0.5
    incl = 80.0*np.pi/180.0
    Gin = 0.1
    Gout = 0.5
    Oin = 0.
    metallicity = 0.02
        
    stars, bins, correct_params = TRES.make_particle_sets([M1,M2,M3], [Ain, Aout], [Ein, Eout], [incl], [Gin, Gout], [Oin])
    
    stellar_code = SeBa()
    stellar_code.parameters.metallicity = metallicity
    secular_code = SecularTriple()

    #make triple object (evolve for small timestep)
    #needs to be bigger then 1e-4|units.Myr for secular code 
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, secular_code, incl, tend=1e-4|units.Myr)
    
    #only evolve the stars until 8Myr 
    stellar_code.evolve_model(8|units.Myr)
    channel_from_stellar = stellar_code.particles.new_channel_to(stars)
    channel_from_stellar.copy()
    #pickup triple evolution at 8Myr
    tr.triple.time = 8|units.Myr
    tr.secular_code.model_time = 8|units.Myr #not redundant!
        
#    continue triple evolution until 9Myr
    tr.evolve_model(9|units.Myr)

#    make modifications to stellar object
#    in this case remove some mass of the envelope (without effect on triple)
    donor_in_stellar_code = stars[0].as_set().get_intersecting_subset_in(stellar_code.particles)[0]
    donor_in_stellar_code.change_mass(-0.2|units.MSun, 0.|units.yr)    
    minimum_time_step = 1.e-9 |units.Myr
    stellar_code.evolve_model(tr.triple.time+minimum_time_step) #to get updated radii    
    channel_from_stellar.copy()

#    continue triple evolution until 10Myr
    tr.evolve_model(10|units.Myr)
    
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)
    
    tr.stellar_code.stop()
    tr.secular_code.stop()
    

def example_10():
    print('TRES example 10 - SSE')
    try:
        from amuse.community.sse.interface import SSE
        tr = TRES.main(stellar_code=SSE())
        print(tr.triple.eccentricity, tr.triple.child2.eccentricity)
    except ImportError:
        SSE = None


def example_11():
    print('TRES example 11 - MESA')
    try:
        from amuse.community.mesa.interface import MESA
        tr = TRES.main(stellar_code=MESA())
        print(tr.triple.eccentricity, tr.triple.child2.eccentricity)
    except ImportError:
        Mesa = None

def example_12():
    print('TRES example 12 - BIN')

    M1 = 1.5|units.MSun
    M2 = 0.6|units.MSun
    semi = 150|units.RSun
    
    bin = BIN.main(primary_mass = M1, secondary_mass = M2, semimajor_axis = semi)
    print(bin.triple.semimajor_axis,bin.triple.eccentricity)

    
   
    
example_1()
#example_2()
#example_3()
#example_4()
#example_5()
#example_6()
#example_7()
#example_8()
#example_9()
#example_10()
#example_11()
#example_12()
