# the default method to run TRES is described in README.md
# in some cases you may want to run TRES from within AMUSE
# aka to run TRES as a sort of community code or to allow for more flexibility
# this document provides a few examples of how to do this 

import numpy as np
from amuse.datamodel import Particles
from amuse.units import units
from amuse.community.seba.interface import SeBa
import TRES as TRES


#simplest way of running TRES
def example_1():
    tr = TRES.main()
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    tr.stellar_code.stop()
    tr.secular_code.stop()


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

    M1 = 1.5|units.MSun
    M2 = 0.6|units.MSun
    M3 = 0.6|units.MSun
    ain = 150|units.RSun
    aout = 15000|units.RSun
    incl = 81.0*np.pi/180.0
    
    tr = TRES.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3,
                        inner_semimajor_axis = ain, outer_semimajor_axis = aout, relative_inclination = incl)
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    tr.stellar_code.stop()
    tr.secular_code.stop()


#simple way of running TRES with adjusting input parameters
# evolves just the stars to 625Myr (without changing the triple), and continues the triple simulation until 630Myr
# useful in case of starting the simulation with evolved stars 
def example_3():

    M1 = 2.5|units.MSun
    tinit = 625|units.Myr
    tend = 630|units.Myr
    
    tr = TRES.main(inner_primary_mass = M1, tinit = tinit, tend = tend)
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    tr.stellar_code.stop()
    tr.secular_code.stop()
           
           
# advanced level
# uses TRES.main_developer() in stead of TRES.main()
# evolve the triple to multiple timestamps: 2Myr, 3Myr, 5Myr
def example_4():    
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
    Oout = Oin - np.pi
    metallicity = 0.02
        
    correct_params, Ein, Eout = TRES.test_initial_parameters(M1,M2,M3, Ain, Aout, Ein, Eout, incl, Gin, Gout, Oin)
    stars = TRES.make_stars(M1,M2,M3)
    bins = TRES.make_bins(stars, Ain, Aout, Ein, Eout, Gin, Gout, Oin, Oout)
    
    stellar_code = SeBa()
    #stellar_code = SeBa(redirection='none')
    stellar_code.parameters.metallicity = metallicity
    stellar_code.particles.add_particles(stars)
    
    #make triple object and evolve unil 2Myr
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, incl, tend=2|units.Myr)
    #continue evolution until 3 myr
    tr.evolve_model(3|units.Myr)
    #continue evolution until 5 myr
    tr.evolve_model(5|units.Myr)

    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)

    tr.stellar_code.stop()
    tr.secular_code.stop()

                   
# advanced level
# uses TRES.main_developer() in stead of TRES.main()
# evolve the triple to 2Myr, 3Myr, 5Myr, then evolves just the stars to 8Myr (without changing the triple), and continues the triple simulation until 9Myr
def example_5():    
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
    Oout = Oin - np.pi
    metallicity = 0.02
        
    correct_params, Ein, Eout = TRES.test_initial_parameters(M1,M2,M3, Ain, Aout, Ein, Eout, incl, Gin, Gout, Oin)
    stars = TRES.make_stars(M1,M2,M3)
    bins = TRES.make_bins(stars, Ain, Aout, Ein, Eout, Gin, Gout, Oin, Oout)
    
    stellar_code = SeBa()
    #stellar_code = SeBa(redirection='none')
    stellar_code.parameters.metallicity = metallicity
    stellar_code.particles.add_particles(stars)
    
    #make triple object and evolve unil 2Myr
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, incl, tend=2|units.Myr)    
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
        
#    continue triple evolution until 9Myr
    tr.evolve_model(9|units.Myr)
    
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)
    
    tr.stellar_code.stop()
    tr.secular_code.stop()


# advanced level
# uses TRES.main_developer() in stead of TRES.main()
# evolve the triple to 2Myr, 3Myr, 5Myr, then evolves just the stars to 8Myr (without changing the triple), and continues the triple simulation until 9Myr
#at 9Myr, some mass is removed from one star (without changing the triple), then the triple is evolved until 10 Myr
def example_6():    
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
    Oout = Oin - np.pi
    metallicity = 0.02
        
    correct_params, Ein, Eout = TRES.test_initial_parameters(M1,M2,M3, Ain, Aout, Ein, Eout, incl, Gin, Gout, Oin)
    stars = TRES.make_stars(M1,M2,M3)
    bins = TRES.make_bins(stars, Ain, Aout, Ein, Eout, Gin, Gout, Oin, Oout)
    
    stellar_code = SeBa()
    #stellar_code = SeBa(redirection='none')
    stellar_code.parameters.metallicity = metallicity
    stellar_code.particles.add_particles(stars)
    
    #make triple object and evolve unil 2Myr
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, incl, tend=2|units.Myr)
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
        
#    continue triple evolution until 9Myr
    tr.evolve_model(9|units.Myr)
    
#    make modifications to stellar object
#    in this case remove some mass of the envelope (without effect on triple)
    donor_in_stellar_code = stars[0].as_set().get_intersecting_subset_in(stellar_code.particles)[0]
    donor_in_stellar_code.change_mass(-0.2|units.MSun, 0.|units.yr)    
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
def example_7():    
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
    Oout = Oin - np.pi
    metallicity = 0.02
        
    correct_params, Ein, Eout = TRES.test_initial_parameters(M1,M2,M3, Ain, Aout, Ein, Eout, incl, Gin, Gout, Oin)
    stars = TRES.make_stars(M1,M2,M3)
    bins = TRES.make_bins(stars, Ain, Aout, Ein, Eout, Gin, Gout, Oin, Oout)
    
    stellar_code = SeBa()
    #stellar_code = SeBa(redirection='none')
    stellar_code.parameters.metallicity = metallicity
    stellar_code.particles.add_particles(stars)
    
    #make triple object (evolve for small timestep)
    #needs to be bigger then 1e-4|units.Myr for secular code 
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, incl, tend=1e-4|units.Myr)
    
    #only evolve the stars until 8Myr 
    stellar_code.evolve_model(8|units.Myr)
    channel_from_stellar = stellar_code.particles.new_channel_to(stars)
    channel_from_stellar.copy()
    #pickup triple evolution at 8Myr
    tr.triple.time = 8|units.Myr
        
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
def example_8():    
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
    Oout = Oin - np.pi
    metallicity = 0.02
        
    correct_params, Ein, Eout = TRES.test_initial_parameters(M1,M2,M3, Ain, Aout, Ein, Eout, incl, Gin, Gout, Oin)
    stars = TRES.make_stars(M1,M2,M3)
    bins = TRES.make_bins(stars, Ain, Aout, Ein, Eout, Gin, Gout, Oin, Oout)
    
    stellar_code = SeBa()
    #stellar_code = SeBa(redirection='none')
    stellar_code.parameters.metallicity = metallicity
    stellar_code.particles.add_particles(stars)
    
    #make triple object (evolve for small timestep)
    #needs to be bigger then 1e-4|units.Myr for secular code 
    tr = TRES.main_developer(stars, bins, correct_params, stellar_code, incl, tend=1e-4|units.Myr)
    
    #only evolve the stars until 8Myr 
    stellar_code.evolve_model(8|units.Myr)
    channel_from_stellar = stellar_code.particles.new_channel_to(stars)
    channel_from_stellar.copy()
    #pickup triple evolution at 8Myr
    tr.triple.time = 8|units.Myr
        
#    continue triple evolution until 9Myr
    tr.evolve_model(9|units.Myr)

#    make modifications to stellar object
#    in this case remove some mass of the envelope (without effect on triple)
    donor_in_stellar_code = stars[0].as_set().get_intersecting_subset_in(stellar_code.particles)[0]
    donor_in_stellar_code.change_mass(-0.2|units.MSun, 0.|units.yr)    
    channel_from_stellar.copy()

#    continue triple evolution until 10Myr
    tr.evolve_model(10|units.Myr)
    
    print(tr.triple.eccentricity, tr.triple.child2.eccentricity)
    
    tr.stellar_code.stop()
    tr.secular_code.stop()
    
    
    
#example_1()
#example_2()
#example_3()
#example_4()
#example_5()
#example_6()
example_7()
#example_8()

