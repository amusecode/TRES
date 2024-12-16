from amuse.units import units, constants
from amuse.io import write_set_to_file
from amuse.io import read_set_from_file
from amuse.support.console import set_printing_strategy
import numpy as np
import io
import sys
import pandas as pd #only necessary when saving as csv
import os
    
minimum_time_step = 1.e-9 |units.Myr

star_type = {   'all': -1,
                'lm_ms': 0,
                'ms': 1, 
                'hg': 2,
                'rgb': 3, 
                'cheb': 4, 
                'eagb': 5, 
                'tpagb': 6, 
                'agb':[5,6], 
                'hems':7, 
                'hehg':8, 
                'hergb':9, 
                'heg':[8,9], 
                'hewd': 10, 
                'cowd': 11, 
                'onewd': 12, 
                'wd': [10,11,12],
                'ns': 13,
                'bh': 14, 
                'sn': 15,
                'unknown': 16, 
                'prems': 17, 
                'planet': 18,
                'bd': 19,
            }
                
bin_type = {    'all': -1,
                'unknown': 0,       
                'merger': 1, 
                'disintegrated': 2, 
                'dynamical_instability': 3, 
                'detached': 4,       
                'contact': 5,    
                'collision': 6,    
                'semisecular': 7,      
                'rlof': 8,   #only used for stopping conditions
                'olof': 9,   #only used for stopping conditions
                'stable_mass_transfer': 10, 
                'common_envelope': 11,     
                'common_envelope_energy_balance': 12,     
                'common_envelope_angular_momentum_balance': 13,
                'double_common_envelope': 14,
            }            

tr_type = {     'all': -1,
                'hierarchical': 0, 
                'dynamical_instability': 1, 
                'semisecular_regime': 2, 
                'error_flag_secular': 3, 
            }

lib_print_style = { 0: "No printing", #default
                    1: "Printing", 
                    2: "Readable format",
                    3: "Full particle"
                }

lib_parameter_style = { 0: "Full", #
                    1: "TRES standard; selected parameters", #default
                    2: "User specified selected parameters",
                }




def print_particle(particle):
        if particle.is_star:
            print(particle)                 
        else:
            print(particle) 
            print_particle(particle.child1)                
            print_particle(particle.child2)                



def create_snapshot_for_dict(triple, parameter_style):
    if parameter_style==0: #all available parameters
        snapshot_props = { 
            'system_ID': triple[0].number,
            'time_Myr': (triple[0].time).value_in(units.Myr),
            'incl': triple[0].relative_inclination,
            'dynamical_instabiity': int(triple[0].dynamical_instability),
            'kozai_type': triple[0].kozai_type,
            'error_flag_secular': triple[0].error_flag_secular,
            'CPU_time': triple[0].CPU_time,
            'delta_e_in': triple[0].delta_e_in,
            'max_delta_e_in': triple[0].max_delta_e_in,
          
            'outer_bin_type': bin_type[triple[0].bin_type], 
            'outer_mt_Stable': int(triple[0].is_mt_stable),
            'outer_semi_RSun': (triple[0].semimajor_axis).value_in(units.RSun),
            'outer_ecc': triple[0].eccentricity,
            'outer_aop': triple[0].argument_of_pericenter,
            'outer_loan': triple[0].longitude_of_ascending_node,
            #'TripleMDotMSunYr': (triple[0].mass_transfer_rate).value_in(units.MSun/units.yr),
    
            'inner_bin_type': bin_type[triple[0].child2.bin_type],
            'inner_mt_stable': int(triple[0].child2.is_mt_stable),
            'inner_semi_RSun': triple[0].child2.semimajor_axis.value_in(units.RSun),
            'inner_ecc': triple[0].child2.eccentricity,
            'inner_aop': triple[0].child2.argument_of_pericenter,
            'inner_loan': triple[0].child2.longitude_of_ascending_node,
    
            'star1_is_donor': int(triple[0].child2.child1.is_donor),
#            'star1_is_OLOF': int(triple[0].child2.child1.is_OLOF),  
            'star1_stellar_type': triple[0].child2.child1.stellar_type.value_in(units.stellar_type),
            'star1_mass_MSun': triple[0].child2.child1.mass.value_in(units.MSun),
            'star1_spin_Myr_inv': triple[0].child2.child1.spin_angular_frequency.value_in(1./units.Myr),
            'star1_rad_RSun': triple[0].child2.child1.radius.value_in(units.RSun),
            'star1_coremass_MSun': triple[0].child2.child1.core_mass.value_in(units.MSun),
            'star1_previous_mass_MSun': triple[0].child2.child1.previous_mass.value_in(units.MSun),
            'star1_luminosity_LSun': triple[0].child2.child1.luminosity.value_in(units.LSun),  
            'star1_temperature_K': triple[0].child2.child1.temperature.value_in(units.K),  
    
            'star2_is_donor': int(triple[0].child2.child2.is_donor),
#            'star2_is_OLOF': int(triple[0].child2.child2.is_OLOF),  
            'star2_stellar_type': triple[0].child2.child2.stellar_type.value_in(units.stellar_type),
            'star2_mass_MSun': triple[0].child2.child2.mass.value_in(units.MSun),
            'star2_spin_Myr_inv': triple[0].child2.child2.spin_angular_frequency.value_in(1./units.Myr),
            'star2_rad_RSun': triple[0].child2.child2.radius.value_in(units.RSun),
            'star2_coremass_MSun': triple[0].child2.child2.core_mass.value_in(units.MSun),
            'star2_previous_mass_MSun': triple[0].child2.child2.previous_mass.value_in(units.MSun),
            'star2_luminosity_LSun': triple[0].child2.child2.luminosity.value_in(units.LSun),  
            'star2_temperature_K': triple[0].child2.child2.temperature.value_in(units.K),  
    
            'star3_is_donor': int(triple[0].child1.is_donor),
#            'star3_is_OLOF': int(triple[0].child1.is_OLOF),  
            'star3_stellar_type': triple[0].child1.stellar_type.value_in(units.stellar_type),
            'star3_mass_MSun': triple[0].child1.mass.value_in(units.MSun),
            'star3_spin_Myr_inv': triple[0].child1.spin_angular_frequency.value_in(1./units.Myr),
            'star3_rad_RSun': triple[0].child1.radius.value_in(units.RSun),
            'star3_coremass_MSun': triple[0].child1.core_mass.value_in(units.MSun),
            'star3_previous_mass_MSun': triple[0].child1.previous_mass.value_in(units.MSun),
            'star3_luminosity_LSun': triple[0].child1.luminosity.value_in(units.LSun),  
            'star3_temperature_K': triple[0].child1.temperature.value_in(units.K),  
    
            }
            
    elif parameter_style==2: 
        #User selected parameters - select by commenting / uncommenting lines
        snapshot_props = {
            'system_ID': triple[0].number,
            'time_Myr': (triple[0].time).value_in(units.Myr),
            'incl': triple[0].relative_inclination,
            'dynamical_instabiity': int(triple[0].dynamical_instability),
            'kozai_type': triple[0].kozai_type,
            'error_flag_secular': triple[0].error_flag_secular,
            'CPU_time': triple[0].CPU_time,
#            'delta_e_in': triple[0].delta_e_in,
#            'max_delta_e_in': triple[0].max_delta_e_in,
          
            'outer_bin_type': bin_type[triple[0].bin_type], 
#            'outer_mt_Stable': int(triple[0].is_mt_stable),
            'outer_semi_RSun': (triple[0].semimajor_axis).value_in(units.RSun),
            'outer_ecc': triple[0].eccentricity,
            'outer_aop': triple[0].argument_of_pericenter,
            'outer_loan': triple[0].longitude_of_ascending_node,
            #'TripleMDotMSunYr': (triple[0].mass_transfer_rate).value_in(units.MSun/units.yr),
    
            'inner_bin_type': bin_type[triple[0].child2.bin_type],
#            'inner_mt_stable': int(triple[0].child2.is_mt_stable),
            'inner_semi_RSun': triple[0].child2.semimajor_axis.value_in(units.RSun),
            'inner_ecc': triple[0].child2.eccentricity,
            'inner_aop': triple[0].child2.argument_of_pericenter,
            'inner_loan': triple[0].child2.longitude_of_ascending_node,
    
            'star1_is_donor': int(triple[0].child2.child1.is_donor),
#            'star1_is_OLOF': int(triple[0].child2.child1.is_OLOF),  
            'star1_stellar_type': triple[0].child2.child1.stellar_type.value_in(units.stellar_type),
            'star1_mass_MSun': triple[0].child2.child1.mass.value_in(units.MSun),
            'star1_spin_Myr_inv': triple[0].child2.child1.spin_angular_frequency.value_in(1./units.Myr),
            'star1_rad_RSun': triple[0].child2.child1.radius.value_in(units.RSun),
            'star1_coremass_MSun': triple[0].child2.child1.core_mass.value_in(units.MSun),
#            'star1_previous_mass_MSun': triple[0].child2.child1.previous_mass.value_in(units.MSun),
#            'star1_luminosity_LSun': triple[0].child2.child1.luminosity.value_in(units.LSun),  
#            'star1_temperature_K': triple[0].child2.child1.temperature.value_in(units.K),  
    
            'star2_is_donor': int(triple[0].child2.child2.is_donor),
#            'star2_is_OLOF': int(triple[0].child2.child2.is_OLOF),  
            'star2_stellar_type': triple[0].child2.child2.stellar_type.value_in(units.stellar_type),
            'star2_mass_MSun': triple[0].child2.child2.mass.value_in(units.MSun),
            'star2_spin_Myr_inv': triple[0].child2.child2.spin_angular_frequency.value_in(1./units.Myr),
            'star2_rad_RSun': triple[0].child2.child2.radius.value_in(units.RSun),
            'star2_coremass_MSun': triple[0].child2.child2.core_mass.value_in(units.MSun),
#            'star2_previous_mass_MSun': triple[0].child2.child2.previous_mass.value_in(units.MSun),
#            'star2_luminosity_LSun': triple[0].child2.child2.luminosity.value_in(units.LSun),  
#            'star2_temperature_K': triple[0].child2.child2.temperature.value_in(units.K),  
    
            'star3_is_donor': int(triple[0].child1.is_donor),
#            'star3_is_OLOF': int(triple[0].child1.is_OLOF),  
            'star3_stellar_type': triple[0].child1.stellar_type.value_in(units.stellar_type),
            'star3_mass_MSun': triple[0].child1.mass.value_in(units.MSun),
            'star3_spin_Myr_inv': triple[0].child1.spin_angular_frequency.value_in(1./units.Myr),
            'star3_rad_RSun': triple[0].child1.radius.value_in(units.RSun),
            'star3_coremass_MSun': triple[0].child1.core_mass.value_in(units.MSun),
#            'star3_previous_mass_MSun': triple[0].child1.previous_mass.value_in(units.MSun),
#            'star3_luminosity_LSun': triple[0].child1.luminosity.value_in(units.LSun),  
#            'star3_temperature_K': triple[0].child1.temperature.value_in(units.K),  
    
            } 
    else: #lib_parameter_style[parameter_style]==1: #selected parameters        
        snapshot_props = {
            'system_ID': triple[0].number,
            'time_Myr': (triple[0].time).value_in(units.Myr),
            'incl': triple[0].relative_inclination,
            'dynamical_instabiity': int(triple[0].dynamical_instability),
            'kozai_type': triple[0].kozai_type,
            'error_flag_secular': triple[0].error_flag_secular,
            'CPU_time': triple[0].CPU_time,
#            'delta_e_in': triple[0].delta_e_in,
#            'max_delta_e_in': triple[0].max_delta_e_in,
          
            'outer_bin_type': bin_type[triple[0].bin_type], 
#            'outer_mt_Stable': int(triple[0].is_mt_stable),
            'outer_semi_RSun': (triple[0].semimajor_axis).value_in(units.RSun),
            'outer_ecc': triple[0].eccentricity,
            'outer_aop': triple[0].argument_of_pericenter,
            'outer_loan': triple[0].longitude_of_ascending_node,
            #'TripleMDotMSunYr': (triple[0].mass_transfer_rate).value_in(units.MSun/units.yr),
    
            'inner_bin_type': bin_type[triple[0].child2.bin_type],
#            'inner_mt_stable': int(triple[0].child2.is_mt_stable),
            'inner_semi_RSun': triple[0].child2.semimajor_axis.value_in(units.RSun),
            'inner_ecc': triple[0].child2.eccentricity,
            'inner_aop': triple[0].child2.argument_of_pericenter,
            'inner_loan': triple[0].child2.longitude_of_ascending_node,
    
            'star1_is_donor': int(triple[0].child2.child1.is_donor),
#            'star1_is_OLOF': int(triple[0].child2.child1.is_OLOF),  
            'star1_stellar_type': triple[0].child2.child1.stellar_type.value_in(units.stellar_type),
            'star1_mass_MSun': triple[0].child2.child1.mass.value_in(units.MSun),
            'star1_spin_Myr_inv': triple[0].child2.child1.spin_angular_frequency.value_in(1./units.Myr),
            'star1_rad_RSun': triple[0].child2.child1.radius.value_in(units.RSun),
            'star1_coremass_MSun': triple[0].child2.child1.core_mass.value_in(units.MSun),
#            'star1_previous_mass_MSun': triple[0].child2.child1.previous_mass.value_in(units.MSun),
#            'star1_luminosity_LSun': triple[0].child2.child1.luminosity.value_in(units.LSun),  
#            'star1_temperature_K': triple[0].child2.child1.temperature.value_in(units.K),  
    
            'star2_is_donor': int(triple[0].child2.child2.is_donor),
#            'star2_is_OLOF': int(triple[0].child2.child2.is_OLOF),  
            'star2_stellar_type': triple[0].child2.child2.stellar_type.value_in(units.stellar_type),
            'star2_mass_MSun': triple[0].child2.child2.mass.value_in(units.MSun),
            'star2_spin_Myr_inv': triple[0].child2.child2.spin_angular_frequency.value_in(1./units.Myr),
            'star2_rad_RSun': triple[0].child2.child2.radius.value_in(units.RSun),
            'star2_coremass_MSun': triple[0].child2.child2.core_mass.value_in(units.MSun),
#            'star2_previous_mass_MSun': triple[0].child2.child2.previous_mass.value_in(units.MSun),
#            'star2_luminosity_LSun': triple[0].child2.child2.luminosity.value_in(units.LSun),  
#            'star2_temperature_K': triple[0].child2.child2.temperature.value_in(units.K),  
    
            'star3_is_donor': int(triple[0].child1.is_donor),
#            'star3_is_OLOF': int(triple[0].child1.is_OLOF),  
            'star3_stellar_type': triple[0].child1.stellar_type.value_in(units.stellar_type),
            'star3_mass_MSun': triple[0].child1.mass.value_in(units.MSun),
            'star3_spin_Myr_inv': triple[0].child1.spin_angular_frequency.value_in(1./units.Myr),
            'star3_rad_RSun': triple[0].child1.radius.value_in(units.RSun),
            'star3_coremass_MSun': triple[0].child1.core_mass.value_in(units.MSun),
#            'star3_previous_mass_MSun': triple[0].child1.previous_mass.value_in(units.MSun),
#            'star3_luminosity_LSun': triple[0].child1.luminosity.value_in(units.LSun),  
#            'star3_temperature_K': triple[0].child1.temperature.value_in(units.K),  
    
            } 
          
                  
    return snapshot_props       






       
#for more info on mass transfer stability, see triple[0].is_mt_stable & triple[0].child2.is_mt_stable
def rdc(file_name_root, parameter_style, print_style, save_all_snapshots, print_init, line_number, inner_primary_star_type, inner_secondary_star_type, outer_star_type, inner_primary_star_type_string, inner_secondary_star_type_string, outer_star_type_string, 
inner_bin_type, outer_bin_type, inner_bin_type_string, outer_bin_type_string, triple_type, triple_type_string):


    file_name = file_name_root + ".hdf"            
    if file_name_root[-4:]==".hdf":
        file_name = file_name_root

    print('Reducing file: \t\t\t', file_name)
    print('Parameter style: \t\t', lib_parameter_style[parameter_style])
    print('Print style: \t\t\t', lib_print_style[print_style])
    print('Save_all_snapshots: \t\t', save_all_snapshots)
    print('')



    triple=read_set_from_file(file_name , "hdf5")
#    counter = list(enumerate(triple.history))[0][1].number 


    if print_init:
        for i, triple in enumerate(triple.history):
            if i == line_number:
                print('amuse TRES.py ', end = '' )
                print(' -M ', triple[0].child2.child1.mass.value_in(units.MSun), ' -m ', triple[0].child2.child2.mass.value_in(units.MSun), ' -l ', triple[0].child1.mass.value_in(units.MSun), end = '')
                print(' -A ', triple[0].child2.semimajor_axis.value_in(units.RSun), ' -a ', triple[0].semimajor_axis.value_in(units.RSun), end = '')
                print(' -E ', triple[0].child2.eccentricity, ' -e ', triple[0].eccentricity, end = '')
                print(' -G ', triple[0].child2.argument_of_pericenter, ' -g ', triple[0].argument_of_pericenter, end = '')
                print(' -I ', triple[0].relative_inclination, end = '\t')
                
                if triple[0].time > minimum_time_step:
                    print('Warning: these parameters do not represent a system at birth (ZAMS).')

        return


#    print(inner_primary_star_type, inner_secondary_star_type, outer_star_type)
#    print(star_type[inner_primary_star_type_string],star_type[inner_secondary_star_type_string],star_type[outer_star_type_string] )
#    print(inner_bin_type, outer_bin_type, triple_type)
#    print(bin_type[inner_bin_type_string], bin_type[outer_bin_type_string], tr_type[triple_type_string])

    if inner_primary_star_type != -1 and star_type[inner_primary_star_type_string] != -1 and inner_primary_star_type != star_type[inner_primary_star_type_string] :
        print('error: two different inner primary star types requested')
        return
    if inner_secondary_star_type != -1 and star_type[inner_secondary_star_type_string] != -1 and inner_secondary_star_type != star_type[inner_secondary_star_type_string] :
        print('error: two different inner secondary star types requested')
        return
    if outer_star_type != -1 and star_type[outer_star_type_string] != -1 and outer_star_type != star_type[outer_star_type_string] :
        print('error: two different outer star types requested')
        return
    if inner_bin_type > -1 and bin_type[inner_bin_type_string] > -1 and inner_bin_type != bin_type[inner_bin_type_string] :
        print('error: two different inner binary types requested')
        return
    if outer_bin_type > -1 and bin_type[outer_bin_type_string] > -1 and outer_bin_type != bin_type[outer_bin_type_string]:
        print('error: two different outer binary types requested')
        return
    if triple_type > -1 and tr_type[triple_type_string] > -1 and triple_type != tr_type[triple_type_string] :
        print('error: two different triple types requested')
        return        
        
    if star_type[inner_primary_star_type_string] != -1:
        inner_primary_star_type = star_type[inner_primary_star_type_string]
    if star_type[inner_secondary_star_type_string] != -1:
        inner_secondary_star_type = star_type[inner_secondary_star_type_string]
    if star_type[outer_star_type_string] != -1:
        outer_star_type = star_type[outer_star_type_string]
    inner_bin_type = max(inner_bin_type, bin_type[inner_bin_type_string])    
    outer_bin_type = max(outer_bin_type, bin_type[outer_bin_type_string]) 
    triple_type = max(triple_type, tr_type[triple_type_string]) 
    

    snapshot_props = []
    previous_triple_number = -1
    previous_correct_system = False
    
    for i, triple in enumerate(triple.history):

        if print_style == 3:
            print_particle(triple[0])
            sys.exit(0)

        triple_number = triple[0].number
        previous_snapshot_props = snapshot_props
        snapshot_props = create_snapshot_for_dict(triple, parameter_style)
        
        if triple_type != tr_type['all']:
            triple_type_snapshot = tr_type['hierarchical']
            if triple[0].dynamical_instability:
                triple_type_snapshot = tr_type['dynamical_instability']
#                if triple[0].semisecular_regime:
#                    triple_type_snapshot = tr_type['semisecular_regime']
            if triple[0].error_flag_secular:
                triple_type_snapshot = tr_type['error_flag_secular']
    
        if (inner_primary_star_type == star_type['all'] or triple[0].child2.child1.stellar_type.value_in(units.stellar_type) in np.array(inner_primary_star_type)) and \
        (inner_secondary_star_type == star_type['all'] or triple[0].child2.child2.stellar_type.value_in(units.stellar_type) in np.array(inner_secondary_star_type)) and \
        (outer_star_type == star_type['all'] or triple[0].child1.stellar_type.value_in(units.stellar_type) in np.array(outer_star_type)) and \
        (inner_bin_type == bin_type['all'] or inner_bin_type == bin_type[triple[0].child2.bin_type]) and \
        (outer_bin_type == bin_type['all'] or outer_bin_type == bin_type[triple[0].bin_type]) and \
        (triple_type == tr_type['all'] or triple_type == triple_type_snapshot):
            correct_system = True

                       
        #which snapshots to save
        if save_all_snapshots and correct_system: #all snapshots    
            SaveDictSet = []           
            SaveDictSet.append(snapshot_props)
            df = pd.DataFrame(SaveDictSet)                   
            df.to_csv(file_name_root+'.csv', mode='a', header=not os.path.exists(file_name_root+'.csv'))       
            if print_style == 1:
                if i==0:
                    print(df.to_string())                
                print(df.to_string(header=False))
        elif triple_number>previous_triple_number: #new system
            #save last line previous system if necessary        
            if previous_correct_system:
                SaveDictSet.append(previous_snapshot_props)
                df = pd.DataFrame(SaveDictSet)                   
                df.to_csv(file_name_root+'.csv', mode='a', header=not os.path.exists(file_name_root+'.csv'))    
                if print_style == 1:
                    if triple_number==0:
                        print(df.to_string())                
                    print(df.to_string(header=False))
    
            #first line of current system        
            #saving for now
            SaveDictSet = []           
            SaveDictSet.append(snapshot_props)
            previous_triple_number = triple_number  
            previous_correct_system = correct_system


    #last line for last system    
    if correct_system and save_all_snapshots==False:
        SaveDictSet.append(snapshot_props)
        df = pd.DataFrame(SaveDictSet)                   
        df.to_csv(file_name_root+'.csv', mode='a', header=not os.path.exists(file_name_root+'.csv'))    
        if print_style == 1:
                print(df.to_string(header=False))

          
    

def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", dest="file_name_root", default = "TRES",
                      help="file name [%default]")                      
    parser.add_option("-P", dest="parameter_style", type="int", default = 1,
                      help="parameter style - which parameters should be recorded[%default]") 
    parser.add_option("-S", dest="print_style", type="int", default = 0,
                      help="print style [%default]") 
    parser.add_option("-F", "--save_all_snapshots", dest="save_all_snapshots", action="store_true", default = False, 
                      help="save every snapshot for specified triple [%default]")
    parser.add_option("--print_init", dest="print_init", action="store_true", default = False, 
                      help="print initial conditions for re running [%default]")
    parser.add_option("-l", dest="line_number", type="int", default = 0,
                      help="line number for printing initial conditions [%default]") #will only do something when print_init = True

    #returns first instance where desired star_type, bin_type & triple_type is reached 
    parser.add_option("--st1", dest="inner_primary_star_type", type="int", default = -1,
                      help="desired stellar type of inner binary primary star (int) [%default]") 
    parser.add_option("--st2", dest="inner_secondary_star_type", type="int", default = -1,
                      help="desired stellar type of inner binary secondary star (int) [%default]") 
    parser.add_option("--st3", dest="outer_star_type", type="int", default = -1,
                      help="desired stellar type of tertiary star (int) [%default]") 
    parser.add_option("--st1str", dest="inner_primary_star_type_string", default = "all",
                      help="desired stellar type of inner binary primary star (int) [%default]") 
    parser.add_option("--st2str", dest="inner_secondary_star_type_string", default = "all",
                      help="desired stellar type of inner binary secondary star (int) [%default]") 
    parser.add_option("--st3str", dest="outer_star_type_string", default = "all",
                      help="desired stellar type of tertiary star [%default]") 
    parser.add_option("--btin", dest="inner_bin_type", type="int", default = -1,
                      help="desired binary type of inner binary (int) [%default]") 
    parser.add_option("--btout", dest="outer_bin_type", type="int", default = -1,
                      help="desired binary type of inner binary (int) [%default]") 
    parser.add_option("--btinstr", dest="inner_bin_type_string", default = "all",
                      help="desired binary type of inner binary (string) [%default]")                      
    parser.add_option("--btoutstr", dest="outer_bin_type_string", default = "all",
                      help="desired binary type of outer binary (string) [%default]")                      
    parser.add_option("--trt", dest="triple_type", type="int", default = -1,
                      help="desired triple type (int) [%default]") 
    parser.add_option("--trtstr", dest="triple_type_string", default = "all",
                      help="desired triple type [%default]")                      

                      
    options, args = parser.parse_args()
    return options.__dict__


if __name__ == '__main__':
    options = parse_arguments()

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")



    print(' ')
    rdc(**options)  
    
