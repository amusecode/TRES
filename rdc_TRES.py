from amuse.units import units, constants
from amuse.io import write_set_to_file
from amuse.io import read_set_from_file
from amuse.support.console import set_printing_strategy
import numpy as np
import io
import sys
    
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

lib_print_style = { 0: "Readable format",
                    1: "Full",
                    2: "TRES standard; selected parameters", 
                    3: "TRES standard; selected parameters - csv style", 
                }#default

def print_particle(particle):
        if particle.is_star:
            print(particle)                 
        else:
            print(particle) 
            print_particle(particle.child1)                
            print_particle(particle.child2)                


def print_to_string(*args, **kwargs):
    output = io.StringIO()
    print(*args, file=output, **kwargs, sep=' ')
    contents = output.getvalue()
    output.close()
    return contents

def print_to_string_csv(*args, **kwargs):
    output = io.StringIO()
    print(*args, file=output, **kwargs, sep=',')
    contents = output.getvalue()
    output.close()
    return contents

       
#for more info on mass transfer stability, see triple[0].is_mt_stable & triple[0].child2.is_mt_stable
def rdc(file_name_root, print_style, print_full, print_init, line_number, inner_primary_star_type, inner_secondary_star_type, outer_star_type, inner_primary_star_type_string, inner_secondary_star_type_string, outer_star_type_string, 
inner_bin_type, outer_bin_type, inner_bin_type_string, outer_bin_type_string, triple_type, triple_type_string):

    file_name = file_name_root + ".hdf"            
    if file_name_root[-4:]==".hdf":
        file_name = file_name_root

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


    correct_system = False
    correct_system_previous = False
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
    
#    print(inner_primary_star_type, inner_secondary_star_type, outer_star_type)
#    print(star_type[inner_primary_star_type_string],star_type[inner_secondary_star_type_string],star_type[outer_star_type_string] )    
#    print(inner_bin_type, outer_bin_type, triple_type)
#    print(bin_type[inner_bin_type_string], bin_type[outer_bin_type_string], tr_type[triple_type_string])

    print('#', lib_print_style[print_style])
    triple_string = ''
    snapshot_string = '' 
    if print_style == 2:
        snapshot_string = snapshot_string + '# '
    if print_style == 2 or print_style == 3:
        snapshot_string = snapshot_string + 'number, time, incl, dyn_inst, kozai_type, error_flag, CPU_time,'
        snapshot_string = snapshot_string + ' bin_type2, a2, e2, p2, l2, bin_type1, a1, e1, p1, l1,'
        snapshot_string = snapshot_string + ' m1_donor, m1_stellar_type, m1_mass, m1_spin, m1_radius, m1_core,'
        snapshot_string = snapshot_string + ' m2_donor, m2_stellar_type, m2_mass, m2_spin, m2_radius, m2_core,'
        snapshot_string = snapshot_string + ' m3_donor, m3_stellar_type, m3_mass, m3_spin, m3_radius, m3_core,'
        snapshot_string = snapshot_string + '\n' 
    triple_number = 0
    previous_triple_number = -1 
    
    for i, triple in enumerate(triple.history):
#        print(triple[0].number, triple[0].time)
#        if triple[0].number == counter:
#            counter += 1    
#            print('\n\n')

        #which snapshots to save
        if print_full and correct_system: #all snapshots    
            triple_string = triple_string + snapshot_string   
        else: #first & last line 
#            print(i, triple[0].number,triple_number, previous_triple_number)
            if triple[0].number>triple_number: # last line & correct_system
                 triple_string = triple_string + snapshot_string
                 triple_number = triple[0].number
                 if correct_system:
                     print(triple_string[:-2]) #prevents empty line
                 triple_string = ''
                 correct_system = False
                 correct_system_previous = False 
            elif triple_number > previous_triple_number: #first line
                 triple_string = triple_string + snapshot_string
                 previous_triple_number = triple_number  
                 if correct_system == True:
                    correct_system_previous = True
            elif correct_system==True and correct_system_previous==False:
                 triple_string = triple_string + snapshot_string
                 previous_triple_number = triple_number  
                 correct_system_previous = True                                                
        snapshot_string = '' 
        
        if triple_type != tr_type['all']:
            triple_type_snapshot = tr_type['hierarchical']
            if triple[0].dynamical_instability:
                triple_type_snapshot = tr_type['dynamical_instability']
#                if triple[0].semisecular_regime:
#                    triple_type_snapshot = tr_type['semisecular_regime']
            if triple[0].error_flag_secular:
                triple_type_snapshot = tr_type['error_flag_secular']
        
        
        if print_style == 0:
#            print(' ')
#            print(i, triple[0].number, triple[0].time, triple[0].relative_inclination, triple[0].dynamical_instability, triple[0].kozai_type, triple[0].error_flag_secular, triple[0].CPU_time)
#                  
#            print( ' bs: ', triple[0].child2.bin_type, triple[0].child2.semimajor_axis, triple[0].child2.eccentricity, triple[0].child2.argument_of_pericenter, triple[0].child2.longitude_of_ascending_node,)# triple[0].child2.mass_transfer_rate,
#            print( '|', triple[0].bin_type, triple[0].semimajor_axis, triple[0].eccentricity, triple[0].argument_of_pericenter, triple[0].longitude_of_ascending_node)#, triple[0].mass_transfer_rate
#            print( ' st: ',  triple[0].child2.child1.is_donor, triple[0].child2.child1.stellar_type, triple[0].child2.child1.mass,  triple[0].child2.child1.spin_angular_frequency, triple[0].child2.child1.radius, triple[0].child2.child1.core_mass,)
#            print( '|', triple[0].child2.child2.is_donor,  triple[0].child2.child2.stellar_type, triple[0].child2.child2.mass, triple[0].child2.child2.spin_angular_frequency, triple[0].child2.child2.radius,triple[0].child2.child2.core_mass, )
#            print( '|', triple[0].child1.is_donor, triple[0].child1.stellar_type, triple[0].child1.mass, triple[0].child1.spin_angular_frequency, triple[0].child1.radius, triple[0].child1.core_mass)

            if i>0:
                snapshot_string = snapshot_string + print_to_string(' ')
            snapshot_string = snapshot_string + print_to_string(i, triple[0].number, triple[0].time, triple[0].relative_inclination, triple[0].dynamical_instability, triple[0].kozai_type, triple[0].error_flag_secular, triple[0].CPU_time)
            snapshot_string = snapshot_string + print_to_string(' bs: ', triple[0].child2.bin_type, triple[0].child2.semimajor_axis, triple[0].child2.eccentricity, triple[0].child2.argument_of_pericenter, triple[0].child2.longitude_of_ascending_node)
            snapshot_string = snapshot_string + print_to_string('|', triple[0].bin_type, triple[0].semimajor_axis, triple[0].eccentricity, triple[0].argument_of_pericenter, triple[0].longitude_of_ascending_node)
            snapshot_string = snapshot_string + print_to_string(' st: ',  triple[0].child2.child1.is_donor, triple[0].child2.child1.stellar_type, triple[0].child2.child1.mass,  triple[0].child2.child1.spin_angular_frequency, triple[0].child2.child1.radius, triple[0].child2.child1.core_mass)
            snapshot_string = snapshot_string + print_to_string('|', triple[0].child2.child2.is_donor,  triple[0].child2.child2.stellar_type, triple[0].child2.child2.mass, triple[0].child2.child2.spin_angular_frequency, triple[0].child2.child2.radius,triple[0].child2.child2.core_mass)
            snapshot_string = snapshot_string + print_to_string('|', triple[0].child1.is_donor, triple[0].child1.stellar_type, triple[0].child1.mass, triple[0].child1.spin_angular_frequency, triple[0].child1.radius, triple[0].child1.core_mass)

        elif print_style == 1:
            print_particle(triple[0])
        
#            print('triple particle', triple[0])
#            print('child1', triple[0].child1)
#            print('child2',triple[0].child2)
#            print('child2.child1',triple[0].child2.child1)
#            print('child2.child2',triple[0].child2.child2)
            sys.exit(0)
        else:
#            print(triple[0].number, triple[0].time.value_in(units.Myr), triple[0].relative_inclination, int(triple[0].dynamical_instability), int(triple[0].kozai_type), int(triple[0].error_flag_secular), triple[0].CPU_time, end = '\t')
#            print(bin_type[triple[0].child2.bin_type], triple[0].child2.semimajor_axis.value_in(units.RSun), triple[0].child2.eccentricity, triple[0].child2.argument_of_pericenter, triple[0].child2.longitude_of_ascending_node, end = '\t')
#            print(bin_type[triple[0].bin_type], triple[0].semimajor_axis.value_in(units.RSun), triple[0].eccentricity, triple[0].argument_of_pericenter, triple[0].longitude_of_ascending_node, end = '\t')
#            print(int(triple[0].child2.child1.is_donor), triple[0].child2.child1.stellar_type.value_in(units.stellar_type), triple[0].child2.child1.mass.value_in(units.MSun),  triple[0].child2.child1.spin_angular_frequency.value_in(1./units.Myr), triple[0].child2.child1.radius.value_in(units.RSun), triple[0].child2.child1.core_mass.value_in(units.MSun), end = '\t')
#            print(int(triple[0].child2.child2.is_donor),  triple[0].child2.child2.stellar_type.value_in(units.stellar_type), triple[0].child2.child2.mass.value_in(units.MSun), triple[0].child2.child2.spin_angular_frequency.value_in(1./units.Myr), triple[0].child2.child2.radius.value_in(units.RSun),triple[0].child2.child2.core_mass.value_in(units.MSun), end = '\t' )
#            print(int(triple[0].child1.is_donor), triple[0].child1.stellar_type.value_in(units.stellar_type), triple[0].child1.mass.value_in(units.MSun), triple[0].child1.spin_angular_frequency.value_in(1./units.Myr), triple[0].child1.radius.value_in(units.RSun), triple[0].child1.core_mass.value_in(units.MSun))
  
  
            if print_style == 3:
                snapshot_string = snapshot_string + print_to_string_csv(triple[0].number, triple[0].time.value_in(units.Myr), triple[0].relative_inclination, int(triple[0].dynamical_instability), int(triple[0].kozai_type), int(triple[0].error_flag_secular), triple[0].CPU_time, end = ',')
                snapshot_string = snapshot_string + print_to_string_csv(bin_type[triple[0].child2.bin_type], triple[0].child2.semimajor_axis.value_in(units.RSun), triple[0].child2.eccentricity, triple[0].child2.argument_of_pericenter, triple[0].child2.longitude_of_ascending_node, end = ',')
                snapshot_string = snapshot_string + print_to_string_csv(bin_type[triple[0].bin_type], triple[0].semimajor_axis.value_in(units.RSun), triple[0].eccentricity, triple[0].argument_of_pericenter, triple[0].longitude_of_ascending_node, end = ',')
                snapshot_string = snapshot_string + print_to_string_csv(int(triple[0].child2.child1.is_donor), triple[0].child2.child1.stellar_type.value_in(units.stellar_type), triple[0].child2.child1.mass.value_in(units.MSun),  triple[0].child2.child1.spin_angular_frequency.value_in(1./units.Myr), triple[0].child2.child1.radius.value_in(units.RSun), triple[0].child2.child1.core_mass.value_in(units.MSun), end = ',')
                snapshot_string = snapshot_string + print_to_string_csv(int(triple[0].child2.child2.is_donor),  triple[0].child2.child2.stellar_type.value_in(units.stellar_type), triple[0].child2.child2.mass.value_in(units.MSun), triple[0].child2.child2.spin_angular_frequency.value_in(1./units.Myr), triple[0].child2.child2.radius.value_in(units.RSun),triple[0].child2.child2.core_mass.value_in(units.MSun), end = ',')
                snapshot_string = snapshot_string + print_to_string_csv(int(triple[0].child1.is_donor), triple[0].child1.stellar_type.value_in(units.stellar_type), triple[0].child1.mass.value_in(units.MSun), triple[0].child1.spin_angular_frequency.value_in(1./units.Myr), triple[0].child1.radius.value_in(units.RSun), triple[0].child1.core_mass.value_in(units.MSun))

            else:
                snapshot_string = snapshot_string + print_to_string(triple[0].number, triple[0].time.value_in(units.Myr), triple[0].relative_inclination, int(triple[0].dynamical_instability), int(triple[0].kozai_type), int(triple[0].error_flag_secular), triple[0].CPU_time, end = '\t')
                snapshot_string = snapshot_string + print_to_string(bin_type[triple[0].child2.bin_type], triple[0].child2.semimajor_axis.value_in(units.RSun), triple[0].child2.eccentricity, triple[0].child2.argument_of_pericenter, triple[0].child2.longitude_of_ascending_node, end = '\t')
                snapshot_string = snapshot_string + print_to_string(bin_type[triple[0].bin_type], triple[0].semimajor_axis.value_in(units.RSun), triple[0].eccentricity, triple[0].argument_of_pericenter, triple[0].longitude_of_ascending_node, end = '\t')
                snapshot_string = snapshot_string + print_to_string(int(triple[0].child2.child1.is_donor), triple[0].child2.child1.stellar_type.value_in(units.stellar_type), triple[0].child2.child1.mass.value_in(units.MSun),  triple[0].child2.child1.spin_angular_frequency.value_in(1./units.Myr), triple[0].child2.child1.radius.value_in(units.RSun), triple[0].child2.child1.core_mass.value_in(units.MSun), end = '\t')
                snapshot_string = snapshot_string + print_to_string(int(triple[0].child2.child2.is_donor),  triple[0].child2.child2.stellar_type.value_in(units.stellar_type), triple[0].child2.child2.mass.value_in(units.MSun), triple[0].child2.child2.spin_angular_frequency.value_in(1./units.Myr), triple[0].child2.child2.radius.value_in(units.RSun),triple[0].child2.child2.core_mass.value_in(units.MSun), end = '\t')
                snapshot_string = snapshot_string + print_to_string(int(triple[0].child1.is_donor), triple[0].child1.stellar_type.value_in(units.stellar_type), triple[0].child1.mass.value_in(units.MSun), triple[0].child1.spin_angular_frequency.value_in(1./units.Myr), triple[0].child1.radius.value_in(units.RSun), triple[0].child1.core_mass.value_in(units.MSun))
#            snapshot_string = snapshot_string + '\n'

        if (inner_primary_star_type == star_type['all'] or triple[0].child2.child1.stellar_type.value_in(units.stellar_type) in np.array(inner_primary_star_type)) and \
        (inner_secondary_star_type == star_type['all'] or triple[0].child2.child2.stellar_type.value_in(units.stellar_type) in np.array(inner_secondary_star_type)) and \
        (outer_star_type == star_type['all'] or triple[0].child1.stellar_type.value_in(units.stellar_type) in np.array(outer_star_type)) and \
        (inner_bin_type == bin_type['all'] or inner_bin_type == bin_type[triple[0].child2.bin_type]) and \
        (outer_bin_type == bin_type['all'] or outer_bin_type == bin_type[triple[0].bin_type]) and \
        (triple_type == tr_type['all'] or triple_type == triple_type_snapshot):
            correct_system = True

        if i==0 and print_full==False:
            triple_string = triple_string + snapshot_string
            #in case triple[0].number in file doesn't start at 0
            triple_number = triple[0].number
            previous_triple_number = triple_number
            if correct_system == True: 
                correct_system_previous = True                         
    
    triple_string = triple_string + snapshot_string
    if correct_system:
        print(triple_string, end='')


def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", dest="file_name_root", default = "TRES",
                      help="file name [%default]")                      
    parser.add_option("-S", dest="print_style", type="int", default = 2,
                      help="print style [%default]") 
    parser.add_option("-F", dest="print_full", action="store_true", default = False, 
                      help="print every snapshot for specified triple [%default]")
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
    
