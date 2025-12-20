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
#keeping because of error flag
tr_type = {     'all': -1,
                'hierarchical': 0, 
                'dynamical_instability': 1, 
                'semisecular_regime': 2, 
                'error_flag_secular': 3, 
            }

lib_print_style = { 0: "Readable format",
                    1: "Full dataset",
                    2: "TRES standard; selected parameters", #default 
                }

def print_particle(particle):
        if particle.is_star:
            print(particle)                 
        else:
            print(particle) 
            print_particle(particle.child1)                
            print_particle(particle.child2)                

#def print_to_string_old(*args, **kwargs):
#    output = io.StringIO()
#    print(*args, file=output, **kwargs, sep=' ')
#    contents = output.getvalue()
#    output.close()
#    return contents

def print_to_string(array, sep):
    
    string = ''
    for i,k in enumerate(array):
        string += str(k) + sep
    
    return string

def print_header(start=''):

    string = start + ' '
    string += 'line number, binary number, time, -, dyn_inst, kozai_type, error_flag, CPU_time,'
    string += ' bin_type2, a2, e2, p2, l2,'
    string += ' m1_donor, m1_stellar_type, m1_mass, m1_spin, m1_radius, m1_core,'
    string += ' m2_donor, m2_stellar_type, m2_mass, m2_spin, m2_radius, m2_core,'
    string += '\n'
    
    string += '#units: time in Myr, masses/coremasses in MSun, separations & radii in RSun, spin in 1/Myr \n' 
    
    return string

       
#for more info on mass transfer stability, see triple[0].is_mt_stable & triple[0].is_mt_stable
def rdc(file_name_root, output_file, print_style, save_every_snapshot, print_init, line_number, primary_star_type, secondary_star_type, primary_star_type_string, secondary_star_type_string, 
bin_bin_type, bin_type_string,  triple_type, triple_type_string):

    if file_name_root[-4:]==".hdf":
        file_name = file_name_root
    else:
        file_name = file_name_root + ".hdf"            
    triple=read_set_from_file(file_name, "hdf5")

    if output_file in ['sys.stdout', 'std', 'screen', 'terminal']:  
        output_stream = sys.stdout
    else:
        output_stream = open(output_file, 'a')
        
    sep = ' '
    if output_file[-4:] == '.csv':
        sep = ', ' 
      
    if print_init:
        for i, triple in enumerate(triple.history):
            if i == line_number:
                print('amuse TRES.py ', end = '', file=output_stream)
                print(' -M ', triple[0].child1.mass.value_in(units.MSun), end = '', file=output_stream)
                print(' -m ', triple[0].child2.mass.value_in(units.MSun), end = '', file=output_stream)
                print(' -a ', triple[0].semimajor_axis.value_in(units.RSun), end = '', file=output_stream)
                print(' -e ', triple[0].eccentricity, end = '', file=output_stream)
                print(' -g ', triple[0].argument_of_pericenter, end = '', file=output_stream)
                print('', file=output_stream)
                
                if triple[0].time > minimum_time_step:
                    print('Warning: these parameters do not represent a system at birth (ZAMS).')
        return


    correct_system = False
    correct_system_previous = False

    if primary_star_type != -1 and star_type[primary_star_type_string] != -1 and primary_star_type != star_type[primary_star_type_string] :
        print('error: two different inner primary star types requested')
        return
    if secondary_star_type != -1 and star_type[secondary_star_type_string] != -1 and secondary_star_type != star_type[secondary_star_type_string] :
        print('error: two different inner secondary star types requested')
        return
    if bin_bin_type > -1 and bin_type[bin_type_string] > -1 and bin_bin_type != bin_type[bin_type_string] :
        print('error: two different inner binary types requested')
        return
    if triple_type > -1 and tr_type[triple_type_string] > -1 and triple_type != tr_type[triple_type_string] :
        print('error: two different triple types requested')
        return        
        
    if star_type[primary_star_type_string] != -1:
        primary_star_type = star_type[primary_star_type_string]
    if star_type[secondary_star_type_string] != -1:
        secondary_star_type = star_type[secondary_star_type_string]
    bin_bin_type = max(bin_bin_type, bin_type[bin_type_string])    
    triple_type = max(triple_type, tr_type[triple_type_string]) 

    print('# Reducing file: \t\t', file_name)
    print('# Print style: \t\t\t', lib_print_style[print_style])
    print('# Save every snapshot: \t\t', save_every_snapshot)

    triple_string = ''    
    if print_style == 1:
        snapshot_string = ''
    elif output_file[-4:]=='.csv':
        snapshot_string = print_header()
    else:
        snapshot_string = print_header('#')
    
    triple_number = 0
    previous_triple_number = -1 
    
    for i, triple in enumerate(triple.history):

        #which snapshots to save
        if save_every_snapshot and correct_system: #all snapshots    
            triple_string = triple_string + snapshot_string   
        else: #first & last line 
            if triple[0].number>triple_number: # last line & correct_system
                 triple_string = triple_string + snapshot_string
                 triple_number = triple[0].number
                 if correct_system:
                     print(triple_string[:-2], file=output_stream) #prevents empty line
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

        #move to dict keys when dicts are fixed in amuse
        keys_intro = [triple[0].number, 
            triple[0].time.value_in(units.Myr), 
            0, 
            triple[0].dynamical_instability, 
            triple[0].kozai_type, 
            triple[0].error_flag_secular, 
            triple[0].CPU_time]
                        
        keys_bin = [triple[0].bin_type, 
            triple[0].semimajor_axis.value_in(units.RSun), 
            triple[0].eccentricity, 
            triple[0].argument_of_pericenter, 
            triple[0].longitude_of_ascending_node]
            
        keys_prim = [triple[0].child1.is_donor, 
            triple[0].child1.stellar_type.value_in(units.stellar_type), 
            triple[0].child1.mass.value_in(units.MSun),  
            triple[0].child1.spin_angular_frequency.value_in(1./units.Myr), 
            triple[0].child1.radius.value_in(units.RSun), 
            triple[0].child1.core_mass.value_in(units.MSun)]
        
        keys_sec = [triple[0].child2.is_donor,  
            triple[0].child2.stellar_type.value_in(units.stellar_type), 
            triple[0].child2.mass.value_in(units.MSun), 
            triple[0].child2.spin_angular_frequency.value_in(1./units.Myr),     
            triple[0].child2.radius.value_in(units.RSun),
            triple[0].child2.core_mass.value_in(units.MSun)]
        
        if output_file[-4:]=='.txt':
            index_keys_intro_int = [0,3,4,5]
            for j in range(len(index_keys_intro_int)):
                keys_intro[index_keys_intro_int[j]] = int(keys_intro[index_keys_intro_int[j]])
            keys_bin[0] = bin_type[keys_bin[0]]
            keys_prim[0] = int(keys_prim[0])
            keys_prim[1] = int(keys_prim[1])
            keys_sec[0] = int(keys_sec[0])
            keys_sec[1] = int(keys_sec[1])


        if print_style == 0:
            snapshot_string += '\n'
            snapshot_string += 'tr: \t\t' + str(i) + sep + print_to_string(keys_intro, sep) + '\n'
            snapshot_string += 'bin in: \t' + print_to_string(keys_bin, sep) + '\n' 
            snapshot_string += 'st prim: \t' + print_to_string(keys_prim, sep) + '\n' 
            snapshot_string += '   sec: \t' + print_to_string(keys_sec, sep) + '\n' 
            
        elif print_style == 1:
            print(triple[0], file=output_stream)
            sys.exit(0)
        else: 
            snapshot_string += str(i) + sep + print_to_string(keys_intro, sep)
            snapshot_string += print_to_string(keys_bin, sep)  
            snapshot_string += print_to_string(keys_prim, sep)  
            snapshot_string += print_to_string(keys_sec, sep)  
            snapshot_string += '\n'


        if (primary_star_type == star_type['all'] or triple[0].child1.stellar_type.value_in(units.stellar_type) in np.array(primary_star_type)) and \
        (secondary_star_type == star_type['all'] or triple[0].child2.stellar_type.value_in(units.stellar_type) in np.array(secondary_star_type)) and \
        (bin_bin_type == bin_type['all'] or bin_bin_type == bin_type[triple[0].bin_type]) and \
        (triple_type == tr_type['all'] or triple_type == triple_type_snapshot):
            correct_system = True

        if i==0 and save_every_snapshot==False:
            triple_string += snapshot_string
            #in case triple[0].number in file doesn't start at 0
            triple_number = triple[0].number
            previous_triple_number = triple_number
            if correct_system == True: 
                correct_system_previous = True                         
    
    triple_string += snapshot_string
    if correct_system:
        print(triple_string, end='', file=output_stream)


def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", dest="file_name_root", default = "BIN",
                      help="file name [%default]")                      
    parser.add_option("-F", dest="output_file", type ="string", default = 'screen',
                      help="output file[%default]")
    parser.add_option("-S", dest="print_style", type="int", default = 2,
                      help="print style [%default]") 
    parser.add_option("--save_every_snapshot", dest="save_every_snapshot", action="store_true", default = False, 
                      help="save every snapshot for specified triple [%default]")
    parser.add_option("--print_init", dest="print_init", action="store_true", default = False, 
                      help="print initial conditions for re running [%default]")
    parser.add_option("-l", dest="line_number", type="int", default = 0,
                      help="line number for printing initial conditions [%default]") #will only do something when print_init = True

    #returns first instance where desired star_type, bin_type & triple_type is reached 
    parser.add_option("--st1", dest="primary_star_type", type="int", default = -1,
                      help="desired stellar type of inner binary primary star (int) [%default]") 
    parser.add_option("--st2", dest="secondary_star_type", type="int", default = -1,
                      help="desired stellar type of inner binary secondary star (int) [%default]") 
    parser.add_option("--st1str", dest="primary_star_type_string", default = "all",
                      help="desired stellar type of inner binary primary star (int) [%default]") 
    parser.add_option("--st2str", dest="secondary_star_type_string", default = "all",
                      help="desired stellar type of inner binary secondary star (int) [%default]") 
    parser.add_option("--btin", dest="bin_bin_type", type="int", default = -1,
                      help="desired binary type of inner binary (int) [%default]") 
    parser.add_option("--btinstr", dest="bin_type_string", default = "all",
                      help="desired binary type of inner binary (string) [%default]")                      
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

    rdc(**options)  
    
