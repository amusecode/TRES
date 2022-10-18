# reducing TRES output for production work.
# Assuming output file names are called TRES_0.hdf, TRES_1.hdf, TRES_2.hdf etc. 
# More specifically TRES_[i].pdf, TRES_[i+1].hdf, ..., TRES_[j].hdf
# -n sets value of i
# -N sets value of j

from optparse import OptionParser
import rdc_TRES as rdc

def rdc_TPS(startnr, nr_of_files):
    file_name_root = "TRES"
    print_style = 0
    print_full = True
    
    #default settings, adjust to your liking
    print_init = False
    line_number = 0 
    inner_primary_star_type =-1
    inner_secondary_star_type =-1
    outer_star_type =-1
    inner_primary_star_type_string="all"
    inner_secondary_star_type_string="all"
    outer_star_type_string="all"
    inner_bin_type=-1
    outer_bin_type=-1
    inner_bin_type_string="all"
    outer_bin_type_string="all"
    triple_type=-1
    triple_type_string="all"
    
    for i in range(nr_of_files):
        print(file_name_root+'_'+str(i+startnr))
        rdc.rdc(file_name_root+'_'+str(i+startnr), print_style, print_full, print_init, line_number, inner_primary_star_type, inner_secondary_star_type, outer_star_type, inner_primary_star_type_string, inner_secondary_star_type_string, outer_star_type_string, inner_bin_type, outer_bin_type, inner_bin_type_string, outer_bin_type_string, triple_type, triple_type_string)




def parse_arguments():
    parser = OptionParser()
    parser.add_option("-n", "--startnr", dest="startnr", type="int", default = 0,
                      help="number of first file [%default]")
    parser.add_option("-N", "--nr_of_files", dest="nr_of_files", type="int", default = 50,
                      help="total number of files [%default]")



    options, args = parser.parse_args()
    return options.__dict__

if __name__ == '__main__':
    options = parse_arguments()

    rdc_TPS(**options)

