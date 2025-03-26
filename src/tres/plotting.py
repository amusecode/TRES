import matplotlib.pyplot as plt
import amuse.plot as aplt
from amuse.units import units
import numpy as np

class plot_data_container():
    def __init__(self):
        return

def plot_function(triple, dir_plots):
    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
    t_max_Myr = max(times_array_Myr)
    a_in_array_AU = triple.plot_data.a_in_array.value_in(units.AU)
    g_in_array = triple.plot_data.g_in_array
    e_in_array = triple.plot_data.e_in_array
    i_relative_array = triple.plot_data.i_relative_array
    o_in_array = triple.plot_data.o_in_array
    a_out_array_AU = triple.plot_data.a_out_array.value_in(units.AU)
    g_out_array = triple.plot_data.g_out_array
    e_out_array = triple.plot_data.e_out_array
    o_out_array = triple.plot_data.o_out_array
    m1_array = triple.plot_data.m1_array.value_in(units.MSun)
    m2_array = triple.plot_data.m2_array.value_in(units.MSun)
    m3_array = triple.plot_data.m3_array.value_in(units.MSun)
    spin1_array = triple.plot_data.spin1_array
    spin2_array = triple.plot_data.spin2_array
    spin3_array = triple.plot_data.spin3_array
    r1_array = triple.plot_data.r1_array
    r2_array = triple.plot_data.r2_array
    r3_array = triple.plot_data.r3_array
    moi1_array = triple.plot_data.moi1_array
    moi2_array = triple.plot_data.moi2_array
    moi3_array = triple.plot_data.moi3_array
    RL1_array = triple.plot_data.RL1_array
    RL2_array = triple.plot_data.RL2_array
    RL3_array = triple.plot_data.RL3_array
    delta_e_in_array = triple.plot_data.delta_e_in_array
    
    f = open(triple.file_name[:-4]+'.txt','w')
    f.write('#' + str(t_max_Myr) + '\n')
    for i_p in range(len(times_array_Myr)):
        f.write(str(times_array_Myr[i_p]) + '\t')
        f.write(str(m1_array[i_p] ) + '\t')
        f.write(str(m2_array[i_p] ) + '\t')
        f.write(str(m3_array[i_p] ) + '\t')
        f.write(str(spin1_array[i_p] ) + '\t')
        f.write(str(spin2_array[i_p] ) + '\t')
        f.write(str(spin3_array[i_p] ) + '\t')
        f.write(str(r1_array[i_p] ) + '\t')
        f.write(str(r2_array[i_p] ) + '\t')
        f.write(str(r3_array[i_p] ) + '\t')
        f.write(str(RL1_array[i_p] ) + '\t')
        f.write(str(RL2_array[i_p] ) + '\t')
        f.write(str(RL3_array[i_p] ) + '\t')
        f.write(str(moi1_array[i_p] ) + '\t')
        f.write(str(moi2_array[i_p] ) + '\t')
        f.write(str(moi3_array[i_p] ) + '\t')
        f.write(str(i_relative_array[i_p] ) + '\t')
        f.write(str(a_in_array_AU[i_p] ) + '\t')
        f.write(str(g_in_array[i_p] ) + '\t')
        f.write(str(e_in_array[i_p] ) + '\t')
        f.write(str(o_in_array[i_p] ) + '\t')
        f.write(str(a_out_array_AU[i_p] ) + '\t')
        f.write(str(g_out_array[i_p] ) + '\t')
        f.write(str(e_out_array[i_p] ) + '\t')
        f.write(str(o_out_array[i_p] ) + '\t')
        f.write(str(delta_e_in_array[i_p] ) + '\t')
        f.write('\n')
    f.close()
    
#    for i_s in range(len(times_array_Myr)):
#        print(times_array_Myr[i_s], end = ' ')
#        print( a_in_array_AU[i_s], end = ' ')
#        print( g_in_array[i_s], end = ' ')
#        print( e_in_array[i_s], end = ' ')
#        print( i_relative_array[i_s], end = ' ')
#        print( o_in_array[i_s], end = ' ')
#        print( a_out_array_AU[i_s], end = ' ')
#        print( g_out_array[i_s], end = ' ')
#        print( e_out_array[i_s], end = ' ')
#        print( o_out_array[i_s], end = ' ')
#        print( m1_array[i_s], end = ' ')
#        print( m2_array[i_s], end = ' ')
#        print( m3_array[i_s], end = ' ')
#        print( spin1_array[i_s], end = ' ')
#        print( spin2_array[i_s], end = ' ')  
#        print( spin3_array[i_s], end = ' ')  
#        print( r1_array[i_s], end = ' ')    
#        print( r2_array[i_s], end = ' ')    
#        print( r3_array[i_s], end = ' ')
#        print( RL1_array[i_s], end = ' ')
#        print( RL2_array[i_s], end = ' ')
#        print( RL3_array[i_s], end = ' ')
#        print( moi1_array[i_s], end = ' ')   
#        print( moi2_array[i_s], end = ' ')   
#        print( moi3_array[i_s], end = ' ')
#        print( delta_e_in_array[i_s])
    
    
            
    ### plots to test secular code ###
    import amuse.plot as aplt
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 16})    
    
#    generic_name = '_M'+str(m1_array[0]) + '_m'+str(m2_array[0]) +'_n'+str(m3_array[0]) + '_a'+str(a_in_array_AU[0]) + '_A'+str(a_out_array_AU[0]) + '_e'+str(e_in_array[0]) + '_E'+str(e_out_array[0]) + '_i'+str(i_relative_array[0]/np.pi*180.0) + '_g'+str(g_in_array[0]) + '_G'+str(g_out_array[0]) + '_o'+str(o_in_array[0]) + '_O'+str(o_out_array[0]) + '_t'+str(t_max_Myr) + '_maxdr'+str(triple.maximum_radius_change_factor)+'_edr'+str(error_dr)
    generic_name = ''

    figure = plt.figure(figsize=(10,13), tight_layout=True)
    N_subplots = 4
    
    plot_e = figure.add_subplot(N_subplots,1,1)
    plot_i_relative = figure.add_subplot(N_subplots,1,2)
    plot_a_in = figure.add_subplot(N_subplots,1,3)
    plot_a_out = figure.add_subplot(N_subplots,1,4)
    
    plot_e.plot(times_array_Myr,e_in_array, label= '$e_\mathrm{in}$')
    plot_e.plot(times_array_Myr,e_out_array, label= '$e_\mathrm{out}$')
    plot_e.set_xlim(0,t_max_Myr)
    plot_e.set_xlabel('$t/\mathrm{Myr}$')
    plot_e.set_ylabel('$e$')
    plot_e.legend(loc=0)
    
    plot_i_relative.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative.set_xlim(0,t_max_Myr)
    plot_i_relative.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')
    
    plot_a_in.plot(times_array_Myr,a_in_array_AU)
    plot_a_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_in.set_ylabel('$a_\mathrm{in}$')

    plot_a_out.plot(times_array_Myr,a_out_array_AU)
    plot_a_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_out.set_ylabel('$a_\mathrm{out}$')
    
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)
    plt.savefig(dir_plots+'TRES'+generic_name+'.pdf')
#    plt.show()
    plt.close()




    figure = plt.figure(figsize=(10,13), tight_layout=True)
    N_subplots = 4
    
    plot_e_in = figure.add_subplot(N_subplots,1,1)
    plot_e_out = figure.add_subplot(N_subplots,1,2)
    plot_i_relative = figure.add_subplot(N_subplots,1,3)
    plot_a_in = figure.add_subplot(N_subplots,1,4)
    
    plot_e_in.plot(times_array_Myr,e_in_array, label= '$e_\mathrm{in}$')
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_e_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_in.set_ylabel('$e$')
    plot_e_in.legend(loc=0)

    plot_e_out.plot(times_array_Myr,e_out_array, label= '$e_\mathrm{out}$')
    plot_e_out.set_xlim(0,t_max_Myr)
    plot_e_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_out.set_ylabel('$e$')
    plot_e_out.legend(loc=0)

    
    plot_i_relative.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative.set_xlim(0,t_max_Myr)
    plot_i_relative.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')
    
    plot_a_in.plot(times_array_Myr,a_in_array_AU)
    plot_a_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_in.set_ylabel('$a_\mathrm{in}$')

    
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)
    plt.savefig(dir_plots+'TRES2'+generic_name+'.pdf')
#    plt.show()
    plt.close()





    figure = plt.figure(figsize=(10,13), tight_layout=True)
    N_subplots = 4

    plot_e_in = figure.add_subplot(N_subplots,1,1)
    plot_i_relative = figure.add_subplot(N_subplots,1,2)
    plot_e_in_g_in = figure.add_subplot(N_subplots,1,3)
    plot_a_in = figure.add_subplot(N_subplots,1,4)


    plot_e_in.plot(times_array_Myr,e_in_array)
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_e_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_in.set_ylabel('$e_\mathrm{in}$')

    plot_i_relative.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative.set_xlim(0,t_max_Myr)
    plot_i_relative.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')

    plot_e_in_g_in.plot(np.cos(g_in_array),e_in_array)
    plot_e_in_g_in.set_xlabel('$\cos(g_\mathrm{in})$')
    plot_e_in_g_in.set_ylabel('$e_\mathrm{in}$')

    plot_a_in.plot(times_array_Myr,a_in_array_AU)
    plot_a_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_in.set_ylabel('$a_\mathrm{in}$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    plt.savefig(dir_plots+'TRES_inner_orbit'+generic_name+'.pdf')
#    plt.show()
    plt.close()


    dyn_inst =  2.8 / (1-e_out_array) * (1-0.3*i_relative_array/np.pi) * ((1+m3_array/(m1_array+m2_array))*(1+e_out_array)/(np.sqrt(1-e_out_array)))**0.4 / (a_out_array_AU/a_in_array_AU)
    oct = (m1_array-m2_array)/(m1_array+m2_array) * a_in_array_AU/a_out_array_AU * e_out_array/(1-e_out_array**2)
    semiseq = (5.*np.pi* m3_array/(m1_array+m2_array) * (a_in_array_AU/a_out_array_AU/(1-e_out_array))**3)/ np.sqrt(1-e_in_array) 

    plt.semilogy(times_array_Myr, dyn_inst, '.')
    plt.semilogy(times_array_Myr, oct, '.')
    plt.semilogy(times_array_Myr, semiseq, '.')
    plt.xlabel('t (Myr)')
    plt.ylabel('stability,oct and semiseq')
    plt.legend(loc=0)
    plt.savefig(dir_plots+'dyn_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    plt.plot(times_array_Myr,e_in_array)
    plt.plot(times_array_Myr,e_in_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$e_\mathrm{in}$')
    plt.savefig(dir_plots+'e_in_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


    plt.plot(times_array_Myr,e_in_array)
    plt.plot(times_array_Myr,e_in_array, '.')
    plt.plot(times_array_Myr,e_out_array)
    plt.plot(times_array_Myr,e_out_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$e_\mathrm{in}$')
    plt.savefig(dir_plots+'e_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()
    
    plt.plot(times_array_Myr,g_in_array)
    plt.plot(times_array_Myr,g_in_array, '.')
    plt.plot(times_array_Myr,g_out_array)
    plt.plot(times_array_Myr,g_out_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$g_\mathrm{in}$')
    plt.savefig(dir_plots+'g_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()
    
    
    plt.plot(e_in_array,(g_in_array%np.pi)/np.pi*180)
    plt.plot(e_in_array,(g_in_array%np.pi)/np.pi*180, '.')
#    plt.plot(e_out_array,g_out_array)
#    plt.plot(e_out_array,g_out_array, '.')
#    plt.xlim(0,1)
    plt.xlabel('$e_\mathrm{in}$')
    plt.ylabel('$g_\mathrm{in}$')
    plt.savefig(dir_plots+'g_e_inner'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    

    plt.plot(times_array_Myr,o_in_array)
    plt.plot(times_array_Myr,o_in_array, '.')
    plt.plot(times_array_Myr,o_out_array)
    plt.plot(times_array_Myr,o_out_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$o_\mathrm{in}$')
    plt.savefig(dir_plots+'o_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


#    a_in_final_theory =  a_in_array_AU[0] * (m1_array[0] * m2_array[0] / m1_array / m2_array)**2 #stable mt
    a_in_final_theory =  a_in_array_AU[0] * (m1_array[0] + m2_array[0]) / (m1_array + m2_array)#wind 
    plt.plot(times_array_Myr,a_in_array_AU)
    plt.plot(times_array_Myr,a_in_array_AU, '.')
    plt.plot(times_array_Myr,a_in_final_theory)
    plt.plot(times_array_Myr,a_in_final_theory, '.')
    
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$a_\mathrm{in}$')
    plt.savefig(dir_plots+'semi_in_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()



    constants = 6283728.92847 # constants.G #in au and solar mass per Myr
    corot_spin_inner = 1./np.sqrt(a_in_array_AU**3/(m1_array+m2_array))*constants
    corot_spin_outer = 1./np.sqrt(a_out_array_AU**3/(m1_array+m2_array+m3_array))*constants
    RSun_in_AU = 0.00464913034382
    critical_rot_1 = constants*np.sqrt(m1_array/(RSun_in_AU*r1_array)**3)
    critical_rot_2 = constants*np.sqrt(m2_array/(RSun_in_AU*r2_array)**3)
    critical_rot_3 = constants*np.sqrt(m3_array/(RSun_in_AU*r3_array)**3)

    plt.plot(times_array_Myr,spin1_array, 'b-')
    plt.plot(times_array_Myr,spin1_array, 'b.')
    plt.plot(times_array_Myr,spin2_array, 'g-')
    plt.plot(times_array_Myr,spin2_array, 'g.')
    plt.plot(times_array_Myr,spin3_array, 'r-')
    plt.plot(times_array_Myr,spin3_array, 'r.')

    plt.plot(times_array_Myr,corot_spin_inner, 'c-')
    plt.plot(times_array_Myr,corot_spin_inner, 'c,')
    plt.plot(times_array_Myr,corot_spin_outer, 'm-')
    plt.plot(times_array_Myr,corot_spin_outer, 'm,')

    plt.plot(times_array_Myr,critical_rot_1, 'y-')
    plt.plot(times_array_Myr,critical_rot_1, 'y,')
    plt.plot(times_array_Myr,critical_rot_2, 'k-')
    plt.plot(times_array_Myr,critical_rot_2, 'k,')
    plt.plot(times_array_Myr,critical_rot_3, 'k-')
    plt.plot(times_array_Myr,critical_rot_3, 'ko')

    
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$spin$')
    plt.savefig(dir_plots+'spin_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()



#        ms = 1|units.MSun
#        rs = 1|units.AU
#        J2 = ms**3 * rs * constants.G
#        J = np.sqrt(J2)
#        print(J)  #2.9071938904e+11 [RSun**2 * MSun * Myr**-1]       

#        rs = 1|units.RSun
#        print(J)  #19822565357. [RSun**2 * MSun * Myr**-1]   

    constants_Jorb = 2.9071938904e+11 #[RSun**2 * MSun * Myr**-1]
    J_orb2 = m1_array**2 * m2_array**2 / (m1_array+m2_array) * a_in_array_AU * ( 1-e_in_array**2)
    J_orb = np.sqrt(J_orb2)*constants_Jorb

    J_spin1 =  spin1_array * moi1_array
    J_spin2 =  spin2_array * moi2_array
    J_spin3 =  spin3_array * moi3_array

    plt.plot(times_array_Myr, J_orb)
    plt.plot(times_array_Myr,J_orb, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J orbit$')
    plt.savefig(dir_plots+'Jorbit_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    plt.plot(times_array_Myr,J_spin1)
    plt.plot(times_array_Myr,J_spin1, '.')
    plt.plot(times_array_Myr,J_spin2)
    plt.plot(times_array_Myr,J_spin2, '.')
    plt.plot(times_array_Myr,J_spin3)
    plt.plot(times_array_Myr,J_spin3, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J spin$')
    plt.savefig(dir_plots+'Jspin_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


    plt.plot(times_array_Myr, J_orb)
    plt.plot(times_array_Myr,J_orb, '.')
    plt.plot(times_array_Myr,J_spin1)
    plt.plot(times_array_Myr,J_spin1, '.')
    plt.plot(times_array_Myr,J_spin2)
    plt.plot(times_array_Myr,J_spin2, '.')
    plt.plot(times_array_Myr,J_spin3)
    plt.plot(times_array_Myr,J_spin3, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J spin$')
    plt.savefig(dir_plots+'Js_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    plt.semilogy(times_array_Myr,moi1_array)
    plt.semilogy(times_array_Myr,moi1_array, '.')
    plt.semilogy(times_array_Myr,moi2_array)
    plt.semilogy(times_array_Myr,moi2_array, '.')
    plt.semilogy(times_array_Myr,moi3_array)
    plt.semilogy(times_array_Myr,moi3_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$moi$')
    plt.savefig(dir_plots+'moment_of_inertia_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()

    

    plt.semilogy(times_array_Myr,r1_array)
    plt.semilogy(times_array_Myr,r1_array, '.')
    plt.semilogy(times_array_Myr,r2_array)
    plt.semilogy(times_array_Myr,r2_array, '.')
    plt.semilogy(times_array_Myr,r3_array)
    plt.semilogy(times_array_Myr,r3_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$radius$')
    plt.savefig(dir_plots+'radius_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()
    
    dr1_array =r1_array[1:]-r1_array[:-1]
    dr2_array =r2_array[1:]-r2_array[:-1]
    dr3_array =r3_array[1:]-r3_array[:-1]
    dt_array = times_array_Myr[1:] - times_array_Myr[:-1]
    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array)
    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array, 'b.')
    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array*-1., 'b*')
    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array)
    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array, 'g.')
    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array*-1., 'g*')
    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array)
    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array, 'r.')
    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array*-1., 'r*')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$dr/dt$')
    plt.savefig(dir_plots+'drdt_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()




    plt.semilogy(times_array_Myr[1:], dr1_array/r1_array[1:])
    plt.semilogy(times_array_Myr[1:], dr1_array/r1_array[1:], 'b.')
    plt.semilogy(times_array_Myr[1:], -1.*dr1_array/r1_array[1:], 'b*')
    plt.semilogy(times_array_Myr[1:], dr2_array/r2_array[1:])
    plt.semilogy(times_array_Myr[1:], dr2_array/r2_array[1:], 'g.')
    plt.semilogy(times_array_Myr[1:], -1.*dr2_array/r2_array[1:], 'g*')
    plt.semilogy(times_array_Myr[1:], dr3_array/r3_array[1:])
    plt.semilogy(times_array_Myr[1:], dr3_array/r3_array[1:], 'r.')
    plt.semilogy(times_array_Myr[1:], -1.*dr3_array/r3_array[1:], 'r*')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$dr/r$')
    plt.savefig(dir_plots+'dr_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()



    dm1_array =m1_array[:-1]-m1_array[1:]
    dm2_array =m2_array[:-1]-m2_array[1:]
    dm3_array =m3_array[:-1]-m3_array[1:]
    
    plt.plot(times_array_Myr[1:], m1_array[1:])
    plt.plot(times_array_Myr[1:], m1_array[1:], '.')
    plt.plot(times_array_Myr[1:], m2_array[1:])
    plt.plot(times_array_Myr[1:], m2_array[1:], '.')
    plt.plot(times_array_Myr[1:], m3_array[1:])
    plt.plot(times_array_Myr[1:], m3_array[1:], '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$mass$')
    plt.savefig(dir_plots+'mass_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


#    plt.semilogy(times_array_Myr[1:], dm1_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm1_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dm2_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm2_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dm3_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm3_array/dt_array, '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dm/dt$')
#    plt.savefig(dir_plots+'dmdt_time'+generic_name+'.pdf')
#    plt.show()
#    plt.close()

#
#    plt.semilogy(times_array_Myr[1:], dm1_array/m1_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm1_array/m1_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dm2_array/m2_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm2_array/m2_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dm3_array/m3_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm3_array/m3_array[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dm/m$')
#    plt.savefig(dir_plots+'dm_time'+generic_name+'.pdf')
#    plt.show()
#    plt.close()

 
    plt.semilogy(times_array_Myr[1:], dt_array)
    plt.semilogy(times_array_Myr[1:], dt_array, '.')
    plt.semilogy(times_array_Myr[1:], dt_array)
    plt.semilogy(times_array_Myr[1:], dt_array, '.')
    plt.semilogy(times_array_Myr[1:], dt_array)
    plt.semilogy(times_array_Myr[1:], dt_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$dt$')
    plt.savefig(dir_plots+'dt_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


    RL1_frac = r1_array/RL1_array
    RL2_frac = r2_array/RL2_array
    RL3_frac = r3_array/RL3_array
    plt.plot(times_array_Myr,RL1_frac)
    plt.plot(times_array_Myr,RL1_frac, '.', label='Primary')
    plt.plot(times_array_Myr,RL2_frac)
    plt.plot(times_array_Myr,RL2_frac, '.', label='Secondary')
    plt.plot(times_array_Myr,RL3_frac)
    plt.plot(times_array_Myr,RL3_frac, '.', label='Tertiary')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('r/RL')
    plt.legend(fontsize=12)
    plt.savefig(dir_plots+'rl_frac_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()
    
    
    plt.plot(times_array_Myr,delta_e_in_array)
    plt.plot(times_array_Myr,delta_e_in_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$\Delta e_{\mathrm{in}}$')
    plt.savefig(dir_plots+'delta_e_in_time'+generic_name+'.pdf')
#    plt.show()
    plt.close()


#   wind a = ai * Mti/Mt
#    Mtot = m1_array+m2_array    
#    plt.plot(times_array_Myr,a_in_array_AU)
#    plt.plot(times_array_Myr,a_in_array_AU, '.')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot, '.')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:])
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{in}$')
#    plt.savefig(dir_plots+'semi_inner_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_inner_rel_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()


#   cons mt a = ai * (m1i*m2i*/m1/m2)**2
#    plt.plot(times_array_Myr,a_in_array_AU, 'b-')
#    plt.plot(times_array_Myr,a_in_array_AU, 'b.')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2, 'g-')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2, 'g.')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*(m1_array[:-1]*m2_array[:-1]/m1_array[1:]/m2_array[1:])**2, 'r-')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*(m1_array[:-1]*m2_array[:-1]/m1_array[1:]/m2_array[1:])**2, 'r.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{in}$')
#    plt.savefig(dir_plots+'semi_inner_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_inner_rel_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()



#    dm = (m1_array[1:] - m1_array[:-1] )
#    dt = (times_array_Myr[1:] - times_array_Myr[:-1])
#    dmdt = (m1_array[1:] - m1_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])

#    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
    
    
    
    
    
#    plt.plot(times_array_Myr, g_in_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$g_\mathrm{in}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(g_in_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(g_\mathrm{in})$')
#    plt.show()
#    plt.close()
#
#
#    plt.plot(times_array_Myr, o_in_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$o_\mathrm{in}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(o_in_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(o_\mathrm{in})$')
#    plt.show()
#    plt.close()


    #outer binary
    figure = plt.figure(figsize=(10,13))
    N_subplots = 4

    plot_e_out = figure.add_subplot(N_subplots,1,1)
    plot_i_relative2 = figure.add_subplot(N_subplots,1,2)
    plot_e_out_g_out = figure.add_subplot(N_subplots,1,3)
    plot_a_out = figure.add_subplot(N_subplots,1,4)

#    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
#    t_max_Myr = max(times_array_Myr)
#    i_relative_array = triple.plot_data.i_relative_array

    plot_e_out.plot(times_array_Myr,e_out_array)
    plot_e_out.set_xlim(0,t_max_Myr)
    plot_e_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_out.set_ylabel('$e_\mathrm{out}$')

    plot_i_relative2.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative2.set_xlim(0,t_max_Myr)
    plot_i_relative2.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative2.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative2.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')

    plot_e_out_g_out.plot(np.cos(g_out_array),e_out_array)
    plot_e_out_g_out.set_xlabel('$\cos(g_\mathrm{out})$')
    plot_e_out_g_out.set_ylabel('$e_\mathrm{out}$')

    plot_a_out.plot(times_array_Myr,a_out_array_AU)
    plot_a_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_out.set_ylabel('$a_\mathrm{out}$')

    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)
    plt.savefig(dir_plots+'TRES_outer_orbit'+generic_name+'.pdf')
#    plt.show()
    plt.close()



#    plt.plot(times_array_Myr,e_out_array)
#    plt.xlim(0,t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$e_\mathrm{out}$')
#    plt.savefig(dir_plots+'e_out_time'+generic_name+'.pdf')
#    plt.show()
#    plt.close()


#    Mtott = m1_array+m2_array+m3_array    
#    plt.plot(times_array_Myr,a_out_array_AU)
#    plt.plot(times_array_Myr,a_out_array_AU, '.')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott, '.')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:])
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{out}$')
#    plt.savefig(dir_plots+'semi_outer_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_outer_rel_wind'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#
#    m_in_array = m1_array+m2_array
#    plt.plot(times_array_Myr,a_out_array_AU, 'b-')
#    plt.plot(times_array_Myr,a_out_array_AU, 'b.')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2, 'g-')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2, 'g.')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*(m_in_array[:-1]*m3_array[:-1]/m_in_array[1:]/m3_array[1:])**2, 'r-')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*(m_in_array[:-1]*m3_array[:-1]/m_in_array[1:]/m3_array[1:])**2, 'r.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{out}$')
#    plt.savefig(dir_plots+'semi_outer_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig(dir_plots+'semi_outer_rel_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    plt.close()
#

#    dm = (m3_array[1:] - m3_array[:-1] )
#    dmdt = (m3_array[1:] - m3_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])

#    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#    plt.close()
#
#    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#    plt.close()



#    plt.plot(np.cos(g_out_array), e_in_array)
#    plt.xlabel('$\cos(g_\mathrm{out})$')
#    plt.ylabel('$e_\mathrm{in}$')
#    plt.show()
#    plt.close()
#
#    plt.plot(np.cos(g_in_array), np.cos(g_out_array))
#    plt.xlabel('$\cos(g_\mathrm{in})$')
#    plt.ylabel('$\cos(g_\mathrm{out})$')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, g_out_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$g_\mathrm{out}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(g_out_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(g_\mathrm{out})$')
#    plt.show()
#    plt.close()
#
#    plt.plot(times_array_Myr, o_out_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$o_\mathrm{out}$')
#    plt.show()
#    plt.close()

#    plt.plot(times_array_Myr, np.cos(o_out_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(o_\mathrm{out})$')
#    plt.show()
#    plt.close()
#
#    aplt.plot(times_array_Myr, m1_array)
#    aplt.plot(times_array_Myr, m1_array, '.')
#    aplt.plot(times_array_Myr, m2_array)
#    aplt.plot(times_array_Myr, m2_array, '.')
#    aplt.plot(times_array_Myr, m3_array)
#    aplt.plot(times_array_Myr, m3_array, '.')
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$M/\mathrm{MSun}$')
#    plt.show()
#    plt.close()
#    
    