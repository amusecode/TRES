# TPS:         Triple population synthesis.
##              computes the evolution of a population of given triples   
##              given any initial conditions (M, q, m, A, a, E, e, i, G, g, O, o).
## Options:   --M_max    upper limit for the inner primary mass [100 Msun]
##            --M_min    lower limit for the inner primary mass [0.1 Msun]
##            --M_distr  mass function option: 
lib_inner_primary_mass_distr = {0: "Kroupa", #default
                            1: "Scalo",
                            2: "Miller & Scalo",
                            3: "Salpeter",
                            4: "Logarithmically flat",
                            5: "Eggleton",
                            6: "Kroupa for massive stars M>0.5 powerlaw with exp=-2.3",}
##            --Q_max    upper limit for the inner mass ratio [1.]
##            --Q_min    lower limit for the inner mass ratio [0.]
##            --Q_distr  inner mass ratio option: 
lib_inner_mass_ratio_distr = {0: "Uniform distribution", #default
                            1: "Kroupa IMF", #draws from mass distribution instead of mass ratio distribution
                            2: "Galicher et al. 2016 powerlaw (M^-1.31)", } #draws from mass distribution instead of mass ratio distribution, appropriate for planets
##            --q_max    upper limit for the outer mass ratio [1.]
##            --q_min    lower limit for the mass of the outer star [0.]
##            --q_distr  outer mass ratio option: 
lib_outer_mass_ratio_distr = {0: "Uniform distribution", #default
                            1: "Kroupa IMF", #draws from mass distribution instead of mass ratio distribution
                            2: "Galicher et al. 2016 powerlaw (M^-1.31)", } #draws from mass distribution instead of mass ratio distribution, appropriate for planets
##            --A_max    upper limit for the inner semi-major axis [5e6 RSun]
##            --A_min    lower limit for the inner semi-major axis [5]
##            --A_distr  inner semi-major axcis option: 
lib_inner_semi_distr = {0: "Log Uniform distribution", #default
                   1: "Constant semi-major axis",
                   2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
                   3: "Lognormal mu = 10^3.5d, sigma = 2.3",
                   4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
                   5: "Sana et al. 2012",
                   6: "flat distribution",
                   7: "Galicher et al. 2016 powerlaw (a^-0.61)",} #appropriate for planets
##            --a_max    upper limit for the outer semi-major axis [5e6 RSun]
##            --a_min    lower limit for the outer semi-major axis [5 RSun]
##            --a_distr  outer semi-major axis option: 
lib_outer_semi_distr = {0: "Log Uniform distribution", #default
                   1: "Constant semi-major axis",
                   2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
                   3: "Lognormal mu = 10^3.5d, sigma = 2.3",
                   4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
                   5: "Sana et al. 2012",
                   6: "flat distribution",
                   7: "Galicher et al. 2016 powerlaw (a^-0.61)",} #appropriate for planets
##            --E_max    upper limit for the inner eccentricity [0.9]
##            --E_min    lower limit for the inner eccentricity [0.]
##            --E_distr  inner eccentricity option: 
lib_inner_ecc_distr = {0: "Thermal", #default
                 1: "Constant eccentricity",
                 2: "Sana et al. 2012 e^-0.45", #-> close binaries
                 3: "Flat distribution",
                 4: "Powerlaw e^0.5",
                 5: "Bowler et al. 2020 Beta distribution",}  #appropriate for planets                                  
##            --e_max    upper limit for the outer eccentricity [0.9]
##            --e_min    lower limit for the outer eccentricity [0.]
##            --e_distr  outer eccentricity option: 
lib_outer_ecc_distr = {0: "Thermal", #default
                 1: "Constant eccentricity",
                 2: "Sana et al. 2012 e^-0.45", #-> close binaries
                 3: "Flat distribution",
                 4: "Powerlaw e^0.5",
                 5: "Bowler et al. 2020 Beta distribution",}  #appropriate for planets                                                      
##            --i_max    upper limit for the relative inclination [pi]
##            --i_min    lower limit for the relative inclination [0]
##            --i_distr  relative inclination option: 
lib_incl_distr = {0: "Circular uniform distribution", #default
                 1: "Constant inclination",}     
##            --G_max    upper limit for the inner argument of pericenter [pi]
##            --G_min    lower limit for the inner argument of pericenter [-pi]
##            --G_distr  inner argument of pericenter option: r
lib_inner_aop_distr = {0: "Uniform distribution", #default
                 1: "Constant argument of pericenter",}     
##            --g_max    upper limit for the outer argument of pericenter [pi]
##            --g_min    lower limit for the outer argument of pericenter [-pi]
##            --g_distr  outer argument of pericenter option: 
lib_outer_aop_distr = {0: "Uniform distribution", #default
                 1: "Constant argument of pericenter",}     
##             outer longitude of ascending nodes = inner - pi               
##            --O_max    upper limit for the inner longitude of ascending node [pi]
##            --O_min    lower limit for the inner longitude of ascending node [-pi]
##            --O_distr  inner longitude of ascending node option: 
lib_inner_loan_distr = {0: "Circular niform distribution", 
                 1: "Constant longitude of ascending nodes",} #default
##            -T or -t   binary end time. [13500 Myr]
##            -z         metallicity of stars  [0.02 Solar] 
##            -n         total number of systems to be simulated.  [1]
##            -N         number ID of first system.  [0]
##            --no_stop_at_mass_transfer                    stopping condition at mass transfer 
##            --no_stop_at_init_mass_transfer               stopping condition at mass transfer at initialisation
##            --no_stop_at_outer_mass_transfer              stopping condition at mass transfer in outer binary 
##            --stop_at_stable_mass_transfer                stopping condition at stable mass transfer 
##            --stop_at_eccentric_stable_mass_transfer      stopping condition at eccentric stable mass transfer 
##            --stop_at_unstable_mass_transfer              stopping condition at unstable mass transfer 
##            --stop_at_eccentric_unstable_mass_transfer    stopping condition at eccentric unstable mass transfer 
##            --stop_at_no_CHE                              stopping condition if no chemically homogeneous evolution 

##            --no_stop_at_merger                           stopping condition at merger 
##            --no_stop_at_disintegrated                    stopping condition at disintegration 
##            --no_stop_at_inner_collision                  stopping condition at collision in inner binary 
##            --no_stop_at_outer_collision                  stopping condition at collision involving tertiary star 
##            --no_stop_at_dynamical_instability            stopping condition at dynamical instability 
##            --stop_at_semisecular_regime                  stopping condition at semisecular regime 
##            --stop_at_SN                                  stopping condition at supernova 
lib_SN_kick_distr = {0: "No kick",
                    1: "Hobbs", #Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
                    2: "Arzoumanian", #Arzoumanian ea 2002, 568, 289
                    3: "Hansen", #Hansen & Phinney 1997, 291, 569
                    4: "Paczynski", #Paczynski 1990, 348, 485
                    5: "Verbunt", #Verbunt, Igoshev & Cator, 2017, 608, 57
                    } #default


lib_CE = {  0: "alpha-ce + alpha-dce",
            1: "gamma-ce + alpha-dce", 
            2: "seba style; combination of gamma-ce, alpha-ce & alpha-dce", 
}


         
#not implemented yet
##            -s         random seed


import TRES as TRES
from amuse.community.seba.interface import SeBa
from seculartriple_TPS.interface import SecularTriple

secular_code = SecularTriple()
import sys
from amuse.units.optparse import OptionParser
from amuse.units import units, constants
from amuse.support.console import set_printing_strategy
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import beta as beta_distribution

from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.ic.scalo import new_scalo_mass_distribution
from amuse.ic.millerscalo import new_miller_scalo_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.flatimf import new_flat_mass_distribution

from TRES_options import REPORT_TPS, \
                         REPORT_USER_WARNINGS_TPS, \
                         EXCLUDE_SSO, \
                         precision,  min_mass, absolute_min_mass, absolute_max_mass                         

def flat_distr(lower, upper):
    return np.random.uniform(lower, upper)

def log_flat_distr(lower, upper):
    x= np.random.uniform(np.log10(lower), np.log10(upper))
    return 10**x
    
def eggleton_mass_distr(lower_mass, upper_mass):
    turnover_mass = 0.3|units.MSun
    power = 0.85
    
    y_max = (upper_mass/turnover_mass) **(1/power)
    upper = y_max / (1+y_max)
    y_min = (lower_mass/turnover_mass) **(1/power)
    lower = y_min / (1+y_min)

    x = np.random.uniform(lower, upper)
    y=turnover_mass * (x/(1-x))**power
    return turnover_mass * (x/(1-x))**power

def powerlaw_distr(m_min, m_max, slope):    
    if slope == -1 or slope == 0:
        sys.exit('slope of powerlap distribution incorrect')
    slope1 = slope + 1
    factor = (pow(m_max / m_min, slope1) - 1.0 )
    x = np.random.uniform(0,1)
    return m_min * (1.0 + factor*x) ** (1.0 / slope1)
    
    
def beta_distr_SSOs(lower, upper, mass):    # (Bowler et al. 2020)
    min_mass_BD = 16 |units.MJupiter        # brown dwarf boundary
    if absolute_min_mass < mass <= min_mass_BD :
        A, B = 30, 200      # for exoplanets 
    elif min_mass_BD < mass < min_mass :
        A, B = 2.30, 1.65   # for brown dwarfs 
    else:
        sys.exit('You may be using a SSOs distribution for a stellar object, exiting')

    e_sample = beta_distribution.rvs( A, B)
    if lower <= e_sample <= upper:
        return e_sample
        
    return beta_distr_SSOs(lower, upper, mass)    # pick another sample within given bounds
    

class Generate_initial_triple:
    #-------
    #setup stellar system
    def __init__(self, inner_primary_mass_max, inner_primary_mass_min, 
                        inner_secondary_mass_max, inner_secondary_mass_min, 
                        outer_mass_min, outer_mass_max,
                        inner_mass_ratio_max, inner_mass_ratio_min, 
                        outer_mass_ratio_max, outer_mass_ratio_min, 
                        inner_semi_max, inner_semi_min, outer_semi_max, outer_semi_min, 
                        inner_semi_latus_rectum_min, outer_semi_latus_rectum_min,
                        inner_semi_latus_rectum_max, outer_semi_latus_rectum_max,
                        inner_ecc_max, inner_ecc_min, outer_ecc_max, outer_ecc_min, 
                        incl_max, incl_min,
                        inner_aop_max, inner_aop_min, outer_aop_max, outer_aop_min,
                        inner_loan_max, inner_loan_min, 
                        inner_primary_mass_distr, inner_mass_ratio_distr, outer_mass_ratio_distr,
                        inner_semi_distr,  outer_semi_distr, inner_ecc_distr, outer_ecc_distr, incl_distr,
                        inner_aop_distr, outer_aop_distr, inner_loan_distr):
                        
                        if inner_primary_mass_distr == 5:
                            mass_convergence = False
                            while mass_convergence == False:
                                mass_convergence = self.generate_mass_and_semi_eggleton(inner_primary_mass_max, inner_primary_mass_min,        
                                inner_secondary_mass_min,outer_mass_min, outer_mass_max,
                                inner_semi_max, inner_semi_min, outer_semi_max, outer_semi_min)                            
    
                            #Does not use boolean inner/outer _semi_latus_rectum_ min/max
                            self.inner_ecc = self.generate_ecc_1d(inner_ecc_max, inner_ecc_min, inner_ecc_distr, self.inner_secondary_mass)
                            self.outer_ecc = self.generate_ecc_1d(outer_ecc_max, outer_ecc_min, outer_ecc_distr, self.outer_mass)
                                

                        else:    
                            self.generate_mass(inner_primary_mass_max, inner_primary_mass_min, 
                                inner_secondary_mass_max,inner_secondary_mass_min,outer_mass_min,outer_mass_max,
                                inner_mass_ratio_max, inner_mass_ratio_min,
                                outer_mass_ratio_max, outer_mass_ratio_min, 
                                inner_primary_mass_distr, inner_mass_ratio_distr, outer_mass_ratio_distr)
                                

                            orbit_convergence = False
                            while orbit_convergence == False:
                                orbit_convergence =  self.generate_semi_and_ecc(inner_semi_max, inner_semi_min, 
                                    outer_semi_max, outer_semi_min,
                                    inner_semi_distr,  outer_semi_distr,
                                    inner_semi_latus_rectum_min, outer_semi_latus_rectum_min, 
                                    inner_semi_latus_rectum_max, outer_semi_latus_rectum_max, 
                                    inner_ecc_max, inner_ecc_min, 
                                    outer_ecc_max, outer_ecc_min,
                                    inner_ecc_distr, outer_ecc_distr, 
                                    self.inner_secondary_mass, self.outer_mass)


                        self.generate_incl(incl_max, incl_min, incl_distr)

                        self.generate_aop(inner_aop_max, inner_aop_min, 
                            outer_aop_max, outer_aop_min,
                            inner_aop_distr, outer_aop_distr)

                        self.generate_loan(inner_loan_max, inner_loan_min, inner_loan_distr)

    #-------                        
    def generate_mass(self, inner_primary_mass_max, inner_primary_mass_min, 
                        inner_secondary_mass_max,inner_secondary_mass_min,outer_mass_min, outer_mass_max,
                        inner_mass_ratio_max, inner_mass_ratio_min,
                        outer_mass_ratio_max, outer_mass_ratio_min,  
                        inner_primary_mass_distr, inner_mass_ratio_distr, outer_mass_ratio_distr):
        if REPORT_TPS:
            print('generate_mass')

        if inner_primary_mass_max == inner_primary_mass_min:
            self.inner_primary_mass = inner_primary_mass_min
        else:
            if inner_primary_mass_distr == 1: #Scalo 1986
                self.inner_primary_mass= new_scalo_mass_distribution(1, inner_primary_mass_max)[0]
                while self.inner_primary_mass < inner_primary_mass_min:
                        self.inner_primary_mass = new_scalo_mass_distribution(1, inner_primary_mass_max)[0]
            elif inner_primary_mass_distr == 2:#Miller & Scale 1979
                self.inner_primary_mass= new_miller_scalo_mass_distribution(1, inner_primary_mass_max)[0]
                while self.inner_primary_mass < inner_primary_mass_min:
                        self.inner_primary_mass = new_miller_scalo_mass_distribution(1, inner_primary_mass_max)[0]
            elif inner_primary_mass_distr == 3: #Salpeter with slope 2.35
                self.inner_primary_mass= new_salpeter_mass_distribution(1, inner_primary_mass_min, inner_primary_mass_max)[0]
            elif inner_primary_mass_distr == 4: # Flat in log space
                self.inner_primary_mass= new_flat_mass_distribution(1, inner_primary_mass_min, inner_primary_mass_max)[0]
            elif inner_primary_mass_distr == 5: # Eggleton 2009, 399, 1471, Salpeter-like with turnover at low masses
                self.inner_primary_mass= eggleton_mass_distr(inner_primary_mass_min, inner_primary_mass_max)
            elif inner_primary_mass_distr == 6: #Salpeter with slope 2.3 -> Kroupa for M>0.5
                self.inner_primary_mass = powerlaw_distr(inner_primary_mass_min, inner_primary_mass_max, -2.3)
            else: #Kroupa 2001
                self.inner_primary_mass = new_kroupa_mass_distribution(1, mass_min = inner_primary_mass_min, mass_max = inner_primary_mass_max)[0]


        if inner_mass_ratio_max == inner_mass_ratio_min:
            inner_mass_ratio = inner_mass_ratio_min
            self.inner_secondary_mass = inner_mass_ratio * self.inner_primary_mass
        else: 
            if inner_mass_ratio_distr == 1:# Kroupa 2001 
                self.inner_secondary_mass = new_kroupa_mass_distribution(1, mass_min=inner_secondary_mass_min, mass_max=inner_secondary_mass_max)[0]
#                self.inner_secondary_mass = new_kroupa_mass_distribution(1, mass_min=inner_secondary_mass_min, mass_max=inner_primary_mass_max)[0]
#                self.inner_secondary_mass = new_kroupa_mass_distribution(1, mass_min=inner_secondary_mass_min, mass_max=inner_primary_mass)[0]
            elif inner_mass_ratio_distr == 2:   # Galicher et al 2016  
                self.inner_secondary_mass = powerlaw_distr( m_min= inner_secondary_mass_min, m_max= inner_secondary_mass_max, slope= -1.31)            
            else: # flat distribution  
               inner_mass_ratio = flat_distr(max(inner_mass_ratio_min, inner_secondary_mass_min/self.inner_primary_mass), min(inner_mass_ratio_max, inner_secondary_mass_max/self.inner_primary_mass))
               self.inner_secondary_mass = inner_mass_ratio * self.inner_primary_mass        
       


        if outer_mass_ratio_max == outer_mass_ratio_min:
            outer_mass_ratio = outer_mass_ratio_min
            self.outer_mass = outer_mass_ratio * (self.inner_primary_mass + self.inner_secondary_mass)
        else: 
            if outer_mass_ratio_distr == 1:# Kroupa 2001 
                self.outer_mass = new_kroupa_mass_distribution(1, mass_min=outer_mass_min, mass_max=outer_mass_max)[0]
            elif outer_mass_ratio_distr == 2:   # Galicher et al 2016  
                self.outer_mass = powerlaw_distr( m_min= outer_mass_min, m_max= outer_mass_max, slope= -1.31)            
            else: # flat distribution
               inner_mass_tot = self.inner_primary_mass + self.inner_secondary_mass
               outer_mass_ratio = flat_distr(max(outer_mass_ratio_min,outer_mass_min/inner_mass_tot), min(outer_mass_ratio_max, outer_mass_max/inner_mass_tot))
               self.outer_mass = outer_mass_ratio * inner_mass_tot



    def generate_semi_and_ecc(self, 
                        inner_semi_max_orig, inner_semi_min_orig, 
                        outer_semi_max_orig, outer_semi_min_orig,
                        inner_semi_distr,  outer_semi_distr,
                        inner_semi_latus_rectum_min, outer_semi_latus_rectum_min, 
                        inner_semi_latus_rectum_max, outer_semi_latus_rectum_max, 
                        inner_ecc_max, inner_ecc_min, 
                        outer_ecc_max, outer_ecc_min,
                        inner_ecc_distr, outer_ecc_distr, 
                        inner_secondary_mass, outer_mass):
        if REPORT_TPS:
            print('generate_semi_and_ecc')
                                
        self.inner_ecc = self.generate_ecc_1d(inner_ecc_max, inner_ecc_min, inner_ecc_distr, inner_secondary_mass)
        self.outer_ecc = self.generate_ecc_1d(outer_ecc_max, outer_ecc_min, outer_ecc_distr, outer_mass)

        inner_semi_min = inner_semi_min_orig
        if inner_semi_latus_rectum_min:
            inner_semi_min = inner_semi_min_orig /(1-self.inner_ecc**2)                       
        outer_semi_min = outer_semi_min_orig 
        if outer_semi_latus_rectum_min:
            outer_semi_min = outer_semi_min_orig /(1-self.outer_ecc**2)                       

        inner_semi_max = inner_semi_max_orig
        if inner_semi_latus_rectum_max:
            inner_semi_max = inner_semi_max_orig /(1-self.inner_ecc**2)                       
        outer_semi_max = outer_semi_max_orig 
        if outer_semi_latus_rectum_max:
            outer_semi_max = outer_semi_max_orig /(1-self.outer_ecc**2)                       

            
                        
        if inner_semi_max == inner_semi_min:
            self.inner_semi = inner_semi_min
        else:
            if inner_semi_distr == 1: #Constant 
                if REPORT_USER_WARNINGS_TPS:
                    print('TPS::generate_semi: unambiguous choice of constant semi-major axis')
                    print('--A_min option to set the value of the semi-major axis in the inner binary')                
                self.inner_semi = inner_semi_min
            elif inner_semi_distr == 2: #Tokovinin Lognormal mu=10^5d, sigma=2.3
                self.inner_semi = 0.|units.RSun
                while (self.inner_semi < inner_semi_min or self.inner_semi > inner_semi_max):
                   self.inner_ecc = self.generate_ecc_1d(inner_ecc_max, inner_ecc_min, inner_ecc_distr, inner_secondary_mass)
                   inner_semi_min = inner_semi_min_orig
                   if inner_semi_latus_rectum_min:
                       inner_semi_min = inner_semi_min_orig /(1-self.inner_ecc**2)                       

                   logP = np.random.normal(5, 2.3, 1)
                   P = (10**logP[0])|units.day
                   self.inner_semi = ((P/2./np.pi)**2 * constants.G* (self.inner_primary_mass + self.inner_secondary_mass))**(1./3.)  
                   if logP < -0.3 or logP > 10:#truncation of Gaussian wings
                        self.inner_semi = 0.|units.RSun
            elif inner_semi_distr == 3: #Lognormal mu=10^3.5d, sigma=2.3
                self.inner_semi = 0.|units.RSun
                while (self.inner_semi < inner_semi_min or self.inner_semi > inner_semi_max):
                   self.inner_ecc = self.generate_ecc_1d(inner_ecc_max, inner_ecc_min, inner_ecc_distr, inner_secondary_mass)
                   inner_semi_min = inner_semi_min_orig
                   if inner_semi_latus_rectum_min:
                       inner_semi_min = inner_semi_min_orig /(1-self.inner_ecc**2)                       
                   inner_semi_max = inner_semi_max_orig
                   if inner_semi_latus_rectum_max:
                       inner_semi_max = inner_semi_max_orig /(1-self.inner_ecc**2)                       

                   logP = np.random.normal(3.5, 2.3, 1)
                   P = (10**logP[0])|units.day
                   self.inner_semi = ((P/2./np.pi)**2 * constants.G* (self.inner_primary_mass + self.inner_secondary_mass))**(1./3.)  
                   if logP < -0.3 or logP > 10:#truncation of Gaussian wings
                        self.inner_semi = 0.|units.RSun
            elif inner_semi_distr == 4: #Rizzuto et al 2013, 436, 1694, Lognormal mu=10^0.95AU, sigma=1.35 
                self.inner_semi = 0.|units.RSun
                while (self.inner_semi < inner_semi_min or self.inner_semi > inner_semi_max):
                   self.inner_ecc = self.generate_ecc_1d(inner_ecc_max, inner_ecc_min, inner_ecc_distr, inner_secondary_mass)
                   inner_semi_min = inner_semi_min_orig
                   if inner_semi_latus_rectum_min:
                       inner_semi_min = inner_semi_min_orig /(1-self.inner_ecc**2)                       
                   inner_semi_max = inner_semi_max_orig
                   if inner_semi_latus_rectum_max:
                       inner_semi_max = inner_semi_max_orig /(1-self.inner_ecc**2)                       

                   logAU = np.random.normal(0.95, 1.35, 1)
                   self.inner_semi = (10**logAU[0])|units.AU
                   if self.inner_semi < 0.5|units.RSun or self.inner_semi > 5e8|units.RSun:#truncation of Gaussian wings
                        self.inner_semi = 0.|units.RSun

            elif inner_semi_distr == 5: #Sana
               self.inner_semi = 0.|units.RSun
               # (logP)^-0.55
               while (self.inner_semi < inner_semi_min or self.inner_semi > inner_semi_max):               
                   self.inner_ecc = self.generate_ecc_1d(inner_ecc_max, inner_ecc_min, inner_ecc_distr, inner_secondary_mass)
                   inner_semi_min = inner_semi_min_orig
                   if inner_semi_latus_rectum_min:
                       inner_semi_min = inner_semi_min_orig /(1-self.inner_ecc**2)                       
                   inner_semi_max = inner_semi_max_orig
                   if inner_semi_latus_rectum_max:
                       inner_semi_max = inner_semi_max_orig /(1-self.inner_ecc**2)                       

                   random_nr = flat_distr(0, 1)
                   logP_min = 0.15
                   logP_max = 8.5
                   c_s = (logP_max**0.45 - logP_min**0.45)
                   logP = (random_nr*c_s +logP_min**0.45)**(1./0.45)
                   P0 = 10**logP|units.day
                   M_inner = self.inner_primary_mass + self.inner_secondary_mass
                   self.inner_semi = ((P0/2./np.pi)**2 * M_inner*constants.G ) ** (1./3.)                

            elif inner_semi_distr == 6:     # flat distr (uniform)
                self.inner_semi = flat_distr( inner_semi_min.value_in(units.RSun), inner_semi_max.value_in(units.RSun))|units.RSun
            elif inner_semi_distr == 7:     # Galicher 2016: powerlaw, slope -0.61
                self.inner_semi = powerlaw_distr( inner_semi_min, inner_semi_max, slope= -0.61)

            else: # log flat distribution
                maximal_semi = min(inner_semi_max, outer_semi_max)
                if inner_semi_min > maximal_semi: #possible for extreme eccentricities
                    return False
                self.inner_semi = log_flat_distr(inner_semi_min.value_in(units.RSun), maximal_semi.value_in(units.RSun))|units.RSun
                        
        
                        
        if outer_semi_max == outer_semi_min:
            self.outer_semi = outer_semi_min
        else:
            if outer_semi_distr == 1: #Constant 
                if REPORT_USER_WARNINGS_TPS:
                    print('TPS::generate_semi: unambiguous choise of constant semi-major axis')
                    print('--a_min option to set the value of the semi-major axis in the outer binary')                
                self.outer_semi = outer_semi_min
            elif outer_semi_distr == 2: #Tokovinin Lognormal mu=10^5.5d, sigma=2.3
                self.outer_semi = 0.|units.RSun
                while (self.outer_semi < outer_semi_min or self.outer_semi > outer_semi_max):
                    self.outer_ecc = self.generate_ecc_1d(outer_ecc_max, outer_ecc_min, outer_ecc_distr, outer_mass)
                    outer_semi_min = outer_semi_min_orig
                    if outer_semi_latus_rectum_min:
                       outer_semi_min = outer_semi_min_orig /(1-self.outer_ecc**2)                       
                    outer_semi_max = outer_semi_max_orig
                    if outer_semi_latus_rectum_max:
                       outer_semi_max = outer_semi_max_orig /(1-self.outer_ecc**2)                       

                    logP_out = np.random.normal(5, 2.3, 1)
                    P_out = (10**logP_out[0])|units.day
                    self.outer_semi = ((P_out/2./np.pi)**2 * constants.G* (self.inner_primary_mass + self.inner_secondary_mass + self.outer_mass))**(1./3.)                    
                    if logP_out < -0.3 or logP_out > 10:#truncation of Gaussian wings
                        self.outer_semi = 0.|units.RSun
                    if logP_out < 3: # no bifurcation
                        self.outer_semi = 0.|units.RSun

            elif outer_semi_distr == 3: #Lognormal mu=10^3.5d, sigma=2.3
                self.outer_semi = 0.|units.RSun
                while (self.outer_semi < outer_semi_min or self.outer_semi > outer_semi_max):
                    self.outer_ecc = self.generate_ecc_1d(outer_ecc_max, outer_ecc_min, outer_ecc_distr, outer_mass)
                    outer_semi_min = outer_semi_min_orig
                    if outer_semi_latus_rectum_min:
                       outer_semi_min = outer_semi_min_orig /(1-self.outer_ecc**2)                       
                    outer_semi_max = outer_semi_max_orig
                    if outer_semi_latus_rectum_max:
                       outer_semi_max = outer_semi_max_orig /(1-self.outer_ecc**2)                       

                    logP_out = np.random.normal(5, 2.3, 1)
                    P_out = (10**logP_out[0])|units.day
                    self.outer_semi = ((P_out/2./np.pi)**2 * constants.G* (self.inner_primary_mass + self.inner_secondary_mass + self.outer_mass))**(1./3.)                    
                    if logP_out < -0.3 or logP_out > 10:#truncation of Gaussian wings
                        self.outer_semi = 0.|units.RSun
                    if logP_out < 3: # no bifurcation
                        self.outer_semi = 0.|units.RSun
                            
            elif outer_semi_distr == 4: #Rizzuto et al 2013, 436, 1694, Lognormal mu=10^0.95AU, sigma=1.35 
                self.outer_semi = 0.|units.RSun
                while (self.outer_semi < outer_semi_min or self.outer_semi > outer_semi_max):
                   self.outer_ecc = self.generate_ecc_1d(outer_ecc_max, outer_ecc_min, outer_ecc_distr, outer_mass)
                   outer_semi_min = outer_semi_min_orig
                   if outer_semi_latus_rectum_min:
                       outer_semi_min = outer_semi_min_orig /(1-self.outer_ecc**2)                       
                   outer_semi_max = outer_semi_max_orig
                   if outer_semi_latus_rectum_max:
                       outer_semi_max = outer_semi_max_orig /(1-self.outer_ecc**2)                       

                   logAU = np.random.normal(0.95, 1.35, 1)
                   self.outer_semi = (10**logAU[0])|units.AU
                   if self.outer_semi < 0.5|units.RSun or self.outer_semi > 5e8|units.RSun:#truncation of Gaussian wings
                        self.outer_semi = 0.|units.RSun

            elif outer_semi_distr == 5: #Sana
               self.outer_semi = 0.|units.RSun
               # (logP)^-0.55
               while (self.outer_semi < outer_semi_min or self.outer_semi > outer_semi_max):
                   self.outer_ecc = self.generate_ecc_1d(outer_ecc_max, outer_ecc_min, outer_ecc_distr, outer_mass)
                   outer_semi_min = outer_semi_min_orig
                   if outer_semi_latus_rectum_min:
                       outer_semi_min = outer_semi_min_orig /(1-self.outer_ecc**2)                       
                   outer_semi_max = outer_semi_max_orig
                   if outer_semi_latus_rectum_max:
                       outer_semi_max = outer_semi_max_orig /(1-self.outer_ecc**2)                       

                   random_nr = flat_distr(0, 1)
                   logP_min = 0.15
                   logP_max = 8.5
                   c_s = (logP_max**0.45 - logP_min**0.45)
                   logP = (random_nr*c_s +logP_min**0.45)**(1./0.45)
                   P0 = 10**logP|units.day
                   mass_tot = self.inner_primary_mass + self.inner_secondary_mass + self.outer_mass
                   self.outer_semi = ((P0/2./np.pi)**2 * mass_tot*constants.G) ** (1./3.)                
                                                         
            elif outer_semi_distr == 6:     # flat distr (uniform)
                self.outer_semi = flat_distr( outer_semi_min.value_in(units.RSun), outer_semi_max.value_in(units.RSun))|units.RSun
            elif outer_semi_distr == 7:     # Galicher 2016: powerlaw, slope -0.61
                self.outer_semi = powerlaw_distr( outer_semi_min, outer_semi_max, slope= -0.61)

            else: # log flat distribution
                if outer_semi_min > outer_semi_max: #possible for extreme eccentricities
                    return False
                self.outer_semi = log_flat_distr(outer_semi_min.value_in(units.RSun), outer_semi_max.value_in(units.RSun))|units.RSun
                       
                     
        if inner_semi_distr == outer_semi_distr and inner_semi_min_orig == outer_semi_min_orig and inner_semi_max == outer_semi_max and self.outer_semi < self.inner_semi and inner_ecc_distr == outer_ecc_distr and inner_ecc_min == outer_ecc_min and inner_ecc_max == outer_ecc_max:
            swap = self.outer_semi
            self.outer_semi = self.inner_semi
            self.inner_semi = swap
            swap = self.outer_ecc
            self.outer_ecc = self.inner_ecc
            self.inner_ecc = swap
            return True
        elif self.outer_semi < self.inner_semi:
            return False
        return True


    def generate_ecc_1d(self, ecc_max, ecc_min, ecc_distr, mass):
        if REPORT_TPS:
            print('generate_ecc_1d')
                        
        if ecc_max == ecc_min:
            return ecc_min
        else:
            if ecc_distr == 1: #Constant 
                if REPORT_USER_WARNINGS_TPS:
                    print('TPS::generate_ecc: unambiguous choise of constant eccentricity')
                    print('--E_min and --e_min option to set the value of the eccentricity')                
                return ecc_min
            elif ecc_distr == 2: #Sana
                return powerlaw_distr(ecc_min+precision, ecc_max, -0.45)
            elif ecc_distr == 3: #flat distribution
                return flat_distr(ecc_min, ecc_max)
            elif ecc_distr == 4: #Powerlaw
                return powerlaw_distr(ecc_min+precision, ecc_max, 0.5)
            elif ecc_distr == 5: # Beta distribution    
                return beta_distr_SSOs(ecc_min, ecc_max, mass)
            else: #Thermal distribution
                 return np.sqrt(np.random.uniform(ecc_min*ecc_min, ecc_max*ecc_max))
                 
      


    def generate_incl(self, incl_max, incl_min, incl_distr):
        if REPORT_TPS:
            print('generate_incl')
        if incl_max == incl_min:
            self.incl = incl_min           
        else:
            if incl_distr == 1: #Constant 
                if REPORT_USER_WARNINGS_TPS:
                    print('TPS::generate_incl: unambiguous choise of constant relative inclination')
                    print('--i_min option to set the value of the relative inclination in the inner triple')                
                self.incl = incl_min
            else: #Circular uniform distribution 
                 self.incl = np.arccos(np.random.uniform(np.cos(incl_min), np.cos(incl_max)))
                 
    def generate_aop(self,
                        inner_aop_max, inner_aop_min, 
                        outer_aop_max, outer_aop_min,
                        inner_aop_distr, outer_aop_distr):
        if REPORT_TPS:
            print('generate_aop')

        if inner_aop_max == inner_aop_min:
            self.inner_aop = inner_aop_min
        else:
            if inner_aop_distr == 1: #Constant 
                if REPORT_USER_WARNINGS_TPS:
                    print('TPS::generate_aop: unambiguous choise of constant argument of pericenter')
                    print('--G_min option to set the value of the argument of pericenter of the inner binary')                
                self.inner_aop = inner_aop_min
            else: #Uniform distribution 
                 self.inner_aop = np.random.uniform(inner_aop_min, inner_aop_max)
                     
                 
        if outer_aop_max == outer_aop_min:
            self.outer_aop = outer_aop_min
        else:
            if outer_aop_distr == 1: #Constant 
                if REPORT_USER_WARNINGS_TPS:
                    print('TPS::generate_aop: unambiguous choise of constant argument of pericenter')
                    print('--g_min option to set the value of the argument of pericenter of the outer binary')                
                self.outer_aop = outer_aop_min
            else: #Uniform distribution 
                 self.outer_aop = np.random.uniform(outer_aop_min, outer_aop_max)
                

    def generate_loan(self, inner_loan_max, inner_loan_min, inner_loan_distr):                                
        if REPORT_TPS:
            print('generate_loan')

        if inner_loan_max == inner_loan_min:
            self.inner_loan = inner_loan_min
        else:
            if inner_loan_distr == 0: #Circular uniform distribution
                self.inner_loan = np.arccos(np.random.uniform(np.cos(inner_loan_min), np.cos(inner_loan_max)))
            else: #Constant
                if REPORT_USER_WARNINGS_TPS:
                    print('TPS::generate_loan: unambiguous choise of constant longitude of ascending nodes')
                    print('--O_min option to set the value of the argument of pericenter of the inner binary')                
                self.inner_loan = inner_loan_min
                    

#-------

# Eggleton 2009, 399, 1471
    def generate_mass_and_semi_eggleton(self, inner_primary_mass_max, inner_primary_mass_min,  
                        inner_secondary_mass_min,outer_mass_min,outer_mass_max,
                        inner_semi_max, inner_semi_min, outer_semi_max, outer_semi_min):
        if REPORT_TPS:
            print('generate_mass_and_semi_eggleton')
        U0_mass = [0., .01, .09, .32, 1., 3.2, 11, 32, np.inf]#solar mass
        U0_l0 = [0.40, 0.40, 0.40, 0.40, 0.50, 0.75, 0.88, 0.94, 0.96]
        U0_l1 = [0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.20, 0.60, 0.80]
        U0_l2 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 0.33, 0.82, 0.90]
        U0_l3 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
              
        Mt = eggleton_mass_distr(min_mass, absolute_max_mass)
        U = np.random.uniform(0, 1)
        f_l0 = interp1d (U0_mass, U0_l0)
        U0 = f_l0(Mt.value_in(units.MSun)) #messy, but otherwise cluster crashes    

        while U >= U0:
            U = np.random.uniform(0, 1)
            f_l0 = interp1d (U0_mass, U0_l0)
            U0 = f_l0(Mt.value_in(units.MSun))     
        
        if U < U0:
            V = np.random.uniform(0, 1)
            P0 = 1.e5 * V**2 / (1-V)**2.5 |units.day
            while P0 > 1e10|units.day:             
                V = np.random.uniform(0, 1)
                P0 = 1.e5 * V**2 / (1-V)**2.5|units.day

            x_p = np.random.uniform(0, 1)
            if P0 > 25|units.day or x_p > 0.25: 
                Q0 = ( (U0-U)/U0 )**0.8
            else:
                Q0 = 0.9+0.09*(U0-U)/U0 
            if Q0 < 0.01:
                Q0 = 0.01                                
                
            M1 = Mt / (1+Q0)
            M2 = M1 * Q0
            f_l1 = interp1d(U0_mass, U0_l1)
            
            U1 = np.random.uniform(0, 1)
            U1_0 = f_l1(M1.value_in(units.MSun))     
            U2 = np.random.uniform(0, 1)
            U2_0 = f_l1(M2.value_in(units.MSun))
            
            #M1 bifurcutas and M2 not
            if U1< U1_0 and U2>=U2_0:
                M_bin = M1  
                U_bin = U1 
                U0_bin = U1_0
                M_comp = M2
            #M2 bifurcutas and M1 not
            elif U1>= U1_0 and U2<U2_0:
                M_bin = M2
                U_bin = U2
                U0_bin = U2_0
                M_comp = M1
            else: #two bifurcations -> higher order multiplicity
#                print(U1, U1_0, U2, U2_0)
#                exit(1)
                 return False
                
            V_bin = np.random.uniform(0, 1)
            P_bin = 0.2 * P0 * 10**(-5*V_bin)

            x_pb = np.random.uniform(0, 1)
            if P_bin > 25|units.day or x_pb > 0.25: 
                Q_bin = ( (U0_bin-U_bin)/U0_bin )**0.8
            else:
                Q_bin = 0.9+0.09*(U0_bin-U_bin)/U0_bin 
            if Q_bin < 0.01:
                Q_bin = 0.01                                

            M1_bin = M_bin / (1+Q_bin)
            M2_bin = M1_bin * Q_bin

            self.inner_primary_mass = M1_bin
            self.inner_secondary_mass = M2_bin
            self.outer_mass = M_comp
            
            self.inner_semi = ((P_bin/2./np.pi)**2 * M_bin*constants.G ) ** (1./3.)
            self.outer_semi = ((P0/2./np.pi)**2 * Mt*constants.G ) ** (1./3.)

            if self.inner_primary_mass < inner_primary_mass_min or self.inner_primary_mass > inner_primary_mass_max:
                return False 
            if self.inner_secondary_mass < inner_secondary_mass_min:
                return False
            if self.outer_mass < outer_mass_min or self.outer_mass > outer_mass_max:
                return False
            if self.inner_semi < inner_semi_min or self.inner_semi > inner_semi_max:
                return False
            if self.outer_semi < outer_semi_min or self.outer_semi > outer_semi_max:
                return False                         

            return True
        else:
            sys.exit('not possible in eggleton distribution')

                 

#-------
        
#-------
    def print_triple(self):
        print('\nTriple - ')
        print('m =', self.inner_primary_mass, self.inner_secondary_mass, self.outer_mass)
        print('a =', self.inner_semi, self.outer_semi)
        print('e =', self.inner_ecc, self.outer_ecc)
        print('i =', self.incl)
        print('g =', self.inner_aop, self.outer_aop)
        print('o =', self.inner_loan, self.inner_loan -np.pi)

    def print_triple_short(self):
        print( self.inner_primary_mass.value_in(units.MSun), self.inner_secondary_mass.value_in(units.MSun), self.outer_mass.value_in(units.MSun), end=" ")
        print( self.inner_semi.value_in(units.RSun), self.outer_semi.value_in(units.RSun), end=" ")
        print( self.inner_ecc, self.outer_ecc, end=" ")
        print( self.incl, end=" ")
        print( self.inner_aop, self.outer_aop, end=" ")
        print( self.inner_loan, self.inner_loan -np.pi, end=" ")
        
#-------

#-------

def evolve_model(inner_primary_mass_max, inner_primary_mass_min,inner_secondary_mass_max, 
                        inner_secondary_mass_min,outer_mass_min,outer_mass_max,
                        inner_mass_ratio_max, inner_mass_ratio_min, 
                        outer_mass_ratio_max, outer_mass_ratio_min, 
                        inner_semi_max, inner_semi_min, outer_semi_max, outer_semi_min, 
                        inner_semi_latus_rectum_min, outer_semi_latus_rectum_min,
                        inner_semi_latus_rectum_max, outer_semi_latus_rectum_max,
                        inner_ecc_max, inner_ecc_min, outer_ecc_max, outer_ecc_min, 
                        incl_max, incl_min,
                        inner_aop_max, inner_aop_min, outer_aop_max, outer_aop_min,
                        inner_loan_max, inner_loan_min, 
                        inner_primary_mass_distr, inner_mass_ratio_distr, outer_mass_ratio_distr,
                        inner_semi_distr,  outer_semi_distr, inner_ecc_distr, outer_ecc_distr, incl_distr,
                        inner_aop_distr, outer_aop_distr, inner_loan_distr,                                                                      
                        metallicity, tend, number, initial_number, seed,
                        stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
                        stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
                        stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
                        stop_at_no_CHE, include_CHE, include_circ,
                        stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
                        stop_at_dynamical_instability, stop_at_semisecular_regime,  
                        stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes,fallback_kick_for_black_holes,
                        stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots):


    i_n = 0
    nr_ids = 0 #number of systems that is dynamically unstable at initialisation
    nr_iss = 0 #number of systems that is in the semisecular regime at initialisation
    nr_imt = 0 #number of systems that has mass transfer at initialisation
    nr_cp = 0 #number of systems with incorrect parameters
    while i_n < number:
        triple_system = Generate_initial_triple(inner_primary_mass_max, inner_primary_mass_min, 
                    inner_secondary_mass_max, inner_secondary_mass_min, outer_mass_min, outer_mass_max,
                    inner_mass_ratio_max, inner_mass_ratio_min, 
                    outer_mass_ratio_max, outer_mass_ratio_min, 
                    inner_semi_max, inner_semi_min, outer_semi_max, outer_semi_min, 
                    inner_semi_latus_rectum_min, outer_semi_latus_rectum_min,
                    inner_semi_latus_rectum_max, outer_semi_latus_rectum_max,
                    inner_ecc_max, inner_ecc_min, outer_ecc_max, outer_ecc_min, 
                    incl_max, incl_min,
                    inner_aop_max, inner_aop_min, outer_aop_max, outer_aop_min,
                    inner_loan_max, inner_loan_min, 
                    inner_primary_mass_distr, inner_mass_ratio_distr, outer_mass_ratio_distr,
                    inner_semi_distr,  outer_semi_distr, inner_ecc_distr, outer_ecc_distr, incl_distr,
                    inner_aop_distr, outer_aop_distr, inner_loan_distr)
                
        if REPORT_TPS:
           triple_system.print_triple()

        if (min_mass > triple_system.inner_primary_mass):
            if REPORT_TPS:
                    print('non-star primary included: ', triple_system.inner_primary_mass)
            continue
        
        if EXCLUDE_SSO: 
            if (min_mass > triple_system.inner_secondary_mass) or (min_mass > triple_system.outer_mass):
                if REPORT_TPS:
                    print('non-star secondary & tertiary included: ', triple_system.inner_secondary_mass, triple_system.outer_mass)
                continue

        number_of_system = initial_number + i_n
        if REPORT_TPS:
            print('number of system = ', number_of_system)

        #do not use main_developer in TPS.py
        #memory of SeBa needs to be cleaned, in particular SeBa time
        #otherwise use evolve_for for particles indivicually -> many calls 
        tr = TRES.main(inner_primary_mass = triple_system.inner_primary_mass, 
                inner_secondary_mass = triple_system.inner_secondary_mass, 
                outer_mass = triple_system.outer_mass, 
                inner_semimajor_axis = triple_system.inner_semi, 
                outer_semimajor_axis = triple_system.outer_semi, 
                inner_eccentricity = triple_system.inner_ecc, 
                outer_eccentricity = triple_system.outer_ecc,
                relative_inclination = triple_system.incl, metallicity = metallicity, tend = tend, number = number_of_system,                     
                stop_at_mass_transfer = stop_at_mass_transfer, stop_at_init_mass_transfer = stop_at_init_mass_transfer,
                stop_at_outer_mass_transfer = stop_at_outer_mass_transfer, 
                stop_at_stable_mass_transfer = stop_at_stable_mass_transfer, 
                stop_at_eccentric_stable_mass_transfer = stop_at_eccentric_stable_mass_transfer, 
                stop_at_unstable_mass_transfer = stop_at_unstable_mass_transfer, 
                stop_at_eccentric_unstable_mass_transfer = stop_at_eccentric_unstable_mass_transfer, 
                stop_at_merger = stop_at_merger, stop_at_disintegrated = stop_at_disintegrated,
                stop_at_inner_collision = stop_at_inner_collision, stop_at_outer_collision = stop_at_outer_collision,
                stop_at_dynamical_instability = stop_at_dynamical_instability, 
                stop_at_semisecular_regime = stop_at_semisecular_regime,  
                stop_at_SN = stop_at_SN, SN_kick_distr = SN_kick_distr, 
                impulse_kick_for_black_holes = impulse_kick_for_black_holes, 
                fallback_kick_for_black_holes = fallback_kick_for_black_holes,
                which_common_envelope = which_common_envelope,
                stop_at_CPU_time = stop_at_CPU_time,
                max_CPU_time = max_CPU_time, file_name = file_name, file_type = file_type, dir_plots = dir_plots, secular_code = secular_code)
    
        if tr.correct_params == False:
            if REPORT_TPS:
                print('Incorrect parameters')            
            nr_cp += 1
        elif tr.semisecular_regime_at_initialisation == True:
            nr_iss +=1
        elif tr.dynamical_instability_at_initialisation == True:
            nr_ids +=1
        elif tr.mass_transfer_at_initialisation == True:
           if tr.has_tertiary_donor():
                nr_imt +=1
           elif include_CHE:
                nr_imt +=1
                # todo reset so that no olof
           else: 
                nr_imt += 1  
                
                if include_circ:
                    i_ecc = 0
                    max_nr_tries_ecc = 10
                    while(i_ecc < max_nr_tries_ecc and tr.mass_transfer_at_initialisation):
                        i_ecc += 1
                        new_ecc = triple_system.generate_ecc_1d(triple_system.inner_ecc, inner_ecc_min, inner_ecc_distr, triple_system.inner_secondary_mass)
                        #resetting semi-major axis creates too many short orbit systems - for now only eccentricity is reset
#                            tr.triple.child2.semimajor_axis *= (1- tr.triple.child2.eccentricity**2)/(1-new_ecc**2) 
#                            triple_system.inner_semi = tr.triple.child2.semimajor_axis
                        triple_system.inner_ecc = new_ecc
                        tr.triple.child2.eccentricity = new_ecc
                        tr.check_RLOF()
                        if not tr.has_donor():
                            i_ecc = max_nr_tries_ecc+1
                            i_n += 1  
                            nr_imt -= 1

                            tr = TRES.main(inner_primary_mass = triple_system.inner_primary_mass, 
                                        inner_secondary_mass = triple_system.inner_secondary_mass, 
                                        outer_mass = triple_system.outer_mass, 
                                        inner_semimajor_axis = triple_system.inner_semi, 
                                        outer_semimajor_axis = triple_system.outer_semi, 
                                        inner_eccentricity = triple_system.inner_ecc, 
                                        outer_eccentricity = triple_system.outer_ecc,  
                                        relative_inclination = triple_system.incl, metallicity = metallicity, tend = tend, number = number_of_system,                     
                                        stop_at_mass_transfer = stop_at_mass_transfer, stop_at_init_mass_transfer = stop_at_init_mass_transfer,
                                        stop_at_outer_mass_transfer = stop_at_outer_mass_transfer, 
                                        stop_at_stable_mass_transfer = stop_at_stable_mass_transfer, 
                                        stop_at_eccentric_stable_mass_transfer = stop_at_eccentric_stable_mass_transfer, 
                                        stop_at_unstable_mass_transfer = stop_at_unstable_mass_transfer, 
                                        stop_at_eccentric_unstable_mass_transfer = stop_at_eccentric_unstable_mass_transfer, 
                                        stop_at_merger = stop_at_merger, stop_at_disintegrated = stop_at_disintegrated,
                                        stop_at_inner_collision = stop_at_inner_collision, stop_at_outer_collision = stop_at_outer_collision,
                                        stop_at_dynamical_instability = stop_at_dynamical_instability, 
                                        stop_at_semisecular_regime = stop_at_semisecular_regime,  
                                        stop_at_SN = stop_at_SN, SN_kick_distr = SN_kick_distr, 
                                        impulse_kick_for_black_holes = impulse_kick_for_black_holes, 
                                        fallback_kick_for_black_holes = fallback_kick_for_black_holes,
                                        which_common_envelope = which_common_envelope,
                                        stop_at_CPU_time = stop_at_CPU_time,
                                        max_CPU_time = max_CPU_time, file_name = file_name, file_type = file_type, dir_plots = dir_plots, secular_code = secular_code)

        else:
            i_n += 1            

        del tr

    if REPORT_TPS:
      print(number, i_n, nr_iss, nr_ids, nr_imt, nr_cp)                              
    secular_code.stop()


def print_distr(inner_primary_mass_max, inner_primary_mass_min, 
                        inner_secondary_mass_max, inner_secondary_mass_min, outer_mass_min, outer_mass_max, 
                        inner_mass_ratio_max, inner_mass_ratio_min, 
                        outer_mass_ratio_max, outer_mass_ratio_min, 
                        inner_semi_max, inner_semi_min, outer_semi_max, outer_semi_min, 
                        inner_semi_latus_rectum_min, outer_semi_latus_rectum_min,
                        inner_semi_latus_rectum_max, outer_semi_latus_rectum_max,
                        inner_ecc_max, inner_ecc_min,  outer_ecc_max, outer_ecc_min, 
                        incl_max, incl_min,
                        inner_aop_max, inner_aop_min, outer_aop_max, outer_aop_min,
                        inner_loan_max, inner_loan_min, 
                        inner_primary_mass_distr, inner_mass_ratio_distr, outer_mass_ratio_distr,
                        inner_semi_distr,  outer_semi_distr, inner_ecc_distr, outer_ecc_distr, incl_distr,
                        inner_aop_distr, outer_aop_distr, inner_loan_distr,                    
                        metallicity, tend, number, initial_number, seed,
                        stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
                        stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
                        stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
                        stop_at_no_CHE, include_CHE, include_circ,
                        stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
                        stop_at_dynamical_instability, stop_at_semisecular_regime,  
                        stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes,fallback_kick_for_black_holes,
                        stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots):

    print('Based on the following distributions:')        
    print('Primary mass: \t\t',                   inner_primary_mass_distr, ' ',lib_inner_primary_mass_distr[inner_primary_mass_distr] )        
    print('Inner mass ratio: \t',               inner_mass_ratio_distr, ' ',lib_inner_mass_ratio_distr[inner_mass_ratio_distr] )        
    print('Outer mass ratio: \t',               outer_mass_ratio_distr, ' ',lib_outer_mass_ratio_distr[outer_mass_ratio_distr] )        
    print('Inner semi-major axis: \t',          inner_semi_distr, ' ',lib_inner_semi_distr[inner_semi_distr] )        
    print('Outer semi-major axis: \t',          outer_semi_distr, ' ',lib_outer_semi_distr[outer_semi_distr] )        
    print('Inner eccentricity: \t',             inner_ecc_distr, ' ',lib_inner_ecc_distr[inner_ecc_distr] )        
    print('Outer eccentricity: \t',             outer_ecc_distr, ' ',lib_outer_ecc_distr[outer_ecc_distr] )        
    print('Inclination: \t\t',                  incl_distr, ' ',lib_incl_distr[incl_distr] )        
    print('Inner aop: \t\t',                    inner_aop_distr, ' ',lib_inner_aop_distr[inner_aop_distr] )        
    print('Outer aop: \t\t',                    outer_aop_distr, ' ',lib_outer_aop_distr[outer_aop_distr] )        
    print('Inner loan: \t\t',                   inner_loan_distr, ' ',lib_inner_loan_distr[inner_loan_distr] )        
    print('Common envelope model: \t',          which_common_envelope, ' ', lib_CE[which_common_envelope])
    print('SN kick distr: \t\t',                SN_kick_distr, ' ', lib_SN_kick_distr[SN_kick_distr])
    print('Metallicity: \t\t',                  '-', ' ', metallicity.value_in(units.none))
    print('\n')


    print('Based on the following assumptions:')
    print('Include CHE: \t\t',                  include_CHE) 
    print('Include circularisation during pre-MS: \t\t',                  include_circ) 
    print('\n')
        
    print('Based on the following stopping conditions:')
    print(stop_at_mass_transfer, '\t Stop at mass transfer')
    print(stop_at_init_mass_transfer, '\t Stop at mass transfer initially')
    print(stop_at_outer_mass_transfer, '\t Stop at outer mass transfer')
    print(stop_at_stable_mass_transfer, '\t Stop at stable mass transfer')
    print(stop_at_eccentric_stable_mass_transfer, '\t Stop at eccentric stable mass transfer')
    print(stop_at_unstable_mass_transfer, '\t Stop at unstable mass transfer')
    print(stop_at_eccentric_unstable_mass_transfer, '\t Stop at eccentric unstable mass transfer')
    print(stop_at_no_CHE, '\t Stop if no chemically homogeneous evolution')
    print(stop_at_merger, '\t Stop at merger')
    print(stop_at_disintegrated, '\t Stop at disintegration')
    print(stop_at_inner_collision, '\t Stop at collision in inner binary')
    print(stop_at_outer_collision, '\t Stop at collision with outer star')
    print(stop_at_dynamical_instability, '\t Stop at dynamical instability')
    print(stop_at_semisecular_regime, '\t Stop at semisecular regime')
    print(stop_at_CPU_time, '\t Stop at maximum CPU time')
    print('\n')


def test_initial_parameters(inner_primary_mass_max, inner_primary_mass_min, 
                        inner_secondary_mass_max, inner_secondary_mass_min, outer_mass_min, outer_mass_max, 
                        inner_mass_ratio_max, inner_mass_ratio_min, 
                        outer_mass_ratio_max, outer_mass_ratio_min, 
                        inner_semi_max, inner_semi_min, outer_semi_max, outer_semi_min, 
                        inner_semi_latus_rectum_min, outer_semi_latus_rectum_min,
                        inner_semi_latus_rectum_max, outer_semi_latus_rectum_max,
                        inner_ecc_max, inner_ecc_min,  outer_ecc_max, outer_ecc_min, 
                        incl_max, incl_min,
                        inner_aop_max, inner_aop_min, outer_aop_max, outer_aop_min,
                        inner_loan_max, inner_loan_min, 
                        inner_primary_mass_distr, inner_mass_ratio_distr, outer_mass_ratio_distr,
                        inner_semi_distr,  outer_semi_distr, inner_ecc_distr, outer_ecc_distr, incl_distr,
                        inner_aop_distr, outer_aop_distr, inner_loan_distr,                    
                        metallicity, tend, number, initial_number, seed,
                        stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
                        stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
                        stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
                        stop_at_no_CHE, include_CHE, include_circ,
                        stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
                        stop_at_dynamical_instability, stop_at_semisecular_regime,  
                        stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes,fallback_kick_for_black_holes,
                        stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots):

    if (inner_primary_mass_min < min_mass) or (inner_primary_mass_max > absolute_max_mass):
        sys.exit("'error: inner primary mass not in allowed range [', min_mass, ',', absolute_max_mass, ']'. min_mass and absolute_max_mass settable in TRES_options.py")
                
    if (inner_secondary_mass_max > absolute_max_mass) :
        sys.exit("'error: inner secondary mass not in allowed range [ < ', absolute_max_mass, ']'. absolute_max_mass settable in TRES_options.py")
    
    if (outer_mass_max > absolute_max_mass) :
        sys.exit("'error: outer mass not in allowed range [ < ', absolute_max_mass, ']'. absolute_max_mass settable in TRES_options.py")
    
    if (inner_secondary_mass_min < absolute_min_mass) :
         sys.exit("'error: inner secondary mass not in allowed range [ >', absolute_min_mass, ']'. absolute_min_mass settable in TRES_options.py")
    if (outer_mass_min < absolute_min_mass)  :
         sys.exit("'error: outer mass not in allowed range [>', absolute_min_mass, ']'. absolute_min_mass settable in TRES_options.py")
    
    if (inner_primary_mass_max < inner_primary_mass_min):
        sys.exit('error: maximum inner primary mass smaller than minimum in primary mass')

    if (inner_secondary_mass_max < inner_secondary_mass_min):
        sys.exit('error: maximum inner secondary mass smaller than minimum in secondary mass')

    if (outer_mass_max < outer_mass_min):
        sys.exit('error: maximum outer mass smaller than minimum in outer mass')


    if (inner_mass_ratio_min < 0.) or (inner_mass_ratio_max > 1.):
        sys.exit('error: inner mass ratio not in allowed range')
        
    if (inner_mass_ratio_max < inner_mass_ratio_min):
        sys.exit('error: maximum inner mass ratio smaller than minimum mass ratio')


    if (outer_mass_ratio_min < 0.) or (outer_mass_ratio_max > 1.):
        sys.exit('error: outer mass ratio not in allowed range')
        
    if (outer_mass_ratio_max < outer_mass_ratio_min):
        sys.exit('error: maximum outer mass ratio smaller than minimum mass ratio')
       

    if (inner_semi_min < 0.5|units.RSun):
        sys.exit('error: inner separation not in allowed range >5 RSun')

    if (outer_semi_min < 0.5|units.RSun):
        sys.exit('error: outer separation not in allowed range >5 RSun')
        
    if (inner_semi_max < inner_semi_min):
        sys.exit('error: maximum inner separation smaller than minimum in separation')

    if (outer_semi_max < outer_semi_min):
        sys.exit('error: maximum outer separation smaller than minimum outer separation')
        
    if (inner_semi_min > outer_semi_max):
        sys.exit('error: maximum outer separation smaller than minimum inner separation - no overlap for inner and outer orbit')
        

    if (inner_ecc_min < 0.) or (inner_ecc_max > 1.):
        sys.exit('error: inner eccentricity not in allowed range [0,1]')

    if (outer_ecc_min < 0.) or (outer_ecc_max > 1.):
        sys.exit('error: outer eccentricity not in allowed range [0,1]')

    if (inner_ecc_max < inner_ecc_min):
        sys.exit('error: maximum inner eccentricity smaller than minimum ecc')

    if (outer_ecc_max < outer_ecc_min):
        sys.exit('error: maximum outer eccentricity smaller than minimum ecc')


    if (incl_min < 0.) or (incl_max > np.pi):
        sys.exit('error: relative inclination not in allowed range [0, pi]')

    if (incl_max < incl_min):
        sys.exit('error: maximum relative inclination smaller than minimum relative inclination')



    if (inner_aop_min < -np.pi) or (inner_aop_max > np.pi):
        sys.exit('error: inner argument of pericenter not in allowed range [-pi,pi]')

    if (outer_aop_min < -np.pi) or (outer_aop_max > np.pi):
        sys.exit('error: outer argument of pericenter not in allowed range [-pi,pi]')

    if (inner_aop_max < inner_aop_min):
        sys.exit('error: maximum inner argument of pericenter smaller than minimum argument of pericenter')

    if (outer_aop_max < outer_aop_min):
        sys.exit('error: maximum outer argument of pericenter smaller than minimum argument of pericenter')


    if (inner_loan_min < -1*np.pi) or (inner_loan_max > np.pi):
        sys.exit('error: inner longitude of ascending node not in allowed range [-pi,pi]')

    if (inner_loan_max < inner_loan_min):
        sys.exit('error: maximum inner longitude of ascending node smaller than minimum argument of pericenter')

    if (number < 1):
        sys.exit('Requested number of systems < 1')

    if (initial_number < 0):
        sys.exit('Initial number of system < 0')


def parse_arguments():
    parser = OptionParser()
    parser.add_option("--M_min", "--M1_min",unit=units.MSun, 
                      dest="inner_primary_mass_min", type="float", default = 1.|units.MSun,
                      help="minimum of inner primary mass [%default]")
    parser.add_option("--M_max", "--M1_max",unit=units.MSun, 
                      dest="inner_primary_mass_max", type="float", default = absolute_max_mass,
                      help="maximum of inner primary mass [%default]")
    parser.add_option("--M_distr", "--M1_distr",dest="inner_primary_mass_distr", type="int", default = 0,
                      help="inner primary mass distribution [Kroupa]")

    parser.add_option("--m_min", "--M2_min",unit=units.MSun, 
                      dest="inner_secondary_mass_min", type="float", default = absolute_min_mass,
                      help="minimum of inner secondary mass [%default]")
    #only used for inner_mass_ratio_distr == 1:# Kroupa 2001    
    parser.add_option("--m_max", "--M2_max",unit=units.MSun, 
                      dest="inner_secondary_mass_max", type="float", default = absolute_max_mass,
                      help="maximum of inner secondary mass [%default]") 
    parser.add_option("--l_min", "--M3_min",unit=units.MSun, 
                      dest="outer_mass_min", type="float", default = absolute_min_mass,
                      help="minimum of outer mass [%default]")
    parser.add_option("--l_max", "--M3_max",unit=units.MSun, 
                      dest="outer_mass_max", type="float", default = absolute_max_mass,
                      help="maximum of outer mass [%default]")
                      
    parser.add_option("--Q_max", "--Qin_max",dest="inner_mass_ratio_max", type="float", default = 1.0,
                      help="maximum of inner mass ratio [%default]")
    parser.add_option("--Q_min", "--Qin_min", dest="inner_mass_ratio_min", type="float", default = 0.,
                      help="minimum of inner mass ratio [%default]")
    parser.add_option("--Q_distr", "--Qin_distr", dest="inner_mass_ratio_distr", type="int", default = 0,
                      help="inner mass ratio distribution [Flat]")

    parser.add_option("--q_max", "--Qout_max", dest="outer_mass_ratio_max", type="float", default = 1.0,
                      help="maximum of outer mass ratio [%default]")
    parser.add_option("--q_min", "--Qout_min", dest="outer_mass_ratio_min", type="float", default = 0.,
                      help="minimum of outer mass ratio [%default]")
    parser.add_option("--q_distr", "--Qout_distr", dest="outer_mass_ratio_distr", type="int", default = 0,
                      help="outer mass ratio distribution [Flat]")

    parser.add_option("--A_min", "--Ain_min", unit=units.RSun,
                      dest="inner_semi_min", type="float", 
                      default = 0.5|units.RSun,
                      help="minimum of inner semi major axis [%default]")
    parser.add_option("--A_max",  "--Ain_max",unit=units.RSun,
                      dest="inner_semi_max", type="float", 
                      default = 5e6|units.RSun,
                      help="maximum of inner semi major axis [%default]")
    parser.add_option("--A_distr",  "--Ain_distr",dest="inner_semi_distr", type="int", default = 0,
                      help="inner semimajor axis distribution [logFlat]")

    parser.add_option("--a_min",  "--Aout_min",unit=units.RSun,
                      dest="outer_semi_min", type="float", 
                      default = 0.5|units.RSun,
                      help="minimum of outer semi major axis [%default]")
    parser.add_option("--a_max",  "--Aout_max", unit=units.RSun,
                      dest="outer_semi_max", type="float", 
                      default = 5e6|units.RSun,
                      help="maximum of outer semi major axis [%default]")
    parser.add_option("--a_distr",  "--Aout_distr", dest="outer_semi_distr", type="int", default = 0,
                      help="outer semimajor axis distribution [logFlat]")

    parser.add_option("--Ar_min",  "--Arin_min", dest="inner_semi_latus_rectum_min", action="store_true", default = False, 
                      help="minimum inner semi latus rectrum  [%default] %unit")
    parser.add_option("--ar_min",  "--Arout_min",dest="outer_semi_latus_rectum_min", action="store_true", default = False, 
                      help="minimum outer semi latus rectrum  [%default] %unit")
    parser.add_option("--Ar_max",  "--Arin_max",dest="inner_semi_latus_rectum_max", action="store_true", default = False, 
                      help="maximum inner semi latus rectrum  [%default] %unit")
    parser.add_option("--ar_max",  "--Arout_max",dest="outer_semi_latus_rectum_max", action="store_true", default = False, 
                      help="maximum outer semi latus rectrum  [%default] %unit")

    parser.add_option("--E_min",  "--Ein_min",
                      dest="inner_ecc_min", type="float", default = 0.,
                      help="minimum of inner eccentricity [%default]")
    parser.add_option("--E_max",   "--Ein_max",
                      dest="inner_ecc_max", type="float", default = 0.9,
                      help="maximum of inner eccentricity [%default]")
    parser.add_option("--E_distr",   "--Ein_distr", dest="inner_ecc_distr", type="int", default = 0,
                      help="inner eccentricity distribution [Thermal]")

    parser.add_option("--e_min",  "--Eout_min",
                      dest="outer_ecc_min", type="float", default = 0.,
                      help="minimum of outer eccentricity [%default]")
    parser.add_option("--e_max",  "--Eout_max",
                      dest="outer_ecc_max", type="float", default = 0.9,
                      help="maximum of outer eccentricity [%default]")
    parser.add_option("--e_distr", "--Eout_distr",dest="outer_ecc_distr", type="int", default = 0,
                      help="outer eccentricity distribution [Thermal]")
                      
                      
    parser.add_option("--i_min, --I_min",
                      dest="incl_min", type="float", default = 0.0,
                      help="minimum of relative inclination [rad] [%default]")
    parser.add_option("--i_max, --I_max",
                      dest="incl_max", type="float", default = np.pi,
                      help="maximum of relative inclination [rad] [%default]")
    parser.add_option("--i_distr, --I_distr", dest="incl_distr", type="int", default = 0,
                      help="relative inclination distribution [Circular uniform]")

                      
    parser.add_option("--G_min",  "--Gin_min",
                      dest="inner_aop_min", type="float", default = -np.pi,
                      help="minimum of inner argument of pericenter [rad] [%default]")
    parser.add_option("--G_max",  "--Gin_max",
                      dest="inner_aop_max", type="float", default = np.pi,
                      help="maximum of inner argument of pericenter [rad] [%default]")
    parser.add_option("--G_distr",   "--Gin_distr",dest="inner_aop_distr", type="int", default = 0,
                      help="inner argument of pericenter distribution [Uniform]")

    parser.add_option("--g_min",  "--Gout_min",
                      dest="outer_aop_min", type="float", default = -np.pi,
                      help="minimum of outer argument of pericenter [rad] [%default]")
    parser.add_option("--g_max", "--Gout_max",
                      dest="outer_aop_max", type="float", default = np.pi,
                      help="maximum of outer argument of pericenter [rad] [%default]")
    parser.add_option("--g_distr",  "--Gout_distr", dest="outer_aop_distr", type="int", default = 0,
                      help="outer argument of pericenter distribution [Uniform]")
                      
    parser.add_option("--O_min",  "--Oin_min",
                      dest="inner_loan_min", type="float", default = -np.pi,
                      help="minimum of inner longitude of ascending node [rad] [%default]")
    parser.add_option("--O_max", "--Oin_max",
                      dest="inner_loan_max", type="float", default = np.pi,
                      help="maximum of inner longitude of ascending node [rad] [%default]")
    parser.add_option("--O_distr",  "--Oin_distr",dest="inner_loan_distr", type="int", default = 1,
                      help="inner longitude of ascending node distribution [Constant]")


    parser.add_option("-z", "-Z",unit=units.none, 
                      dest="metallicity", type="float", default = 0.02|units.none,
                      help="metallicity [%default] %unit")                     
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 13500|units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("-n", dest="number", type="int", default = 10,
                      help="total number of systems to be simulated [%default]")
    parser.add_option("-N", dest="initial_number", type="int", default = 0,
                      help="number ID of first system [%default]")
    parser.add_option("-s", unit=units.none, 
                      dest="seed", type="float", default = 0.|units.none,
                      help="seed [%default] %unit")
#    int actual_seed = srandinter(input_seed);

    parser.add_option("--no_stop_at_mass_transfer", dest="stop_at_mass_transfer", action="store_false", default = True,
                      help="stop at mass transfer [%default] %unit")
    parser.add_option("--no_stop_at_init_mass_transfer", dest="stop_at_init_mass_transfer", action="store_false", default = True,
                      help="stop if initially mass transfer[%default] %unit")
    parser.add_option("--no_stop_at_outer_mass_transfer", dest="stop_at_outer_mass_transfer", action="store_false", default = True,
                      help="stop at outer mass transfer [%default] %unit")

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
    parser.add_option("--CE", dest="which_common_envelope",  type="int", default = 2,
                      help="which common envelope modeling [%default]")                      

    parser.add_option("--stop_at_no_CHE",dest="stop_at_no_CHE", action="store_true", default = False,
                      help="stop if no chemically homogeneous evolution [%default] %unit") 
    parser.add_option("--include_CHE", dest="include_CHE", 
                    action="store_true", default = False, help="include chemically homogeneous evolution in the stellar evolution [%default] %unit")
    parser.add_option("--include_circularisation_during_preMS", dest="include_circ", 
                    action="store_true", default = False, help="include circularisation during pre-MS [%default] %unit")

    parser.add_option("--no_stop_at_merger", dest="stop_at_merger", action="store_false", default = True, 
                      help="stop at merger [%default] %unit")
    parser.add_option("--no_stop_at_disintegrated", dest="stop_at_disintegrated", action="store_false", default = True,
                      help="stop at disintegrated [%default] %unit")
    parser.add_option("--no_stop_at_inner_collision", dest="stop_at_inner_collision", action="store_false",default = True,
                      help="stop at collision in inner binary[%default] %unit")
    parser.add_option("--no_stop_at_outer_collision", dest="stop_at_outer_collision", action="store_false",default = True,
                      help="stop at collision in outer binary[%default] %unit")
    parser.add_option("--no_stop_at_dynamical_instability", dest="stop_at_dynamical_instability", action="store_false", default = True,
                      help="stop at dynamical instability [%default] %unit")
    parser.add_option("--stop_at_semisecular_regime", dest="stop_at_semisecular_regime", action="store_true", default = False,
                      help="stop at semisecular regime [%default] %unit")

    parser.add_option("--stop_at_SN", dest="stop_at_SN", action="store_true", default = False,
                      help="stop at supernova [%default] %unit")
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
                      
                      
    parser.add_option("-f", dest="file_name", type ="string", default = "TRES.hdf",#"TRES.txt"
                      help="file name[%default]")
    parser.add_option("-F", dest="file_type", type ="string", default = "hdf5",#"txt"
                      help="file type[%default]")
    parser.add_option("--dir_plots", dest="dir_plots", type ="string", default = "",#"txt"
                      help="directory for plots for debugging mode [%default]")



    options, args = parser.parse_args()
    return options.__dict__



if __name__ == '__main__':
    options = parse_arguments()
    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    test_initial_parameters(**options)
    print_distr(**options)
    evolve_model(**options)
    
#    stellar_code.stop()
#    secular_code.stop()
    
    print('\nYou have used the TRES triple evolution code. Literature reference:')
    print('** Toonen, Hamers & Portegies Zwart 2016, ComAC, 3, 6T:')
    print('... "The evolution of hierarchical triple star-systems" ')


