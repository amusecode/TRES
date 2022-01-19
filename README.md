# TRES
TRiple Evolution Simulation 

### Description
TRES is a numerical framework for simulating hierarchical triple systems with stellar and planetary components. 
Mass transfer from one star to another and the consequential effect to the orbital dynamics is realized via heuristic recipes.
These recipes are combined with  three-body  dynamics and stellar evolution inluding their mutual influences. 

TRES includes the effects of common-envelope evolution, circularized stable mass transfer, tides, gravitational wave emission and up-to-date stellar evolution through SeBa. Other stellar evolution codes such as SSE can also be used. Coming soon: TRES with MESA, transition to N-body calculations (including stellar evolution and dissipative processes) when the system's evolution is not secular anymore. 

This document contains the following parts:

[Compilation](#Compilation)

[Simple examples](#Simple-examples-of-runs)

[Understanding the TRES output](#Understanding-the-TRES-output)

[Reducing the TRES output](#Reducing-the-TRES-output)

[TRES development team](#TRES-development-team)

[References](#References)



## Compilation

TRES makes use of the Astrophysical Multipurpose Software Environment (AMUSE) See https://amusecode.github.io/ for how to install AMUSE. 
Note that for standard TRES usage, the only necessary community code to install is SeBa. 

After compiling AMUSE, TRES needs to be compiled by means of it Makefile as following:

```

cd seculartriple_TPS
make 

```

## Simple examples of runs

To evolve a single system with the parameters:
primary mass M=1.2 Solar mass, 
secondary mass m=0.5 Solar mass, 
tertiary mass m=0.6 Solar mass, 
inner and outer eccentricity e=0.1 & 0.5, 
inner and outer orbital separation a=200 & 20000 Solar radii, 
metallicity z=0.001, 
and time T=10 Myrs,  you need to run:

```

python TRES.py -M 1.2 -m 0.5 -l 0.6 -E 0.1 -e 0.5 -A 200 -a 20000 -z 0.001 -T 10 

```
assuming AMUSE is loaded in your python. 

### Input parameters 

The full list of possible input parameters is:

```
        parameter                               unit / default
-M      inner_primary_mass                      in Solar mass
-m      inner secondary mass                    in Solar mass 
-l      outer mass                              in Solar mass
-A      inner semi major axis                   in Solar radius
-a      outer semi major axis                   in Solar radius
-E      inner eccentricity 
-e      outer eccentricity 
-i, -I  relative inclination                    in rad  
-G      inner argument of pericenter            in rad
-g      outer argument of pericenter            in rad
-O      inner longitude of ascending node       in rad
    (outer longitude of ascending nodes = inner - pi)               
-z      metallicity                             default = 0.02 (Solar)
-t, -T  end time                                in Myr
-N, -n  integer number asigned to the triple    default = 0


-f      name of output file                     default = TRES
-F      type of output file (hdf5/txt)          default = hdf5
--dir_plots   directory for plots for debugging default = "" (current directory)
        mode  (aka REPORT_DEBUG == True)
--max_CPU_time   maximum CPU time allowed (only works in combination with "stop_at_CPU_time")    
                                                default = 7200 (seconds)

--SN_kick_distr   supernova kick distribution   default = 10
        options:
        0:  No kick 
        1:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
        2:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  scaled down for bh by mass
        3:  Arzoumanian ea 2002, 568, 289
        4:  Arzoumanian ea 2002, 568, 289 scaled down for bh by mass
        5:  Hansen & Phinney 1997, 291, 569
        6:  Hansen & Phinney 1997, 291, 569 scaled down for bh by mass
        7:  Paczynski 1990, 348, 485
        8:  Paczynski 1990, 348, 485 scaled down for bh by mass
        9:  Verbunt, Igoshev & Cator, 2017, 608, 57
        10:  Verbunt, Igoshev & Cator, 2017, 608, 57 scaled down for bh by mass 

```

Additionally, there is a list of stopping conditions that determines whether the simulation of a system should stop at a certain evolutionary phase. 
By default, these stopping conditions are set to True, which means they are in effect. However, the four specific mass transfer cases (stable, unstable, eccentric stable & eccentric unstable) are set to False by default. Once "--no_stop_at_mass_transfer" is set to False, it is possible to set the specific mass transfer cases to True.

```

action items                                    add these to:
--no_stop_at_mass_transfer                      avoid stopping the simulation at the onset of mass transfer 
--no_stop_at_init_mass_transfer                 avoid stopping the simulation if there is mass transfer initially
--no_stop_at_outer_mass_transfer                avoid stopping the simulation when tertiary initiates mass transfer 
                                                methodology is as of yet non-existent
--stop_at_stable_mass_transfer                  avoid stopping the simulation at the onset of stable mass transfer in the inner binary
--stop_at_unstable_mass_transfer                avoid stopping the simulation at the onset of unstable mass transfer in the inner binary (leading to common-envelope evolution)
--stop_at_eccentric_stable_mass_transfer        avoid stopping the simulation at the onset of stable mass transfer in the inner binary if the orbit is eccentric
--stop_at_eccentric_unstable_mass_transfer      avoid stopping the simulation at the onset of unstable mass transfer in the inner binary if the orbit is eccentric

--no_stop_at_merger                             avoid stopping the simulation after a merger
--no_stop_at_inner_collision                    avoid stopping the simulation after a collision in the inner binary
--no_stop_at_outer_collision                    avoid stopping the simulation after a collision involving the outer star
--no_stop_at_disintegrated                      avoid stopping after the system disintegrated into seperate systems
--stop_at_semisecular_regime                    to stop the simulation if the sytem is in the semi secular regime
--stop_at_SN                                    to stop the simulation when a supernova occurs
--stop_at_CPU_time                              to stop the simulation when the computational time exceeds a given value

```


### Multiple systems with specified parameters

If you need to follow the triple evolution for multiple systems with parameters which are already specified you can start TRES multiple times, e.g.
```

python TRES.py -M 1.2 -m 0.5 -l 0.6 -E 0.1 -e 0.5 -A 200 -a 20000 -z 0.001 -T 10 
python TRES.py -M 1.5 -m 1 -l 0.6 -E 0.1 -e 0.5 -A 50 -a 20000 -z 0.001 -T 10 
python TRES.py -M 1.5 -m 1 -l 0.05 -E 0.1 -e 0.5 -A 50 -a 20000 -z 0.02 -T 10 

```
This is probably not handy for more than 5 systems. Although this can be added in e.g. a shell or Python script.


### Random population
A random population can be generated with TPS.py with a Monte Carlo based approach, e.g. 
```

python TPS.py -n 5 
python TPS.py -n 10 --M_max 5 --M_min 4  --M_distr 0 --A_max 2000 --A_min 200 --A_distr 2

```

The full list of options is [default]:
```

--M_max    upper limit for the inner primary mass [100 Msun]
--M_min    lower limit for the inner primary mass [0.1 Msun]
--M_distr  mass function option: 
        0: "Kroupa", #default
        1: "Scalo",
        2: "Miller & Scalo",
        3: "Salpeter",
        4: "Logarithmically flat",
        5: "Eggleton",
        6: "Kroupa for massive stars M>0.5 powerlaw with exp=-2.3",
--Q_max    upper limit for the inner mass ratio [1.]
--Q_min    lower limit for the inner mass ratio [0.]
--Q_distr  inner mass ratio option: 
        0: "Uniform distribution", #default
        1: "Kroupa IMF",
--q_max    upper limit for the outer mass ratio [1.]
--q_min    lower limit for the mass of the outer star [0.]
--q_distr  outer mass ratio option: 
        0: "Uniform distribution", #default
        1: "Kroupa IMF",
--A_max    upper limit for the inner semi-major axis [5e6 RSun]
--A_min    lower limit for the inner semi-major axis [5]
--A_distr  inner semi-major axcis option: 
        0: "Log Uniform distribution", #default
        1: "Constant semi-major axis",
        2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
        3: "Lognormal mu = 10^3.5d, sigma = 2.3",
        4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
        5: "Sana et al. 2012",
--a_max    upper limit for the outer semi-major axis [5e6 RSun]
--a_min    lower limit for the outer semi-major axis [5 RSun]
--a_distr  outer semi-major axis option: 
        0: "Log Uniform distribution", #default
        1: "Constant semi-major axis",
        2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
        3: "Lognormal mu = 10^3.5d, sigma = 2.3",
        4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
        5: "Sana et al. 2012",
--E_max    upper limit for the inner eccentricity [1.]
--E_min    lower limit for the inner eccentricity [0.]
--E_distr  inner eccentricity option: 
        0: "Thermal", #default
        1: "Constant eccentricity",
        2: "Sana et al. 2012 e^-0.45", #-> close binaries
        3: "Flat distribution",
        4: "Powerlaw e^0.5",                                   
--e_max    upper limit for the outer eccentricity [1.]
--e_min    lower limit for the outer eccentricity [0.]
--e_distr  outer eccentricity option: 
        0: "Thermal", #default
        1: "Constant eccentricity",
        2: "Sana et al. 2012 e^-0.45", #-> close binaries
        3: "Flat distribution",
        4: "Powerlaw e^0.5",
--i_max    upper limit for the relative inclination [pi]
--i_min    lower limit for the relative inclination [0]
--i_distr  relative inclination option: 
        0: "Circular uniform distribution", #default
        1: "Constant inclination",
--G_max    upper limit for the inner argument of pericenter [pi]
--G_min    lower limit for the inner argument of pericenter [-pi]
--G_distr  inner argument of pericenter option: r
        0: "Uniform distribution", #default
        1: "Constant argument of pericenter",    
--g_max    upper limit for the outer argument of pericenter [pi]
--g_min    lower limit for the outer argument of pericenter [-pi]
--g_distr  outer argument of pericenter option: 
        0: "Uniform distribution", #default
        1: "Constant argument of pericenter",     
--O_max    upper limit for the inner longitude of ascending node [pi]
--O_min    lower limit for the inner longitude of ascending node [-pi]
--O_distr  inner longitude of ascending node option: 
        0: "Circular niform distribution", 
        1: "Constant longitude of ascending nodes", #default
        (outer longitude of ascending nodes = inner - pi)             
-T or -t   binary end time. [13500 Myr]
-z         metallicity of stars  [0.02 (Solar)] 
-n         number of triples to be simulated  [1]
-N         number of initial triple  [0]

--SN_kick_distr   supernova kick distribution   default = 10
        options:
        0:  No kick 
        1:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
        2:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  scaled down for bh
        3:  Arzoumanian ea 2002, 568, 289
        4:  Arzoumanian ea 2002, 568, 289 scaled down for bh
        5:  Hansen & Phinney 1997, 291, 569
        6:  Hansen & Phinney 1997, 291, 569 scaled down for bh
        7:  Paczynski 1990, 348, 485
        8:  Paczynski 1990, 348, 485 scaled down for bh
        9:  Verbunt, Igoshev & Cator, 2017, 608, 57
        10:  Verbunt, Igoshev & Cator, 2017, 608, 57 scaled down for bh #default

action items                    add these to:
--no_stop_at_merger             avoid stopping the simulation after a merger
--no_stop_at_inner_collision    avoid stopping the simulation after a collision in the inner binary
--no_stop_at_outer_collision    avoid stopping the simulation after a collision involving the outer star
--no_stop_at_disintegrated      avoid stopping after the system disintegrated into seperate systems
--no_stop_at_mass_transfer      avoid stopping the simulation at the onset of mass transfer 
--stop_at_semisecular_regime    to stop the simulation if the sytem is in the semi secular regime
--stop_at_SN                    to stop the simulation when a supernova occurs
--stop_at_CPU_time              to stop the simulation when the computational time exceeds a given value

```


## Understanding the TRES output

Normally TRES adds the evolution history of individual triples in the TRES.hdf file. Every snapshot represents a moment in the evolution of the triple when something interesting happened, for example one of the star transitions from the main-sequence to the hertzsprung gap, or mass transfer starts or stops. 


## Reducing the TRES output

The python script rdc_TRES.py reduce the TRES hdf output to a txt file keeping only selected parameters. These can be adjusted to your liking in the function rdc(). Currently there are 6 lines for every snapshot. The columns represent:

General information on the system: 
```
Line 1: snapshot number, time, triple number, relative_inclination, dynamical_instability, kozai_type, error_flag_secular
```
Orbital information (inner binary | outer binary) :
```
Line 2: 'bs:', binary type, stability, semimajoraxis, eccentricity, argument_of_pericenter, longitude_of_ascending_node 
        | binary type, stability, semimajoraxis, eccentricity, argument_of_pericenter, longitude_of_ascending_node 
```
Stellar information (primary | secondary | tertiary)
```
Line 3: 'st:', is_donor, stellar_type, mass, spin_angular_frequency, radius, core mass
        | is_donor, stellar_type, mass, spin_angular_frequency, radius, core mass
        | is_donor, stellar_type, mass, spin_angular_frequency, radius, core mass
```        

The stellar types in TRES follow the standard terminology of AMUSE:
```
0 deeply or fully convective low mass MS star
1 Main Sequence star
2 Hertzsprung Gap
3 First Giant Branch
4 Core Helium Burning
5 First Asymptotic Giant Branch
6 Second Asymptotic Giant Branch
7 Main Sequence Naked Helium star
8 Hertzsprung Gap Naked Helium star
9 Giant Branch Naked Helium star
10 Helium White Dwarf
11 Carbon/Oxygen White Dwarf
12 Oxygen/Neon White Dwarf
13 Neutron Star
14 Black Hole
15 Massless Supernova
16 Unknown stellar type
17 Pre-main-sequence Star
18 Planet
19 Brown Dwarf
```


## References

See the following publication: [Toonen et al 2016](https://ui.adsabs.harvard.edu/abs/2016ComAC...3....6T/abstract) for more details.




