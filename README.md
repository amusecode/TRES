# TRES
Triple evolution Simulation package

### Description
TRES is a numerical framework for simulating hierarchical triple systems. 
Mass transfer from one star to another and the consequential effect to the orbital dynamics is realized via heuristic recipes.
These recipes are combined with  three-body  dynamics and stellar evolution inluding their mutual influences. 

This document contains the following parts:

[Compilation](#Compilation)

[Simple examples](#Simple-examples-of-runs)

[Understanding the SeBa output](#Understanding-the-SeBa-output)

[References](#References)



## Compilation

TRES can be compiled with its Makefile as following:
in the home directory: make 

## Simple examples of runs

To evolve a single system with the parameters:
primary mass M=1.2 Solar mass, 
secondary mass m=0.5 Solar mass, 
tertiary mass m=0.6 Solar mass, 
inner and outer eccentricity e=0.1 & 0.5, 
inner and outer orbital separation a=200 & 20000 Solar radii, 
metallicity z=0.001, 
and time T=10 Myrs,  you need to run:

python TRES.py -M 1.2 -m 0.5 -l 0.6 -E 0.1 -e 0.5 -A 200 -a 20000 -z 0.001 -T 10 
(assuming AMUSE is loaded in your python)

### Input parameters 

The full list of possible input parameters is:

        parameter                               unit / default                  
-M      inner_primary_mass                      in Solar mass
-m      inner secondary mass                    in Solar mass 
-l      outer mass                              in Solar mass
-A      inner semi major axis                   in Solar radius
-a      outer semi major axis                   in Solar radius
-E      inner eccentricity 
-e      outer eccentricity 
-i, -I  relative inclination                    in rad
-I  
-G      inner argument of pericenter            in rad
-g      outer argument of pericenter            in rad
-O      inner longitude of ascending node       in rad
    (outer longitude of ascending nodes = inner - pi)               
-z      metallicity                             default = 0.02 (Solar)
-t, -T  end time                                Myr
-N, -n  integer number asigned to the triple    default = 0


-f      name of output file                     default = TRES
-F      type of output file (hdf5/txt)          default = hdf5


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



### Multiple systems with specified parameters

If you need to follow the triple evolution for multiple systems with parameters which are already specified you can start TRES multiple times, e.g.
```
python TRES.py -M 1.2 -m 0.5 -l 0.6 -E 0.1 -e 0.5 -A 200 -a 20000 -z 0.001 -T 10 
python TRES.py -M 1.5 -m 1 -l 0.6 -E 0.1 -e 0.5 -A 50 -a 20000 -z 0.001 -T 10 
python TRES.py -M 1.5 -m 1 -l 0.05 -E 0.1 -e 0.5 -A 50 -a 20000 -z 0.02 -T 10 

```
This is probably not handy for more than 5 systems. Although this can be added in e.g. a shell or Python script.


### Random population
A random population can be generated with TPS.py. 

Detailed instructions  TBD. 


## Understanding the TRES output

Normally TRES adds the evolution history of individual triples in the TRES.hdf file. Every line represents a moment in the evolution of the triple when something interesting happened, for example one of the star transitions from the main-sequence to the hertzsprung gap, or mass transfer starts or stops. The meaning of the columns is defined below. 

TBD

## References

See the following publications: [Toonen etal 2016](https://ui.adsabs.harvard.edu/abs/2016ComAC...3....6T/abstract) for more details.




