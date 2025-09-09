# Making Mayonaise 

## Section 1: Introduction 
 
The goal of this project is to study fluids often encountered in the food processing industry. 

Dressings, sauces and mayonaises are typically be described as [emulsions](https://en.wikipedia.org/wiki/Emulsion) obtained by mixing oil and water. Industrial partners in this project are eager to further develop modeling and simulations techniques to obtain a better grip on the batch production and long term storage of emulsions. This grip is obtained by a better understanding of the physical properties of the mixture in terms of the properties of the oil and water components. 

We will consider [mayonnaise](https://en.wikipedia.org/wiki/Mayonnaise) as an emulsions obtained by mixing oil droplets in water. These emulsion can be described as [linear visco-elastic materials](https://en.wikipedia.org/wiki/Viscoelasticity). An important physical property in the study of emulsions is the [dynamic modulus](https://en.wikipedia.org/wiki/Dynamic_modulus). The modulus is defined as the ratio of stress to strain under vibratory conditions. Given that strain is dimensionless, the dynamic modulus has the same units as stress. The dynamic modulus is complex-valued quantify. The real part is called the storage modulus and represents the elastic portion of the material. The imaginary part is called the loss modulus and represent the viscous part of the material. Our interest is understanding how the dynamic modulus of mayonaise dependent on its ingredients and how it evolves over the shell live time. 

## Section 2: Oil Mixture Fraction and Visco-Elastic Properties of Emulsions 

The dynamic modulus of mayonaise depends on (among many others factors) the types and amount of oil used. The dynamic modulus can be both measured experimentally and modeled mathematically.    

### Measuring the Strength of Mayonaise 

The dynamic modulus of mayonaise can be measured using [dynamic mechanical analysis](https://en.wikipedia.org/wiki/Dynamic_mechanical_analysis) (general term) and [rheometry](https://en.wikipedia.org/wiki/Rheometry) (more specific terms). 

1. picture of weak and strong mayonaise;
2. provide figures for $G'(\omega)$ and $G''(\omega)$: explain that sample breaks down at sufficiently large frequencies as sample no longer allowed to restore from exitation, multiple relaxation times).

### Modeling the Strength of Mayonaise  

1. single relaxation time [Kelvin-Voigt model](https://en.wikipedia.org/wiki/Kelvin–Voigt_material) (parallel spring-dashpot)
2. single relaxation time [Maxwell model](https://en.wikipedia.org/wiki/Maxwell_model) (series spring-dashpot)  

### Mason-Scheffold-2014

We suggest to start the project the studying the 2014 paper by Mason and Sheffold. In a particular frequency range, the storage modulus can be consider to be constant. This value is refered to as the plateau elastic modulus. 

We suggest to rephrase 
1. the model (the minimization of sum of total free energy, i.e., the sum of entropic and interfacial energy);
2. the solution approach (finding first order critical points) and;
3. the results obtained (small and large plateau elastic modulus at small and large oil mixture fraction);

see seperate notebook [mason-sheffold-2014.ipynb](mason-sheffold-2014.ipynb).

Assume a mono-disperse oil-in-water emulsion subject to a shear strain with value $]\gamma$. 
1. $N$: number of droplets; 
1. $a$: oil droplet radius, in the order of nano-meters (nano-emulsions) or micro-meter (micro-emulsions);
2. $V_{drop} = 4 \pi \, a^3 / 3$: volume of droplet;  
3. $\phi$: droplet mixture fraction for oil-in-water emulsion;
4. $\phi_c$: threshold value of $\phi$ for jamming. In the range $0.6 \leq \phi_c \leq 0.7$. Paper uses $\phi_c = 0.64$; 
5. $\phi'_c > \phi_c$: new and larger value for the threshold value of jamming;
6. $\phi_d = \phi'_c - \phi_c$: deformation mixture fraction. Precise value is given by imposing equilibrium in energy contributions; 
7. $\sigma$: oil in water [interfacial tension](https://www.youtube.com/watch?v=2BBHl8Zvs2U&t=40s) of the droplets, in the order of miliNewtons/meter. $\sigma/a$ is the Laplace pressure scale. Q: is this value independent of the strain applied to the system?;
8. $G'_p(\phi)$: linear (or plateau) shear [elastic modulus](https://en.wikipedia.org/wiki/Elastic_modulus) (as opposed to Youngs or bulk modulus). Also referred to as effective spring constant. Clearly independent of frequency of shearing force. Assumption is that the emulsion is in the state in which the shear tension responds linear to the shear rate. In this state the definition as second order derivative applies;
9. $\Pi(\phi)$: [osmotic pressure](https://en.wikipedia.org/wiki/Osmotic_pressure)
10. $F(\phi,\gamma)/N$: free energy per particle

### Scaling laws based on [Surface Evolver](https://kenbrakke.com/evolver/evolver.html) 
1. scaling laws for $G'_p(\phi) \sim \sigma/a \, \phi \, (\phi-\phi_c)$;
2. scaling laws for osmotic pressure $\Pi(\phi) \sim \sigma/a \, \phi^2 \, (\phi-\phi_c)$

### Scaling laws based on jamming simulations  
1. scaling laws for $G'_p(\phi) \sim \sigma/a \, [ \phi^{8/3} \, (\phi-\phi_c)^{0.82} + \phi^{5/3} \, (\phi-\phi_c)^{1.82} ]$;

### Free Energy
1. entropic free energy: $ F_{entropic}/N = -3 \, k_B \, T \, \log[ \phi_c + \phi_d - \phi - \alpha \, \gamma^2] $ 
3. interfacial free energy: $ F_{interfacial}/N = 4 \, \pi \, \xi \, a \, \phi_d^2 $  

### Energy Minimization 
1. first order critical point determines the optinmal deformation volume fraction $\phi^*_d$

### Computation of osmotic pressure and plateau shear modulus 

## Section 3: Complement Understanding of the Problem 

We suggest to further improve the understanding of concepts introduced by Mason-Sheffold-2014 by complementing with material from other sources. We suggest consulting e.g. the 1999 book by Russel, Saville and Schowalter entitled “Colloidal Dispersions”. 

## Section 4: Oil Droplet Size and Visco-Elastic Properties

The framework introduced by Mason-Sheffold-2014 is limited to mono-disperse emulsions in which all oil droplets have the same diameter. Motivated by practical applications, we wish to extend the Mason-Sheffold-2014 model to poly-disperse emulsions. 

### Extension of the Energy Minimization Model 

A first alternative to is to first generalize the expression for the entropic and interfacial free energy to a weighted sum of energy contributions valid for a mixture of oil droplet sizes. Subsequently, one needs to determine the critical mixture fraction that determines the equilibrium of the configuration. Finally, one can determine the plateau shear modulus by the second derivative of the total free energy with respect to the strain rate. 

### Particle Based Models

A second alternative is to replace the energy minimization framework by a mechanical model of a number of oil drops interconnected by springs, dashpots and dampers. One wishes to study the stress response of these systems given as oscillatory shear. Possibly one can start by considering interconnections (networks or graphs) of Kelvin-Voight and Maxwell models. Possibly one can borrow ideas from molecular dynamics or other particle-based simulations methods.

<b>Example of a Particle Based Model</b> [computation of the viscosity](https://docs.lammps.org/Howto_viscosity.html) by ensemble averaging of the auto-correlation of the stress/pressure tensor using [LAMMPS](https://docs.lammps.org/Manual.html);

### Continuum Models 

A third alternative is to study computational fluid dynamics models of the rheometer that determines the plateau shear modulus. One thus needs to solve a rotating lid-driven cavity Stokes flow model for a shear-thinning non-Newtonian fluid.  

## Section 5: Rotating-Lid-Driven Cavity Stokes Flow for Shear-Thinning Non-Newtonian Fluid

### Section 1.5: Pre-Processor 

The preprocessor consist of the following two components: 

1. geometry modeling tool;
2. mesh generation tool;  

The first component is the geometry modeler.  This modeler is either implements a form of [constructive solid geometry](https://en.wikipedia.org/wiki/Constructive_solid_geometry) or allows to import a geometry definition from file in standard definition. 

The second component is the mesh generation tool. Here we use [GMSH](https://gmsh.info). 

Here, the preprocessor should generate a mesh on the cylindrical cavity of the rheometer in which the sample resides. 

<b>Coding</b>: see seperate notebook [gmsh_cavity.ipynb](gmsh_cavity.ipynb).

## Section 2.5: Computational Kernel 

The computational kernel assumes the responses of the sample to the excitation of the rheometer to be modeled in terms of physical laws of conservation (conservastion of mass and momentum) and constituive laws of the material (stress-strain relation). Given these models, the computational kernel consists of the following three four components: 

1. read the mesh;
2. discretize the conservation and constitutive equations in space and time;
3. solve the discrete model in space and time;
4. write snapshots for spatial distribution of velocity, strain and stress to file for visualization;
   
Simulations reported here were performed using the [Ferrite](https://ferrite-fem.github.io/Ferrite.jl/stable/) finite element software.

### Section 3.5: Scalar Diffusion 

A elementary proof of concept is to solve a scalar diffusion eqiation on the volume (computational domain) of a cylindrical cavity. 

<b>Coding</b>: see seperate notebook [ferrite_diffusion_cavity.ipynb](ferrite_diffusion_cavity.ipynb).

### Section 4.5: Stationary Stokes Flow of a Newtonian Fluid 

Here we consider the rotation of the lid covering the cylindrical cavity to be modeled by a overly simplified [Stokes flow](https://en.wikipedia.org/wiki/Stokes_flow) model for a [Newtonian fluid](https://en.wikipedia.org/wiki/Newtonian_fluid).  

<b>Coding</b>: see seperate notebook [ferrite_stokes_cavity.ipynb](ferrite_stokes_cavity.ipynb). 

### Section 5.5: Transient Navier-Stokes Flow of a Newtonian Fluid

Here we extend flow Stokes to [Navier-Stokes flow](https://en.wikipedia.org/wiki/Navier–Stokes_equations), thus including non-linear convective terms. We also include time-dependent terms;

<b>Coding</b>: see seperate notebook [ferrite_navier_stokes_cavity.ipynb](ferrite_navier_stokes_cavity.ipynb) (in progress).

### Section 6.6: Post-Processor 

Resultsa are visualized using [paraview](https://www.paraview.org). 

### Section 7.5: Preliminary Results 

Here we show computational results for the magnitude of the velocity for Stokes flow. 

<img src="./stokes_3d_cylinder.png" width=400 />


## Section 6: References 
1. Mason-Sheffold-2014; 
2. Russel, Saville and Schowalter, “Colloidal Dispersions”; 
3. master thesis [Wei Fan](https://repository.tudelft.nl/record/uuid:acbfbd35-3ab1-4ec2-9a7b-b3098bdfcea9);
4. Tabilo-Munizags-2005: provides good overview.
5. Barnes-1994: Rheology of Emulsion - A review: provides good introduction
6. Barnes-2000: [Handbook_of_Rheology](https://ia601206.us.archive.org/4/items/HandbookOfRheology/Handbook_of_Rheology.pdf) - valuable resource; 
7. slides by D. Lahaye;

### Books 
1. Belitz, Grosch and Schieberle, <i>Food Chemistry</i>, Springer, 2009, [link](https://www.google.nl/books/edition/Food_Chemistry/xteiARU46SQC?hl=en&gbpv=0);
2. Hasenhuettl and Hartel, <i>Food Emulsifiers and Their Applications</i>, Springer, 2019, [link](https://www.google.nl/books/edition/Food_Emulsifiers_and_Their_Applications/5Ea9DwAAQBAJ?hl=en&gbpv=0); 


### Videos 
1. [What is an emulsion? by Silverston](https://www.youtube.com/watch?v=mBvKar6t1LY): mixing by high shear to reduce surface tension;  
2. [The emulsification process by Jacob Burton](https://www.youtube.com/watch?v=qnudmk_63r4): viscosity as a stabilizer; 
3. [What is an emulsion by Dow](https://www.youtube.com/watch?v=uWfdU92uPNY) phases separate to find state with lesser energy; 


```julia

```
