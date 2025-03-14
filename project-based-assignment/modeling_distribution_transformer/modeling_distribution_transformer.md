# Future Distribution Grids: A Case Study of Modeling the Magnetic Field, the Thermal Field and the Electrical Circuit Parameters in a Power Transformer

## Domenico Lahaye (TU Delft) and Jeroen Schuddebeurs (STEDIN Rotterdam)

## Section 1: Introduction 

The operation of electrical distrubution grids is changing fast. The infeed of renewables and the loading by battery charging is causing unforeseen challenges. It is expected that the introduction of new sources (e.g. solar panels) and loads (e.g. electrical vehicle charging) in the distribution grid which causes higher frequency harmonics. STEDIN is interested in monitoring the life-time expectancy of assests in the distribution grid. 

The live-time expectancy of distribution grids is closly linked to the degradation of the thermal insulating properties of the transformer and the reduction of the life-time of the transformer over time (change in heat absorption properties of the transformer oil and in the electricallly insulationg paper material) due to higher harmonics.

The finite element modeling of electrical power transformer is indispensable to chart the amplitude, frequency content and location of electromagnetic losses.   

We consider the same three-phase power distribution transformer considered earlier in the master thesis of Max van Dijk. (How to discribe the power rating of the transformer? Need to include more specs here). 

We wish to develop a model consisting of five sub-models. These components are described in the following. 

### Sub-model (1/5): A two-dimensional magnetic field - electrical circuit coupled model

The first submodel is a model for the quasi-stationary electromagnetic field in the transformer. The goal of this submodel is to compute the amplitude of the fundamental (base or ground) and higher harmonics of the B-field and the H-field. These <b>higher harmonics</b> are caused by higher harmonics in the current or voltage excitation, by the distortion of the harmonic fields due to magnetic saturation or by a combination of these effects.  

The quasi-stationary electromagnetic submodel consists of the following three components: 
1. a two-dimensional <b>quasi-stationary magnetic field model</b>. It is assumed that three-dimensional end effects can be neglected due to the lamination of the ferromagnetic core. A two-dimensional perpendicular current model for the magnetic field will therefore be solved. The magnetic field will be solved in terms of the perpendicular component of the vector potential $A_z$ as before in the EE4375 course. The ferromagnetic saturation and the electrical conductivity on the ferromagnetic core will be taken into account;
2. an <b>electrical circuit model</b> for the six voltage-driven (three high voltage and three low voltage) stranded (or winded) conductors;  
3. a model for the <b>field-circuit coupling</b> between the magnetic field and the electrical circuit in terms of the induced voltage. 

Due to the non-linearity of the ferromagnetic core (yoke), a linear superposition principle for the field will be likely to result in accurate results. Both a multi-harmonic (coupling fundamental and higher frequencies) and a transient formulation are expected to be carried out.

Challenges in developing this model include: 
1. the <b>uncertainty</b> in the in-plane  electrical conductivity $\sigma$ and the magnetic permeability $\mu$ of the laminated ferromagnetic core. The electrical conductivity is expected to varry with the temperature of the core. The magnetic permeability is expected to vary with the magnetic field strength to take magnetic saturation into account. [Magnetic hysteresis](https://en.wikipedia.org/wiki/Steinmetz%27s_equation) effects (multi-valued nature of the B-H characteristic) will be neglected in line with the two-dimensional simplification of the model. Out-of-plane plane variations of material parameters are neglected in the two-dimensional model; 
2. the <b>computational cost</b> required to solve the field-circuit coupled model (either in time or in frequency domain); 


### Sub-model (2/5): A core loss model
The second submodel is a model for the <b>ferromagnetic core</b> losses. The losses are caused by eddy current effect in the core. They can be computed from the B-field and the H-field computed in previous sub-model. The core is assumed to be laminated. 

In the laminated ferromagnetic core, losses have two contributions: a first [eddy-current loss](https://en.wikipedia.org/wiki/Magnetic_core#Core_loss) term and an second [iron-loss (Steinmetz)](https://en.wikipedia.org/wiki/Magnetic_core#Core_loss) term. These lossses are computed in two steps. First two-dimensional results (expressed in W/m) are computed in the modeling plane. These results are subsequently are subsequently multiplied with the length perpendicular denoted by $\ell_z$ (to obtain losses expressed in W).  
  
Challenges in developing this model include:
1. the uncertainty in setting the material characteristics for the application of the Steinmetz formula (unclear how sensivite results are w.r.t. these parameters);
2. the lack of experimental validation for the losses (as losses can only indirectly be deduced from measurements of temperature); 
3. challenges in modeling the lamination of the core; 

### Sub-model (3/5): A copper loss model  
The third submodel is model for the <b>copper losses</b> in the primary and secondary winding. The windings are assumed to be three-phase multi-turn windings (no eddy currents). Losses in these windings are caused by the Ohmic resistance of the copper material. They cause heating. These are modeled by terms proportional to $R \, I^2$.

In case that the winding is current driven, the current is known. In this case, Ohmic losses are straightforward to compute. In case that the winding is voltage driven, the current is computed by solving the field-circuit coupled model. Given the current, Ohmic losses can again be computed as before. 

### Sub-model (4/5): A three-dimensional thermal model
The fourth submodel solves for the <b>temperature</b> in the transformer and the vicinity of the transformer. First by solving a diffusion problem for the temperature, next by solving a convection-diffusion probem assuming a given divergence-free velocity again for the temperature, next by solving convervation of mass, momentum and energy assuming ideal gas for density, velocity and temperature of the air. 

Need to take into account the following factors: the heat absorption capabilities of the transformer oil, the cooling fins, the ventilation of the environment in which the transformer is placed, the outside temperature in which the substation is placed (day and night fluctuations) and the thermal insulation capabilities of the substation. 

Results of the temperature field can be compared with measured data. 

### Sub-model (5/5): An aging (life-time) model
A model for the degradation of the heat absorption capabilities of the transformer oil and the insulating capabilities of the transformer. See e.g. the [bathtub curve](https://en.wikipedia.org/wiki/Bathtub_curve) are related references to reliability engineering and deterioration modeling. 

## Work Performed in the Past 
1. Max van Dijk:  model in Comsol Multiphysics assuming linear superpostion to hold, see [master thesis text of Max van Dijk](https://repository.tudelft.nl/islandora/object/uuid%3A15b25b42-e04b-4ff2-a187-773bc170f061?collection=education); 
2. Gijs Lagerweij: start of transition from commercial software to public domain implementation in Julia, [previous assignments by Gijs Lageweij](https://github.com/gijswl/ee4375_fem_ta); 
3. Philip Soliman and Auke Schaap: fast finite element assembly and transient simulation in the presence of ferromagnetic saturation, [previous assignments by Philip Soliman and Auke Schaap](https://github.com/aukeschaap/am-transformers); 
4. Rahul Rane: extension of previous work by Philip Soliman and Auke Schaap with [an extensive set of results](https://github.com/rahulmrane/fem_future_distribution_grids); 

## Section 2: Literature Study as a Start on the Assignment 
Consult references listed below to become familiar with the project. Rephrase asssignment with own words. Find estimates for the copper losses, iron losses and temperature in various operating conditions of the transformer. Check understanding of the assignment with D. Lahaye. 

## Section 3: Two-Dimensional Coupled Magnetic Field Electrical Circuit Model    
Magnetic field computations will be performed to estimate the copper and iron losses in the transformer under in various operating conditions. 

Describe the mathematical model for the magnetic field computation in terms of the $z$-component of the magnetic vector potential. This component is oriented perpendicular to the 2D modeling plane. This formulation is referred to as the perpendicular current formulation. Formulate the PDE in frequency (at single and multiple frequencies) and time domain, the boundary conditions and the initial conditions. Rewrite the strong form of equations into weak form. 

Take electrical conductivity $\sigma$ and magnetic permeability $\mu = \mu_0 \, \mu_r$ of the laminated core into account. What is the effect of the lamination on the value of $\sigma$ and $\mu$. Lowering these values? We will consider two models to represent $\mu$. In the first model, $\mu$ is consider to be constant. No magnetic saturation is not taken into account. In the second model, $\mu$ changes with the strength of the magnetic field normB. The modeling of hysteresis is not taken into account. First use an analytical model for the non-linear material characteristic. Then use a BH-curved measured from data (with some form of spline interpolation between the data points). Compute derivative to arrive at a value for $\mu$. Results by Gijs Lagerweij show that the BH-curve (effective BH-curve, scaling with the RMS value) choosen by Max van Dijk in his thesis might not be the most appropriate ones. This requires further looking into.  

Derive an expression for the skin-depth and the value of the B-field and H-field in a thin layer close to the edge of the material. Compute the skin-depth for various relevant parameter values. $\delta = \sqrt{1/(\omega \, \sigma \, \mu)}$. High frequency will result in small skin-depth with large values of magnetic flux inside a small layer. Take this effect into account when generating the mesh (for both first and second order elements).

Describe primary and secondary windings in terms of number of turns, cross-section of the windings, insulation material. Describe current-driven windings and voltage-driven windings. 

Describe post-processing for B-field, H-field, skin depth in the core, core lossess and winding losses. Describe Steinmetz formula to compute core losses. See [previous work on losses by Gijs Lagerweij](https://github.com/gijswl/ee4375_fem_ta/blob/main/assignment2_fem2d/distribution_transformer/Distribution%20Transformer%20Stedin%20-%20Losses.ipynb).  

Describe operating conditions such as no load, short circuit (requires more text to describe these operating conditions). 

We will extend code developed in the EE4375 course and use Gridap, Ferrite or alternative as a reference for computations.

### 2D Magnetic Field Geometry Definition
Use GMSH (or alternative) to construct a parametric 2D CAD model magnetic field for the power transformer. Beware of introduction of parameter $\ell_z$ that describes length of the model in the $z$-direction. Extend with <b>fillets</b> to avoid corner singularities.  See [project by Gijs Lagerweij](https://github.com/gijswl/ee4375_fem_ta/tree/main/assignment2_fem2d).  

### 2D Mesh Generation 
Construct computational mesh for above 2D CAD model. See [project by Gijs Lagerweij](https://github.com/gijswl/ee4375_fem_ta/tree/main/assignment2_fem2d). Beware of: 
1. structured mesh or boundary refinement to capture eddy current effects in the core;
2. mesh the coil cross sections using rectangular elements and mesh remainder of geometry using triangular elements; 

### Material Parameter Setting 

#### Analytical BH-Curve Models 

#### Data Driven BH-Curve Models 
Given BH-curve, derive expression for $\mu(normB2)$ and $d\,\mu(normB2)/dnormB2$ using cubic splines.  

#### Parameter Curves for Iron Loss Computed 

### Coding Requirements 

#### FE Assemby Coding Requirements 
1. 2D first and second order triangular elements (see Assignments of Rahul and Savion). Ought to be rendered type-stable; 
2. mixed triangular and quadrilateral elements both for first and second order elements;
3. electrical circuit equation (see previous work by Gijs Lagerweij [notebook](https://github.com/gijswl/ee4375_fem_ta/blob/main/assignment2_fem2d/distribution_transformer/Distribution%20Transformer%20Stedin.ipynb)). How is time-step in transient simulations affected?   
4. harmonic balance method (amplitude of the vector potential component at various frequencies); 

#### FE Non-Linear Solve Coding Requirements
1. non-linear solve in time-harmonic computations: non-analytical residual (required), and thus fixed-point iterations; how to take circuit equation into account? See [previous work by Gijs Lagerweij](https://github.com/gijswl/ee4375_fem_ta/blob/main/assignment2_fem2d/distribution_transformer/Distribution%20Transformer%20Stedin.ipynb); 
2. non-linear solve in transient computations; 
Perform FEM assembly of the residual and the Jacobian (how to compute this Jacobian) at every time step.

#### Post-Processing Coding Requirements
1. write data to VTK-file for post-processing (done); 
2. post-processing using Paraview by integrating over surfaces and computing FFT of transient solution; 
3. post-processing in code by integrating computed vector potential component over surfaces of the coils to be able to set-up field-circuit coupled problem;

### Time-Harmonic Magnetic Field Electrical Circuit for the Magnetic Vector Potential Component A$_z$
   
#### Single frequency excitation, linear magnetic material, current driven coils
Assume that 
1. the electromagnetic field varries harmonically at a single frequency. This frequency is the base frequency (50 Hz) or an odd multiple of this frequency;   
2. the coils are current driven, i.e., the applied current density over the cross-section of the coils is given. Confirm that no electrical circuit equations are required. The case of voltage driven coils will be different. See previous work. 

Perform finite element simulations of the magnetic vector potential at various frequencies for various operating conditions (no-load and short-circuit) of the transformer. Tabulate copper and iron losses at various frequencies. 

#### Single frequency excitation, non-linear magnetic material, current driven coils
Next we assume the ferromagnetic material to be non-linear. Magnetic saturation is thus taken into account. Repeat above computation at various frequencies.

#### Single frequency excitation, linear magnetic material, voltage driven coils
Next assume the coils to be voltage driven. Derive expression for the electrical circuit equation including the coupling with the magnetic field. Repeat above computations. Verify how the computed current value differs from the case of current-driven 

#### Multiple frequency excitation
Extend from single frequency (fundamental frequency at 50 Hz) to take higher order harmonics (third, fifth and seventh harmonic, truncated Fourier series) into account. Even though excitation contains 50 Hz component only, higher order are generated by the material non-linearity. Obtain PDE for various amplitudes of the harmonics by testing with the various harmonics.   

Details on the harmonic balance method: 
1. see [wiki on harmonic balance method](https://en.wikipedia.org/wiki/Harmonic_balance). Illustrative example of how non-linearity induces higher harmonics is given; 
2. master thesis of [Alessandro Giust at Pavia University](https://thesis.unipd.it/retrieve/2b82a529-633a-4f25-8e4a-ea9fcd226a5a/Giust_Alessandro_tesi.pdf). Detailed explained. Relevant examples given; 
3. application to electromagnetic energy transformers [paper by De Gersem](http://bib.iem.rwth-aachen.de/IEMpublications/AltesBib/pub_690.pdf) (Krylov method extended to convolution products. Harmonics large due to corner singularities). [paper De Gersem, Van de Sande, Hameyer](https://www-emerald-com.tudelft.idm.oclc.org/insight/content/doi/10.1108/03321640110383852/full/pdf?title=strong-coupled-multiharmonic-finite-element-simulation-package). Application to higher order in transformer; 
4. Julia package [HarmonicBalance.jl](https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/background/harmonic_balance/); 

### Transient Magnetic Field Electrical Circuit for the Magnetic Vector Potential Component A$_z$
Extend above (formulation and finite element computations) from time-harmonic formulation to transient computations (by time-stepping on the vector potential $A_z$). 

#### Single frequency excitation, linear magnetic material, current driven coils
Assume current-driven windings and linear magnetic permeability. Perform time-stepping until effect of the initial conditions has disappeared and a time-oscillating solution is obtained. Use DifferentialEquations.jl. Repeat for various parameter settings, e.g. settings for the magnetic permeability, electrical conductivity and frequency of the applied current density; 

Post-process for the B-field, H-field, skin-depth (by plotting B-field along line perpendicular to the surface of the core), induced voltage (by surface integration over the volume), core loss, winding loss. Expect to see higher B-field peak values at higher frequencies due to the skin-effect; 

#### Single frequency excitation, non-linear magnetic material, current driven coils
Goal: extend from linear magnetic material to non-linear magnetic material (thus taking magnetic saturation into account). 

Coding requirements: place FE assembly inside of function that defines the ODE to be solved; 

Expected results: lower B-field peak values as higher frequencies. Verify whether single low (high) frequency of the current excitation induces high (low) frequencies in the potential.

#### Single frequency excitation, linear magnetic material, voltage driven coils
Goal: extend from current-driven coils (ODE-formulation) to voltage-driven coils (Differential-Algebraic Equation formulation);

Coding requirements: extend definition of the residual for the magnetic field to the residual for the electrical circuit and the field-circuit coupling term; 

Expected results: lowering the Ohmic resistance of the circuit should result in higher current values, higher field values and thus larger losses; 

#### Single frequency excitation,  non-linear magnetic material, voltage driven coils
Goal: extend to non-linear BH-curve. Perform FE assembly in space during time-stepping loop; 

### Compare time-harmonic and transient computations
1. Transform time-harmonic results to time-domain by multiplying with the time-domain part of the solution of the solution;
2. Transform time-domain to the frequency domain by Fast Fourier Transform? 

## Section 4:/ Two-Dimensional Copper and Iron Loss Model  
1. estimate losses from literature or approximate models (loading guide).
2. copper losses in the winding: [copper losses](https://en.wikipedia.org/wiki/Copper_loss); determined sperimentally using [short-circuit test](https://www.electricaleasy.com/2014/04/open-and-short-circuit-test-on-transformer.html); 
3. iron losses in the ferromagnetic core: [core loss](https://en.wikipedia.org/wiki/Magnetic_core#Core_loss): hysteresis losses and eddy-current losses; [open-circuit test](https://www.electricaleasy.com/2014/04/open-and-short-circuit-test-on-transformer.html); using [Steinmetz formula](https://en.wikipedia.org/wiki/Steinmetz%27s_equation); material dependent parameters might be difficult to recover; 
4. both copper losses and core losses expressed in W/m^3 such that value over entire volume is expressed in W;
5. Estimate reduction in losses due to heat absorption by the transformer oil; 
6. Estimate reduction of quality of the oil due to degradation of the pH value of the oil. Bathtub curve for the life-time of the transformer oil; 
7. run simulation of short-circuit and open-circuit test for a set of increasing frequencies, post-process for skin depth in the ferromagnetic core along a vertical line, post-process for relative importance of coil and core losses; 

## Section 5:/ Three-Dimensional Thermal Source Term 
Translate 2D copper and iron losses to equivalent 3D losses by (in first instances) multiplying by the length perpendicular to the model. More refined ways can be developed later.   

## Section 6:/ Three-Dimensional Thermal Field Computations

### 3D Thermal Field Geometry Definition
Construct a parametric 3D CAD model for the thermal field in the insulating material, cooling oil container, cooling fins and the environment (substation) of the transformer. Use GMSH, [Blender](www.blender.org) or any alternative to create STL-files (or alternative) files. 

### 3D Mesh Generation 
Use GMSH to construct computational mesh for above 3D CAD model.

### 3D Thermal Field Finite Element Modeling 
Compute temperature distribution with and without oil container. Look into infrastructure to interpolate from 2D to 3D mesh and to restrict from 3D mesh to 2D mesh. (In Gridap.jl, see the Lagrange multiplier example in the tests). Measured temperature on lateral surface of the oil container is max 200 degrees Celcius. Given no information on largest temperature inside of domain and duration of this temperature. Properties of [transformer oil](https://en.wikipedia.org/wiki/Transformer_oil).

1. document thermal conductivity, density and heat capacity of coil, core, insulating material and air domain;
2. compute thermal field; how to limit computational domain? What boundary conditions should be imposed at the boundaries;
3. how to take the presence of the insulating oil into account; break-down on oil due to overheating; 
4. how to take laminar/turbulent flow of oil due to buoyancy into account; using fluid dynamics?  
5. repeat for various values of the convective contribution; 

### Temperature Dependence of Magnetic Parameters and 2D Coupled Magnetic-Thermal Model
1. temperature dependence of [thermal conductivity](https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity#Temperature_dependence);  
2. temperature dependence of [magnetic permeability](https://www.sentesoftware.co.uk/mint-project/uploads/850554056.pdf); 
 
### Electrical Circuit Model 
1. more complex representation of the driving electrical circuit; allow for single phase and three phase electrical circuit;  

## Section 7:/ Future Work

### Joint Tasks for TA1 and TA2 to be hired 

- Consider three-phase distribution transformer operated by STEDIN
- Document construction of the ferromagnetic core, of primary and secondary winding and of tank for transformer oil: dimensions, materials used, number of turns and insulating material used for windings;
- Document thermal paramaters (heat capacity of the ferromagnetic core) and the thermal dependency of the magnetic permeability, electrical conductivity and properties of transformer oil;
- Document expected results of iron and copper losses and temperature from literature study;
- Study previous work by Gijs Lagerweij;

### Tasks for TA1
- 2D magnetic field (z-component of the magnetic vector potential)- electrical circuit (for voltage driven stranded conductor) simulations assuming load on the transformer; 
- validating previous work in various scenarios (no load, short circuit, other) using literature. Limited thermal field simulations; 
- No 3D simulations: use assumption that due to lamination, field in mainly a two-dimensional effect; the thermal field is, however, full dimnensional 3D field; 
- develop solver for coupled magnetic field, electrical circuit (and thermal field) simulations in both frequency and time domain by combining previously developed components; further develop either own code or make use of existing components in Julia (e.g. Ferrite or Gridap)
- translate 2D results from iron and copper losses into 3D sources for temperature taking into account the insulating capacity of transformer oil (part of heat losses generated will be absorbed by the oil) and cooling fins (improved heat transfer from transformer to surrounding);  
- Run solver in various scenarios;
- Document finding in a Jupyter notebook made available on github; 
 
### Tasks for TA2 
- 3D thermal field simulations assuming heat source modeled by TA1 or values from literature; 
- diffusion of heat, natural convection of heat, forced convection of heat by ventilation; does air humidy play a role here? 
- develop 3D parametric CAD models using e.g. Blender for transformer substation taking walls, floor, ceiling (and ventilation?) into account. Ventilation might require fluid flow simulations? 
- denerate 3D meshes (coarse, intermediate and fine); export STL generated by Blender, import oin GMSH and generate mesh; 
- Impose localized heat source or heat source found by TA1; 
- Run simulation for (e.g. Ferrite or Gridap)
- Run solver in various scenarios;
- Document finding in a Jupyter notebook made available on github; 

### Gridap specific issues 
1. how to accomodate electrical circuit equations such that coupling with other fields becomes possible? Extend function affineoperators() prior to linear solve? Need to be more specific here; 
2. how to interpolate copper and iron losses computed on 2D magnetic field mesh to source term how the 3D thermal field mesh?; 

### Recommendations by Gijs Lagerweij (in Dutch)
De volgende stap hangt heel erg af van het doel.

Als de wens is om met simpele modellen de verliezen (en later evt. temperatuurverdeling) in de trafo te modelleren, dan denk ik dat de volgende stap moet zijn om netjes alle materiaalparameters te bepalen.

De BH curve in mijn model is vrij eenvoudig (gebaseerd op een eenvoudige formule), terwijl deze in werkelijkheid misschien anders is.

Verder zijn de Steinmetz parameters uit Max’ thesis naar mijn idee helemaal niet juist. Als we deze goed kunnen bepalen moet het mogelijk zijn een stuk dichter bij de gemeten verliezen te komen, zonder dat er complexe 3D modellen gemaakt hoeven te worden.
 
Misschien het meest belangrijke, als we Gridap verder willen gebruiken, is de ontwikkeling van een betere interface voor de field-circuit coupling. Op dit moment is het door mijn implementatie namelijk niet mogelijk om gebruik te maken van sommige Gridap-onderdelen (zoals FEOperators), en daarmee wordt het modelleren van multi-field en transiente problemen een stuk lastiger. Dit zou bijvoorbeeld een extensie voor Gridap kunnen worden. De meeste code wordt op dit moment namelijk “verspild” aan de herhaaldelijke definitie van allerlei dingen in Gridap (integralen over de wikkelingen, voor de field-circuit coupling). Dit zou de code voor zowel TA1 als TA2 een stuk leesbaarder maken en fijner om mee te werken.

## References 
- Master thesis [Max van Dijk](https://repository.tudelft.nl/islandora/object/uuid%3A15b25b42-e04b-4ff2-a187-773bc170f061?collection=education)
- TA Ship of master student [Rahul Rane](https://github.com/rahulmrane/fem_future_distribution_grids);
- Comsol Multiphysics [Computing Losses in a Three Phase Transformer](https://www.comsol.com/blogs/computing-losses-in-a-three-phase-power-transformer).
- book [Transformers Practical Design Guide](https://books.google.nl/books?id=XQEWEAAAQBAJ&printsec=frontcover&dq=distribution+transformers&hl=en&sa=X&redir_esc=y#v=onepage&q=distribution%20transformers&f=false) Especially Table 2.1 on page 8 is of interest to us; 
- book [Transformer Engineering: Design, Technology and Diagnostics](https://www.google.nl/books/edition/Transformer_Engineering/aynOBQAAQBAJ?hl=en&gbpv=1&dq=kulkarni+transformer+design&printsec=frontcover) and similar books by same authors; 
- book [Magnetostatic Modelling of Thin Layers](https://link.springer.com/book/10.1007/978-3-319-77985-0) Figure 3.3 on page 27 provides a good illustration of how non-linear effects (due to magnetic saturation) cause higher harmonics in the magnetic field; 


```julia

```
