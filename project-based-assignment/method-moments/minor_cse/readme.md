# Magnetized Metalic Objects - Readme 

## Section 1: Introduction 

When a metalic object is placed inside a magnetic field, the object is magnetized. The objective of this project is to contribute to the development of a new simulator for the magnetization field inside the object. 

(Insert figures here). 

## Section 2: Project Description 

### Section 1.2: Problem Formulation 

<b>Computational Domain</b> Assume the computational domain $\Omega$ to denote a cube with lenght $L$, height $H$ and depth $D$ alligned with the coordinate axes. Then $0 \leq x \leq L$, $0 \leq y \leq H$ and $0 \leq z \leq D$ and $\Omega = [0,L] \times [0,H] \times [0,D]$. Typical values are $L = H = 1 \, \text{m}$ and $D = 0.1 \, \text{m}$. Assume $\mathbf{r} = (x,y,z)$ to denote the position vector inside $\Omega$. 

<b>Material Properties</b> Assume that $\Omega$ has a [magnetic susceptibility](https://en.wikipedia.org/wiki/Magnetic_susceptibility) denoted by $\chi_{mag}$. Then $\Omega$ will magnetize in an external field and thus mimmic a metalic plate. Here it will be suffucient to assume that the plate is homogeneous and that therefore $\chi_{mag}$ is constant. Assume that $200 \leq \chi_{mag} \leq 200,000$ (dimensionless). Non-homogeneous plate can be modeled assuming that $\chi_{mag}$ is piecewise constant. Non-linear magnetization effects are excluded in this project.   

<b>Problem Description</b> Assume $\Omega$ to be placed in external (assumed known or given) [magnetic field](https://en.wikipedia.org/wiki/Magnetic_field). Let the external magnetic field be denoted by ${\mathbf H}_{ext}(\mathbf{r})$ in units Ampere per meter. Then ${\mathbf H}_{ext}(\mathbf{r})$ is a vector field with three components. This can be expressed as ${\mathbf H}_{ext}(\mathbf{r}) = \left( H_{ext,x}(\mathbf{r}), H_{ext,y}(\mathbf{r}), H_{ext,z}(\mathbf{r})\right)$ or as $ \mathbf{H}_{ext}(\mathbf{r}) = H_{ext,x}(\mathbf{r}) \mathbf{i} + H_{ext,y}(\mathbf{r}) \mathbf{j} + H_{ext,z}(\mathbf{r}) \mathbf{k}$. In case that the external magnetic field is constant and aligned with the $x$-direction, we can write that $\mathbf{H}_{ext} = \left( H_{0,ext}, 0,0\right)$ where $H_{0,ext}$ is a constant.

The goal of the project is to compute the [magnetization vector](https://en.wikipedia.org/wiki/Magnetization) in the plate. This vector field is denoted by $\mathbf{M}(\mathbf{r})$ in units Ampere per meter. This is a vector field with three components. The vector field to be computed can thus be written as $\mathbf{M}(\mathbf{r}) = \left( M_x(\mathbf{r}), M_y(\mathbf{r}), M_z(\mathbf{r})\right)$. 

The physical problem can thus be stated as: given the scalar $\chi_{mag}$ in $\Omega$ and given the vector field ${\mathbf H}_{ext}(\mathbf{r})$, compute the vector field ${\mathbf M}(\mathbf{r})$.

<b>Mathematical Model</b> (Provide more explanation on the eqiuation for $\mathbf{M}(\mathbf{r})$. State that total magnetic field ${\mathbf H}_{tot}$ in $\Omega$ is sum of two components. The external field ${\mathbf H}_{ext}$ and the induced field ${\mathbf H}_{M}$. We thus have that ${\mathbf H}_{ext} + {\mathbf H}_{M} = {\mathbf H}_{tot}$ where ${\mathbf H}_{M} =  \text{grad}_{\mathbf{r}} \phi_M = \nabla_{\mathbf{r}} \phi_M$ and where $\phi_M$ can be related to $\mathbf{M}$) We wish to compute $\mathbf{M}(\mathbf{r})$ by solving a partial differential equation in which the external field ${\mathbf H}_{ext}(\mathbf{r})$ acts as a source term (i.e., ${\mathbf H}_{ext}(\mathbf{r})$ appears in the right-hand side of the equation). This differential equation can be written as  

$$
{\cal D} \left[ \mathbf{M}(\mathbf{r}) \right] = {\mathbf H}_{ext}(\mathbf{r}) \text{ on } \Omega 
$$

where ${\cal D}$ is the [differential operator](https://en.wikipedia.org/wiki/Differential_operator) 

$$
{\cal D}\left[ \mathbf{M}(\mathbf{r}) \right] = \frac{1}{\chi_{mag}} \mathbf{M}(\mathbf{r}) - 
\text{grad}_{\mathbf{r}} \, \text{div}_{\mathbf{r}} 
\frac{1}{4 \pi}  \int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' = 
\frac{1}{\chi_{mag}} \mathbf{M}(\mathbf{r}) - 
\nabla_{\mathbf{r}} \, \nabla_{\mathbf{r}} \cdot \frac{1}{4 \pi}
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, . 
$$

This differential equation can be derived from [Biot-Savart Law](https://en.wikipedia.org/wiki/Biot–Savart_law). The differential operator ${\cal D}$ contains an integral term. The differential equation that determines $\mathbf{M}(\mathbf{r})$ is therefore an [integro-differential equation](https://en.wikipedia.org/wiki/Integro-differential_equation). The above integro-differential equation is a system of three coupled equations for the three components of the magnetization $\mathbf{M}(\mathbf{r})$. This equation is solved on $\Omega$ only. The magnetization $\mathbf{M}(\mathbf{r})$ remains zero outside $\Omega$. The kernel in the second term of ${\cal D}$ can be written as 

$$
{\cal K}(\mathbf{r},\mathbf{r}') = \frac{1}{4 \pi}  \frac{1}{\|\mathbf{r}' - \mathbf{r} \|} \, . 
$$

This kernel is singular for $\mathbf{r}' = \mathbf{r}$. The convolution of $\mathbf{M}(\mathbf{r})$ and ${\cal K}(\mathbf{r},\mathbf{r}')$ can be written as 

$$
{\cal K}(\mathbf{r},\mathbf{r}') * \mathbf{M}(\mathbf{r}) = \frac{1}{4 \pi} \int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, . 
$$

For future reference, we write the integro-differential equation for $\mathbf{M}(\mathbf{r})$ in explicit form. The first of the three scalae equations is 

$$
\frac{1}{\chi_{mag}} M_x(\mathbf{r}) - 
\frac{\partial}{\partial x} \, \nabla_{\mathbf{r}} \cdot \frac{1}{4 \pi}
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' = H_{ext,x}(\mathbf{r}) \, , 
$$

where 

$$
\nabla_{\mathbf{r}} \cdot 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' = 
\frac{\partial}{\partial x} \frac{1}{4 \pi} \int_{\Omega} \frac{M_x(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' + 
\frac{\partial}{\partial y} \frac{1}{4 \pi} \int_{\Omega} \frac{M_y(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' 
+
\frac{\partial}{\partial z} \frac{1}{4 \pi} \int_{\Omega} \frac{M_z(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, ,  
$$ 

and simarly for $M_y(\mathbf{r})$ and $M_z(\mathbf{r})$.   

### Section 2.2: Mesh Generation and Shape Functions 

We wish to compute the magnetization vector $\mathbf{M}(\mathbf{r})$ inside the metallic plate $\Omega$. We therefore generate a mesh on $\Omega$. 

<b>Mesh Generation</b> We denote the mesh by $\Omega^h$. We assume that the mesh consists of tetrahedral elements only. We denote the number of nodes, edges, facets and elements of $\Omega^h$ by $N_n$, $N_{ed}$, $N_f$ and $N_e$, respectively. Assume the elements to be denoted by $\left\{ P_{\alpha} | 1 \leq \alpha \leq N_e \right\}$. The union of all elements $P_{\alpha}$ forms the entire domain of computation. That is, we have that $\cup P_{\alpha} | 1 \leq \alpha \leq N_e = \Omega$. We will use this elementary fact in the computation of the matrix $\underline{\underline{A}}$ and the right-hand side vector $\mathbf{b}$.    

We assume that the information to decompose an element $P_{\alpha}$ into a set of facets, a facet into a set of edges and an edge into a set of points to be available (representation of the mesh $\Omega^h$ as a directed a-cyclic graph (DAG) with corresponding operations to find parent nodes and child nodes. See e.g. [wikipedia entry on polygon mesh](https://en.wikipedia.org/wiki/Polygon_mesh).

<b>Shape Functions</b> The mesh $\Omega^h$ allows to define a set of nodal P1 linear Lagrange basis functions $\{ \phi_i(\mathbf{r}) | 1 \leq i \leq N _n \}$ (see classical finite elements method). The functions $\phi_i(\mathbf{r})$ are scalar functions. To be able to cast a vector-valued partial differential equation in weak or variational form, we define $3 \, N_n$ vector-valued shape functions $\boldsymbol{\phi}_{3\,i-2} = \left( \phi_i, 0,0\right)$, 
$\boldsymbol{\phi}_{3i-1} = \left( 0, \phi_i, 0\right)$ and 
$\boldsymbol{\phi}_{3i} = \left( 0, 0, \phi_i \right)$ for $1 \leq i \leq N_n$.  

More information on the mesh generation is provided in the [notebook](./mom_3d_plate_gmsh.ipynb). 

### Section 3.2: Galerkin Method    

Here we describe the [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method) that allows to convert the system of integral partial differential-equations for $\mathbf{M}(\mathbf{r})$ into a weak or variational formulation. This variational formulation will allow a spatial discretization. 

(see also example of the weak formulation of the Poisson equation as [wiki on weak formulation](https://en.wikipedia.org/wiki/Weak_formulation).) 

Assume $g(\mathbf{r})$ and $h(\mathbf{r})$ to be scalar functions on $\Omega$ such that 

$$
g(\mathbf{r}) = \int_{\Omega} \frac{h(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, . 
$$

Then by decomposing $\Omega$ into a set of volumes $P_{\alpha}$, we find that 

$$
\begin{eqnarray}
\int_{\Omega} g(\mathbf{r}) \, d\Omega = 
\sum_{\alpha=1}^{N_e}  \int_{P_{\alpha}} g(\mathbf{r}) \, d\Omega 
& = & \sum_{\alpha,\beta=1}^{N_e}  \int_{P_{\alpha}} \int_{P_{\beta}} \frac{h(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, d\Omega \\ 
& = & \sum_{\alpha,\beta=1}^{N_e}  \int_{\mathbb{R}^3} \int_{\mathbb{R}^3} \frac{h(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, \mathbb{1}_{P_{\beta}} \, d\Omega' \, \mathbb{1}_{P_{\alpha}} \, d\Omega
\end{eqnarray}
$$   

**Expansion of Magnetization Components** The numerical approximation of the magnetization components can be expanded in as 

\begin{eqnarray}
M_x(\mathbf{r}) = \sum_{i=1}^{N_n} c_i \, \phi_i(\mathbf{r}) \\ 
M_y(\mathbf{r}) = \sum_{i=1}^{N_n} d_i \, \phi_i(\mathbf{r}) \\
M_z(\mathbf{r}) = \sum_{i=1}^{N_n} e_i \, \phi_i(\mathbf{r}) 
\end{eqnarray}

Note that components $M_x(\mathbf{r})$, $M_y(\mathbf{r})$ and $M_z(\mathbf{r})$ are each expanded in the same basis with proper expansion coefficients. In case of linear Lagrange shape functions on tetrahedra, the element $P_{\alpha}$ has 4 local degrees of freediom per components, thus in total $3*4=12$ local degrees of freedom. Mapping from local to global degrees of freedom. The expansion coefficients $c_i$, $d_i$ and $e_i$ can stacked into a global vector $\mathbf{u}$ of dimension $3 \, N_n$ where $\mathbf{u}^T = \begin{pmatrix} \mathbf{c}^T \mathbf{d}^T \mathbf{e}^T \end{pmatrix}$. This global vector is the solution of the Galerkin-MoM linear system 

$$
\underline{\underline{A}} \, \mathbf{u} = \mathbf{b}
$$

where the matrix $A$ is a dense $3 \, N_n$-by-$3 \, N_n$ matrix consisting of a stiffness part $A^{(2)}$ and a mass matrix part $A^{(1)}$.  

**Discrete Weak Form**

$$
\int_{\Omega} {\cal D} \left[ \mathbf{M}(\mathbf{r}) \right] \cdot \boldsymbol{\phi}_k \, d\Omega = 
\int_{\Omega} \mathbf{H}_{ext} \cdot \boldsymbol{\phi}_k \, d\Omega 
\hspace{1cm} \forall 1 \leq k \leq N_n 
$$ 

or more explicitly 

$$
\frac{1}{\chi_{mag}} \int_{\Omega} \mathbf{M}(\mathbf{r}) \cdot \boldsymbol{\phi}_k \, d\Omega + 
\int_{\Omega} \left[ \nabla_{\mathbf{r}} \, \nabla_{\mathbf{r}} \cdot 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \right] \cdot \boldsymbol{\phi}_k \, d\Omega = \int_{\Omega} \mathbf{H}_{ext} \cdot \boldsymbol{\phi}_k \, d\Omega
$$ 

and thus after integration by parts on the second term 

$$
\frac{1}{\chi_{mag}} \int_{\Omega} \mathbf{M}(\mathbf{r}) \cdot \boldsymbol{\phi}_k \, d\Omega - 
\int_{\Omega} \left[ \nabla_{\mathbf{r}} \cdot \boldsymbol{\phi}_k \right]
\left[\nabla_{\mathbf{r}} \cdot 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \right] d\Omega = \int_{\Omega} \mathbf{H}_{ext} \cdot \boldsymbol{\phi}_k \, d\Omega
$$

resulting in $3 \, N_e$ equations for the $3 \, N_e$ components of the vector of expansion coefficients $\mathbf{u}$. The mass matrix and load vector can be assembled by a loop over elements as in classical Galerlin FEM. The stiffness matrix requires dedicated attention due to the integral term.

### Section 4.2: Element-by-element Assembly of Stiffness Matrix

**Local Representation**

Assume that $P_{\alpha}$ has $4$ nodes. Assume that $\phi_k({\mathbf r}) = a_k x + b_k y + c_k y + d_k$ is a (scalar) basis function for $1 \leq k \leq N_e^{\alpha}$.   

**Divergence of Vector Shape Function**

We have that 

$$
\nabla_{\mathbf{r}} \cdot \boldsymbol{\phi}_k = a_k \text{ if } \mod{k,3} = 0 \\
\nabla_{\mathbf{r}} \cdot \boldsymbol{\phi}_k = b_k \text{ if } \mod{k,3} = 1 \\
\nabla_{\mathbf{r}} \cdot \boldsymbol{\phi}_k = c_k \text{ if } \mod{k,3} = 2
$$

**Bilinear Form**

$$
\underline{\underline{A}}^{(2)} = \sum_{\alpha} \underline{\underline{A}}_{\alpha}^{(2)} = \sum_{\alpha\beta} \underline{\underline{A}}_{\alpha\beta}^{(2)} \text{ of size } 3 \, N_n \text{ by } 3 \, N_n 
$$

where the sum over $\alpha$ (destination) and $\beta$ (source) discretize the integral over $\mathbf{r}$ and $\mathbf{r}'$, respectively. In case that $\alpha = \beta$, self-interaction terms are considered. 

$$
\underline{\underline{A}}_{loc,\alpha\beta}^{(2)} = 
\begin{pmatrix} 
A_{\alpha\beta}^{xx} & A_{\alpha\beta}^{xy} &  A_{\alpha\beta}^{xz} \\
A_{\alpha\beta}^{yx} & A_{\alpha\beta}^{yy} &  A_{\alpha\beta}^{yz} \\ 
A_{\alpha\beta}^{zx} & A_{\alpha\beta}^{zy} &  A_{\alpha\beta}^{zz}
\end{pmatrix}
\text{ of size } 12 \text{ by } 12
$$

where (row determined by the expansion of the magnetization, columns determined by the test function) 
(assuming tetrahedra and $N_e^{\alpha}=4$). 

$$
A_{\alpha\beta}^{xx} = 
\begin{pmatrix} 
A_{\alpha\beta}^{xx,11} & A_{\alpha\beta}^{xx,12} & A_{\alpha\beta}^{xx,13} & A_{\alpha\beta}^{xx,14} \\
A_{\alpha\beta}^{xx,21} & A_{\alpha\beta}^{xx,22} & A_{\alpha\beta}^{xx,23} & A_{\alpha\beta}^{xx,24} \\
A_{\alpha\beta}^{xx,31} & A_{\alpha\beta}^{xx,32} & A_{\alpha\beta}^{xx,33} & A_{\alpha\beta}^{xx,34} \\
A_{\alpha\beta}^{xx,41} & A_{\alpha\beta}^{xx,42} & A_{\alpha\beta}^{xx,43} & A_{\alpha\beta}^{xx,44} 
\end{pmatrix}
\text{ of size } 4\text{ by } 4 
$$

where 

$$
\begin{eqnarray}
A_{\alpha\beta}^{xx,11} & =  & - \int_{P_{\alpha}} \frac{\partial \phi_1({\mathbf r})}{\partial x} \, \frac{\partial}{\partial x} \int_{P_{\beta}} \frac{\phi_1({\mathbf r})}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, d\Omega \\
& =  & - \int_{\mathbb{R}^3} \frac{\partial}{\partial x} \left[ \phi_1({\mathbf r}) \mathbb{1}_{P_{\alpha}} \right] \frac{\partial}{\partial x} \left[ \int_{\mathbb{R}^3} \frac{\phi_1({\mathbf r})}{\|\mathbf{r}' - \mathbf{r} \|} \, \mathbb{1}_{P_{\beta}} \, d\Omega' \right] \,d\Omega
\end{eqnarray}
$$ 

and similar for other components. Introduce indicator functions on $P_{\alpha}$ and $P_{\beta}$ and integrate over entire $\mathbb{R}^3$. This will allow to use the divergence theorem and the fact that the basis functions centered on $P_{\alpha}$ and $P_{\beta}$ are zero sufficiently far away from $P_{\alpha}$ and $P_{\beta}$.     

**Numerical Computation of Matrix Elements**

Use of hcubature (or alternative) to compute the matrix elements. Possibly regularization in the kernel;  

**Analytical Method of Matrix Elements**

1. [homogeneous functions](https://en.wikipedia.org/wiki/Homogeneous_function#) 
2. Theorem of [Euler](https://en.wikipedia.org/wiki/Homogeneous_function#Euler's_theorem) for homogeneous functions, reduce dimensionality of integration by one and examples;  
3. approach by Rosen/Cormack  
4. Reduce 6D integrals to 4D integrals
5. reduction of inner integral from 3D to 2D: interchange order of differentiation and integration, obtain two terms, apply Euler on one term;
6. reduction of outer integral from 3D to 2D: derivate of linear shape function is a constant. Apply Green-Gauss divergence theorem; 

\begin{equation}
\int_a^b f(x) \, dx = F(x) |_{x=a}^{x=b} = F(b) - F(a)
\end{equation}

## Section 3: Possible Project Roadmaps  

### Section 1.3: Problem Formulation 

<b>Physical Problem</b> Study problem formulation. Study how the magnetization vector $\mathbf{M}(\mathbf{r})$ change with  of the plate, dimensions of the plate, and the direction of the external magnetic field. Use to this end literature or a reference code.

### Section 2.3: Solution Method  

Explore weak formulation, weighted residual method, Galerkin approximation and linear system formulation for a sequence of problems of increasing complexity.  

1. Poisson equation;
2. one-dimensional Fredholm integro-differential equation;
3. see grad-div equation without kernel;

### Section 3.3: Assembly of the Mass Matrix and Load Vector 

1. loop over elements in the mesh to assemble the mass matrix and the lopad vector. See notebook [mom_hcubature](./mom_hcubature.ipynb); 

### Section 4.3: Numerical Assembly of the Stiffness Matrix 

1. use of quadrature implemented in [hcubature.jl](https://github.com/JuliaMath/HCubature.jl) for integration in 1D (possibly singular, look into number of function evaluatiohs), 2D (reference triangle, coordinate transformation, general triangle), 3D (reference tetrahedra, coordinate transformation, general tetrahedra), 4D (by calling hcubature for 2D twice) and 6D (by calling hcubature for 2D twice). See notebook [mom_hcubature](./mom_hcubature.ipynb);
2. formulate and compute 6D integrals for tetra/tetra interaction assuming no parallel facets. $P_{\alpha}$ with nodes ${\mathbf r}_1 = (0,0,0)$, ${\mathbf r}_2 = (1,0,0)$, ${\mathbf r}_3 = (0,1,0)$ and ${\mathbf r}_4 = (0,0,1)$. 
$P_{\beta}$ with nodes ${\mathbf r}_1 = (0,0,0)$, ${\mathbf r}_2 = (1,0,0)$, ${\mathbf r}_3 = (0.5,1,0)$ and ${\mathbf r}_4 = (0.5,0.5,1)$. Verify that no two faces are parallel to each other;
3. extend previous case to parallel faces. $P_{\alpha}$ reference triangle. $P_{\beta}$ shift of $P_{\alpha}$;
4. apply Euler theorem for homogeneous function and compute resulting 4D integrals numerically; 
5. loop over elements over the mesh to compute the load vector. Assume constant external field. Provide more details (expression of the local vector per element) here. See notebook [mom_3d_plate_gmsh](./mom_3d_plate_gmsh.ipynb); 
6. loop over elements over the mesh to compute the mass matrix. Provide more details (expression of the local matrix per element) here;
7. use of hcubature to compute the stiffness matrix elements. Possibly regularization in the kernel;

### Section 5.3: Analytical Assembly of the Stiffness Matrix 

1. implement point-edge interactions (cfr. seperate notebook); 
2. implement edge-edge interactions  (cfr. seperate notebook); 
2. implement point-facet interactions  (cfr. seperate notebook);

## Section 4: The Julia Programming Language

### Introductory Material
- Elementary introduction: [Thinking Julia](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html);
- Course: [Scientific-Programming-in-Julia](https://juliateachingctu.github.io/Scientific-Programming-in-Julia/stable/);  
- Aalto Short Course: [julia-introduction](https://github.com/AaltoRSE/julia-introduction); 
- Video Collection by Chris Rackauckas: [link](https://www.youtube.com/playlist?list=PLCAl7tjCwWyGjdzOOnlbGnVNZk0kB8VSa) 

### Project Specific Packages 

- Numerical Integral using hcubature 
- Mesh Generation using GMSH 
- Variational Formulation and Finite Element Methods: [Minfem.jl](https://minfem.github.io/MinFEM.jl/stable/) 

## References 

1. Morandi (method of moments);
2. [Introduction to Numerical Methods for Variational Problems](https://link.springer.com/book/10.1007/978-3-030-23788-2) by Hans Petter Langtangen and Kent-Andre Mardal. The [book](https://link.springer.com/book/10.1007/978-3-030-23788-2) is freely available; 
3. [Wolfgang Bangerth's video lectures](https://www.math.colostate.edu/~bangerth/videos.html); 
4. [wiki Finite Element Method](https://en.wikipedia.org/wiki/Finite_element_method): Section 3 for the weak form and Section 4 for the finite element discretization;  
5. [Comsol Multiphysics Finite Element Method](https://www.comsol.com/multiphysics/finite-element-method): more information and illustrations; 
6. [Comsol Multiphysics Brief Introduction to the Weak Form](https://www.comsol.com/blogs/brief-introduction-weak-form): good introduction to a theoretical concept that provides a basis for the finite element method; 
7. [Ferrite Introduction to FEM](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/fe_intro/)




```julia

```
