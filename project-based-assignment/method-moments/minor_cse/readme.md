# Magnetized Metalic Objects - Readme 

## Section 1: Introduction 

When a metalic object is placed inside a magnetic field, the object is magnetized. The objective of this project is to contribute to the development of a new simulator for the magnetization field inside the object. 

(Insert figures here). 

## Section 2: Project Description 

### Section 1.2: Problem Formulation 

<b>Physical Problem</b> Assume metalic plate to be placed in external (assumed known or given) [magnetic field](https://en.wikipedia.org/wiki/Magnetic_susceptibility). Let the external magnetic field be denoted by ${\mathbf H}_{ext}(\mathbf{r})$. Then ${\mathbf H}_{ext}(\mathbf{r})$ is a vector field with three components. This can be expressed as ${\mathbf H}_{ext}(\mathbf{r}) = \left( H_{ext,x}(\mathbf{r}), H_{ext,y}(\mathbf{r}), H_{ext,z}(\mathbf{r})\right)$ or as $ \mathbf{H}_{ext}(\mathbf{r}) = H_{ext,x}(\mathbf{r}) \mathbf{i} + H_{ext,y}(\mathbf{r}) \mathbf{j} + H_{ext,z}(\mathbf{r}) \mathbf{k}$. Assume that the plate has a [magnetic susceptibility](https://en.wikipedia.org/wiki/Magnetic_field) denoted by $\chi_{mag}$. The goal is to compute the [magnetization vector](https://en.wikipedia.org/wiki/Magnetization) in the plate denoted by $\mathbf{M}(\mathbf{r})$. This is a vector field with three components. The vector field to be computed can thus be written as $\mathbf{M}(\mathbf{r}) = \left( M_x(\mathbf{r}), M_y(\mathbf{r}), M_z(\mathbf{r})\right)$. 

<b>Mathematical Model</b> (requires reformulation) (see Poisson equation, see one-dimensional Fredholm integro-differential equation, see e.g. paper Morandi) Here we describe how the above problem can be translated (captured) into a mathematical model. 

Assume given a metalic object with magnetic susceptibility $\chi_{mag}$ (dimensionless, value between $20$ and $200.000$) and with volume $\Omega \subset \mathbb{R}^3$ (typical dimensions?) placed in a given external magnetic field $\mathbf{H}_{ext}$. In case that the external magnetic field is in the $x$-direction $\mathbf{H}_{ext} = \left( H_{ext}, 0,0\right)$ (value?, units?). Then the object is magnitized, meaning that the magnetization vector $\mathbf{M}(\mathbf{r})$ is non-zero inside $\Omega$. We wish to compute $\mathbf{M}(\mathbf{r})$ by solving the linear vector-valued grad-div volume integro-differential equation with singular kernel given by 

$$
\mathbf{A} \left[ \mathbf{M}(\mathbf{r}) \right] = {\mathbf H}_{ext}
$$

where 

$$
\mathbf{A}\left[ \mathbf{M}(\mathbf{r}) \right] = \frac{1}{\chi_{mag}} \mathbf{M}(\mathbf{r}) - 
\text{grad}_{\mathbf{r}} \, \text{div}_{\mathbf{r}} 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' = 
\frac{1}{\chi_{mag}} \mathbf{M}(\mathbf{r}) - 
\nabla_{\mathbf{r}} \, \nabla_{\mathbf{r}} \cdot 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega'
$$

The magnetization $\mathbf{M}$ remains zero outside $\Omega$, i.e., on $\mathbb{R}^3 \setminus \Omega$ (the equivalent of homogeneous Dirichlet boundary conditions). The integro-differential equation above is a system of three coupled equations for the three components of the magnetization $\mathbf{M}(\mathbf{r})$. 

For future reference, we write the equation in more explicit form. The first of the three scalae equations is 

$$
\frac{1}{\chi_{mag}} M_x(\mathbf{r}) - 
\frac{\partial}{\partial x} \, \nabla_{\mathbf{r}} \cdot 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' = H_{ext,x} \, , 
$$

where 

$$
\nabla_{\mathbf{r}} \cdot 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' = 
\frac{\partial}{\partial x} \int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' + 
\frac{\partial}{\partial y} \int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' 
+
\frac{\partial}{\partial z} \int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, . 
$$ 

### Section 2.2: Mesh Generation and Shape Functions 

We wish to compute the magnetization vector $\mathbf{M}(\mathbf{r})$ inside the object $\Omega$. We therefore generate a mesh on $\Omega$. 

<b>Mesh Generation</b> We denote the mesh by $\Omega^h$. We assume that the mesh consists of tetrahedral elements only. We denote the number of nodes, edges, facets and elements of $\Omega^h$ by $N_n$, $N_{ed}$, $N_f$ and $N_e$, respectively. Assume the elements to be denoted by $\left\{ P_{\alpha} | 1 \leq \alpha \leq N_e \right\}$. The union of all elements $P_{\alpha}$ forms the entire domain of computation. That is, we have that $\cup P_{\alpha} | 1 \leq \alpha \leq N_e = \Omega$. We will use this elementary fact in the computation of the matrix $\underline{\underline{A}}$ and the right-hand side vector $\mathbf{b}$.    

We assume that the information to decompose an element $P_{\alpha}$ into a set of facets, a facet into a set of edges and an edge into a set of points to be available (representation of the mesh $\Omega^h$ as a directed a-cyclic graph (DAG) with corresponding operations to find parent nodes and child nodes. See e.g. [wikipedia entry on polygon mesh](https://en.wikipedia.org/wiki/Polygon_mesh).

<b>Shape Functions</b> The mesh $\Omega^h$ allows to define a set of nodal P1 linear Lagrange basis functions $\{ \phi_i(\mathbf{r}) | 1 \leq i \leq N _n \}$ (see classical finite elements method). The functions $\phi_i(\mathbf{r})$ are scalar functions. To be able to cast a vector-valued partial differential equation in weak or variational form, we define $3 \, N_e$ vector-valued shape functions $\boldsymbol{\phi}_{3\,i-2} = \left( \phi_i, 0,0\right)$, 
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
\int_{\Omega} \mathbf{A} \left[ \mathbf{M}(\mathbf{r}) \right] \cdot \boldsymbol{\phi}_k \, d\Omega = 
\int_{\Omega} \mathbf{H}_{ext} \cdot \boldsymbol{\phi}_k \, d\Omega 
\hspace{1cm} \forall 1 \leq k \leq N_e 
$$ 

or more explicitly 

$$
\frac{1}{\chi_{mag}} \int_{\Omega} \mathbf{M}(\mathbf{r}) \cdot \boldsymbol{\phi}_k \, d\Omega + 
\int_{\Omega} \left[ \nabla_{\mathbf{r}} \, \nabla_{\mathbf{r}} \cdot 
\int_{\Omega} \frac{\mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \right] \cdot \boldsymbol{\phi}_k \, d\Omega = \int_{\Omega} \mathbf{H}_{ext} \cdot \boldsymbol{\phi}_k \, d\Omega
$$ 

and thus after integration by parts 

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

### Section 2.3: Focus on Numerical Computations

1. use of quadrature implemented in [hcubature.jl](https://github.com/JuliaMath/HCubature.jl) for integration in 1D (possibly singular, look into number of function evaluatiohs), 2D (reference triangle, coordinate transformation, general triangle), 3D (reference tetrahedra, coordinate transformation, general tetrahedra), 4D (by calling hcubature for 2D twice) and 6D (by calling hcubature for 2D twice). See notebook [mom_hcubature](./mom_hcubature.ipynb);
2. formulate and compute 6D integrals for tetra/tetra interaction assuming no parallel facets. $P_{\alpha}$ with nodes ${\mathbf r}_1 = (0,0,0)$, ${\mathbf r}_2 = (1,0,0)$, ${\mathbf r}_3 = (0,1,0)$ and ${\mathbf r}_4 = (0,0,1)$. 
$P_{\beta}$ with nodes ${\mathbf r}_1 = (0,0,0)$, ${\mathbf r}_2 = (1,0,0)$, ${\mathbf r}_3 = (0.5,1,0)$ and ${\mathbf r}_4 = (0.5,0.5,1)$. Verify that no two faces are parallel to each other;
3. extend previous case to parallel faces. $P_{\alpha}$ reference triangle. $P_{\beta}$ shift of $P_{\alpha}$;
4. apply Euler theorem for homogeneous function and compute resulting 4D integrals numerically; 
5. loop over elements over the mesh to compute the load vector. Assume constant external field. Provide more details (expression of the local vector per element) here. See notebook [mom_3d_plate_gmsh](./mom_3d_plate_gmsh.ipynb); 
6. loop over elements over the mesh to compute the mass matrix. Provide more details (expression of the local matrix per element) here;
7. use of hcubature to compute the stiffness matrix elements. Possibly regularization in the kernel;

### Section 3.3: Focus on Analytical Computations  

1. implement point-edge interactions (cfr. seperate notebook); 
2. implement edge-edge interactions  (cfr. seperate notebook); 
2. implement point-facet interactions  (cfr. seperate notebook);

## References 

1. Morandi




```julia

```
