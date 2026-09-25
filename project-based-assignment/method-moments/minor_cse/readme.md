# Magnetized Metalic Objects - Readme 

### How to further extend these notes 

**Grad-Div Equation**

1. See [MFEM weak gradient example](https://mfem.org/fem_weak_form/). 
2. See [Moose grad-div example](https://mooseframework.inl.gov/syntax/MFEM/Grad-Div.html). Not sure whether FENICS or dealII provides similar examples.  
3. Using search term <i>weak gradient form equation</i> and grad-div problem. 

**Singular Reisz Kernel**
1. extend notes on the singularity of the kernel; 

**Mesh Data Structure**
1. struct edge with field node holding global node index (extend what currently exists);
2. struct face with first field holding global node index and second field holding global edge index (extend what currently exists);
3. struct element with first field holding global node index, second field holding global edge index and third field holding global face index; 
4. the mesh has a field Points that holds an array of length npoints of type Point3D, similar for edges, faces and elements; 

**Two Triangular Facets Test Case** 

Adapt mesh structure (assignment for students) to 
1. remove node duplication;  
2. add element - face connectivity, add face - node connectivity, add face normal extracted from GMSH;  


**Extend double loop over elements to 4-fold loop over faces belonging to element or elementp (avoid nested constructions!)** 
1. extend mesh structure with list of structs that define the faces. Use information from GMSH to find global nodes for each face. For each face, compute and store relevant information; 
2. extend mesh structure with list of structs that define the edges; 
3. extend element structure with 4-vector with global index of faces. The loop <i>for fi in element.faces</i> should results the **global** index of the face; 
4. extend double for-loop over elements by four-fold for-loop over faces belonging to elements (using the functions cols maybe?); 

**Post-processing** 
Describe post-processing using VTK files and Paraview. 

<b>To do</b>:
1. introduce [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl) 

wish to obtain single expression for the stiffness matrix per element  

$$
\begin{eqnarray}
A_{\alpha\beta}^{xx,11} & =  & - \int_{P_{\alpha}} \frac{\partial \phi_1({\mathbf r})}{\partial x} \, \frac{\partial}{\partial x} \int_{P_{\beta}} \frac{\phi_1({\mathbf r})}{\|\mathbf{r}' - \mathbf{r} \|} \, d\,P_{\beta} \, d\,P_{\alpha}  \\
& =  & - a_1 \, a_1 \int_{P_{\alpha}} \left[ \int_{P_{\beta}} \frac{\phi_1({\mathbf r})}{\|\mathbf{r}' - \mathbf{r} \|} d\,P_{\beta} \right] d\,P_{\alpha}
\end{eqnarray}
$$ 

## Section 1: Introduction 

When a metalic object is placed inside a magnetic field, a magnetization field is generated in the object. This magnetization opposes the external field. The object is said to be magnetized. The objective of this project is to compute the spatial distribution of the magnetization field inside the object. A coupled system of three partial differential equations for the components of the magnetization vector therefore needs to be solved numerically. These equations contains an integral term that renders this tasks challenging.     

The computation of the magnitization field has numerous practical applications. In the non-destructive testing of metallic object for instance, anomalies in the magnetic field distribution signal cracks or other defects. Other applications can be found in medical imaging, geophysical prospecting and magnetic field sensoring. 

Simulators for the magnization field have a long history of development. Traditional approaches are based on so-called collocation methods. These methods ressemble the finite difference method. The advantage of these methods is that the resulting linear system is straightforward to assemble. In the presence of the integral term in the partial differential equations, however, the resulting coefficient matrix becomes non-symmetric. This renders its large scale deployment cumbersome. 

The objective of this project is to contribute to the development of a novel simulation approach for the magnetization field inside metallic objects. Unlike traditional approaches, our appropach is based on a variational formulation. Our approach thus ressembles a finite element method. The system of integro-differential equations for the components of the magnetization field is discretized by a weighted residual method. The construction of the coefficient matrix of the resulting linear system requires the computation of the six-dimensional singular integrals that represent the source - receiver interactions. These integrals can be computed using an elegant analytical approach based on Euler's integration theorem for homogeneous functions. The linear systems that resuls from our approach are symmetric and positive definite. It is therefore easy to solve. 

The approach we suggest here is expected to render the computation of the magnetization field in realistic applications feasible.    

<b>Assembly of the linear system</b>: The focus of the project is placed on the assembly of the linear system. A computationally efficient procedure to treat the six-dimensional interaction integrals will have to be developed. We will use the Euler Integration Theorem (a generalization of the Green-Gauss divergence theorem) for homogeneous functions.  

<b>Solve of the linear system</b>: Once the linear system is solved, it can be solved by a direct linear solveer for symmetric and positive definite linear system. Such solver are based on a Choleskly decomposition of the coefficient matrix.  

(Insert figures here).

## Section 2: Mathematical Preliminaries 

**Identities** We have that 

$$
\nabla' \cdot \frac{\mathbf{r} - \mathbf{r}'}{\| \mathbf{r} - \mathbf{r}' \|} = \frac{-2}{\| \mathbf{r} - \mathbf{r}' \|}
\text{ and }
\nabla \cdot \frac{\mathbf{r} - \mathbf{r}'}{\| \mathbf{r} - \mathbf{r}' \|} = \frac{2}{\| \mathbf{r} - \mathbf{r}' \|} \, . 
$$

**Application of Euler Integration Theorem for Homogeneous Functions** 

Assume $P_{\alpha}$ to be a tetrahedron bounded by 4 triangular facets $F_{\alpha} \in P_{\alpha}$. Then (add references here)

$$
\int_{P_{\alpha}} \frac{1}{\| {\mathbf r} - {\mathbf r}' \|} \, d\Omega = \frac{1}{2} \sum \int_{F_{\alpha} \in P_{\alpha}} \left[ {\mathbf n}({\mathbf r}) \cdot ({\mathbf r} - {\mathbf r}') \right] \frac{1}{\| {\mathbf r} - {\mathbf r}' \|} \, dS
$$

and 

$$
\int_{P_{\alpha}} \frac{{\mathbf r} - {\mathbf r}'}{\| {\mathbf r} - {\mathbf r}' \|} \, d\Omega = \frac{1}{3} \sum \int_{F_{\alpha} \in P_{\alpha}} \left[ {\mathbf n}({\mathbf r}) \cdot ({\mathbf r} - {\mathbf r}') \right] \frac{{\mathbf r} - {\mathbf r}'}{\| {\mathbf r} - {\mathbf r}' \|} \, dS 
$$

**Integration by parts** for grad-div equations. Assume $\mathbf{V}(\mathbf{r})$ and $\mathbf{W}(\mathbf{r})$ to be two vector fields on $\Omega \subset \mathbb{R}^3$ bounded by a surface $\partial \Omega$ with outward normal ${\mathbf n}$. Then 

$$
\text{div} \left[ \mathbf{V} \,  \text{div} \mathbf{W} \right] = 
\text{div} \mathbf{V} \, \text{div} \mathbf{W} + \mathbf{V} \cdot \text{grad} \left( \text{div} \mathbf{W} \right) \, . 
$$

Then by integration over $\Omega$ and applying the Green-Gauss divergence theorem, we obtain that 

$$
\int_{\Omega} \text{grad} \left( \text{div} \mathbf{W} \right) \cdot \mathbf{V} \, d\Omega = \int_{\partial \Omega} \left[ \mathbf{V} \,  \text{div} \mathbf{W} \right] \cdot \mathbf{n} \, dS - \int_{\Omega} \text{div} \mathbf{V} \, \text{div} \mathbf{W} \, d\Omega \, .  
$$

## Section 3: Project Description 

### Section 1.3: Problem Formulation 

<b>Computational Domain</b> We assume that the computational domain $\Omega$ is a cube with lenght $L$, height $H$ and depth $D$ alligned with the coordinate axes. Then $0 \leq x \leq L$, $0 \leq y \leq H$ and $0 \leq z \leq D$ and $\Omega = [0,L] \times [0,H] \times [0,D]$. Typical values are $L = H = 1 \, \text{m}$ and $0.01 \, \text{m} \leq D \leq 0.1 \, \text{m}$. In thick plate, the magnetization will varry with the $z$-coordinates. In thin plates, this variation is expected to be small. 

We will distinguish between the interior and the boundary of $\Omega$. We therefore assume that $\Omega$ is open (the boundary of $\Omega$ does not belong to $\Omega$). We assume that $\partial \Omega$ (the union of the four quadrilateral facets) denotes the boundary of $\Omega$. We assume that $\overline{\Omega} = \Omega \cup \partial \Omega$ (i.e. the union of the volume and its surface) denotes the closure of $\Omega$. We assume that the exterior to $\Omega$ is the volume ${\mathbb R}^3 \setminus \overline{\Omega}$. 

We assume that $\mathbf{r} = (x,y,z)$ denotes the position vector of a observer (destination) inside of $\overline{\Omega}$ (thus possibly on the boundary). We assume that $\mathbf{r}' = (x',y',z')$ denotes the position vector of a source inside of $\overline{\Omega}$. Introduce notation for gradient and divergence wrt. $\mathbf{r}$ and 
$\mathbf{r}'$. 

<b>Material Properties</b> Assume that $\Omega$ has a [magnetic susceptibility](https://en.wikipedia.org/wiki/Magnetic_susceptibility) denoted by $\chi_{mag}$ and a relative [magnetic permeability](https://en.wikipedia.org/wiki/Permeability_(electromagnetism)) denoted by $\mu_r$. Then $\chi_{mag} = \mu_r -1$. When the volume $\Omega$ is placed inside an external magnetic field, it will magnetize. The volume $\Omega$ thus mimmics a metalic plate. Here it will be suffucient to assume that the plate is homogeneous and that therefore $\chi_{mag}$ is constant. Assume that $10 \leq \chi_{mag} \leq 1000$ (dimensionless, value for common steel types). Non-homogeneous plate can be modeled assuming that $\chi_{mag}$ is piecewise constant. Non-linear magnetization effects are excluded in this project.   

<b>Problem Description</b> To describe the problem, we introduce three [magnetic fields](https://en.wikipedia.org/wiki/Magnetic_field) (in units Ampere per meter). 

The external magnetic field will be denoted by ${\mathbf H}_{ext}$. This field causes the magnetization to be generated in $\Omega$. This field is assumed to be constant and equal in the inside and outside of $\Omega$. Then $\mathbf{H}_{ext}$ is a vector with three components . This can be expressed as $\mathbf{H}_{ext} = \left( H_{ext,x}, H_{ext,y}, H_{ext,z}\right)$ or as $ \mathbf{H}_{ext} = H_{ext,x} \mathbf{i} + H_{ext,y} \mathbf{j} + H_{ext,z} \mathbf{k}$. The orientation of $\Omega$ with respect to $\mathbf{H}_{ext}$ matters. We distinguish an in-plane and a perpendicular plane configuration. In an in-plane configuration, the plate $\Omega$ is alligned with the external field, e.g., with the $x$-axis. We can then write that $\mathbf{H}_{ext} = \left( H_{ext,x}, 0,0\right)$ with $H_{ext,x}$ a constant. In a perpendicular plane configuration, the plate $\Omega$ is oriented perpendicular to $\mathbf{H}_{ext}$. We can then write that $\mathbf{H}_{ext} = \left( 0, 0, H_{ext,z} \right)$ with $H_{ext,z}$ a constant. (Provide a representative value for the strength of the magnetic field. The earth magnetic field for instance is on average $50,000$ nanoTesla). 

The induced magnetic field will be denoted by ${\mathbf H}_{M}(\mathbf{r})$. This field is caused by the presence of a non-zero magnetization ${\mathbf M}(\mathbf{r})$ inside of $\Omega$. The relation between ${\mathbf M}(\mathbf{r})$ and ${\mathbf H}_{M}(\mathbf{r})$ can be expressed through the magnetic scalar potential $\phi_M(\mathbf{r})$ as 

$$
{\mathbf H}_{M}(\mathbf{r}) = - \text{grad}_{\mathbf{r}} \phi_M(\mathbf{r}) = - \nabla_{\mathbf{r}} \phi_M(\mathbf{r}) 
$$

(minus sign by convention). The scalar potential for ${\mathbf M}(\mathbf{r})$ inside of $\Omega$ (volumetric density $\rho_M = - \nabla_{\mathbf{r'}} \cdot \mathbf{M}(\mathbf{r}')$) and ${\mathbf M_{\partial}}(\mathbf{r})$ on $\partial \Omega$ (surface density $\rho_S = \mathbf{M}_{\partial}(\mathbf{r}') \cdot \mathbf{n}(\mathbf{r'})$) is given by 

$$
\phi_M(\mathbf{r}) = \frac{1}{4 \pi}
\int_{\Omega} \frac{- \nabla_{\mathbf{r'}} \cdot \mathbf{M}(\mathbf{r}')}{ \|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' + \frac{1}{4 \pi} \int_{\partial \Omega} \frac{ \mathbf{M}_{\partial}(\mathbf{r}') \cdot \mathbf{n}(\mathbf{r'}) }{\|\mathbf{r}' - \mathbf{r} \|} \, dS' \, . 
$$

(Who is $\mathbf{M}_{\partial}(\mathbf{r})$)? The total magnetic field will be denoted by ${\mathbf H}_{tot}(\mathbf{r})$. This field is the sum of the external and the induced magnetic field. We thus have that ${\mathbf H}_{ext} + {\mathbf H}_{M} = {\mathbf H}_{tot}$. The total field is equal to the external field in the exterior of $\overline{\Omega}$. The total field on $\overline{\Omega}$ is related to the magnetization via the constitutive equation ${\mathbf H}_{tot} = {\mathbf M}/\chi_{mag}$.

<b>Mathematical Model</b> We wish to compute $\mathbf{M}(\mathbf{r})$ by solving a partial differential equation in which the external field ${\mathbf H}_{ext}$ acts as a source term (i.e., ${\mathbf H}_{ext}$ appears in the right-hand side of the equation). This differential equation can be written as  

$$
{\cal D} \left[ \mathbf{M}(\mathbf{r}) \right] = {\mathbf H}_{ext} \text{ on } \overline{\Omega}  
$$

where ${\cal D}$ is the [differential operator](https://en.wikipedia.org/wiki/Differential_operator) acting on $\mathbf{M}(\mathbf{r})$ as 

$$
\begin{eqnarray}
{\cal D}\left[ \mathbf{M}(\mathbf{r}) \right] & = & {\mathbf H}_{tot}(\mathbf{r}) - {\mathbf H}_{M}(\mathbf{r}) \nonumber \\ 
& = & \frac{1}{\chi_{mag}} \mathbf{M}(\mathbf{r}) +
\text{grad}_{\mathbf{r}} \left\{  
\frac{1}{4 \pi}  \int_{\Omega} \frac{- \text{div}_{\mathbf{r}} \mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' + \frac{1}{4 \pi}  \int_{\partial \Omega} \frac{ \mathbf{M}_{\partial}(\mathbf{r}') \cdot \mathbf{n}(\mathbf{r'}) }{\|\mathbf{r}' - \mathbf{r} \|} \, dS' \right\} \nonumber \\
& = & \frac{1}{\chi_{mag}} \mathbf{M}(\mathbf{r}) + 
\nabla_{\mathbf{r}} \left\{ 
\frac{1}{4 \pi} \int_{\Omega} \frac{- \nabla_{\mathbf{r}} \cdot \mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' + \frac{1}{4 \pi}  \int_{\partial \Omega} \frac{ \mathbf{M}_{\partial}(\mathbf{r}') \cdot \mathbf{n}(\mathbf{r'}) }{\|\mathbf{r}' - \mathbf{r} \|} \, dS' \right\} \, . \nonumber 
\end{eqnarray}
$$

The differential operator ${\cal D}$ contains a volume integral term and a surface integral term. The differential equation that determines $\mathbf{M}(\mathbf{r})$ is therefore an [integro-differential equation](https://en.wikipedia.org/wiki/Integro-differential_equation). The above integro-differential equation is a system of three coupled equations for the three components of the magnetization $\mathbf{M}(\mathbf{r})$. This equation is solved on $\overline{\Omega}$ only. The magnetization $\mathbf{M}(\mathbf{r})$ remains zero outside of $\overline{\Omega}$. The kernel in the second term of ${\cal D}$ can be written as 

$$
{\cal K}(\mathbf{r},\mathbf{r}') = \frac{1}{4 \pi}  \frac{1}{\|\mathbf{r}' - \mathbf{r} \|} \, . 
$$

This kernel depends on the distance $\|\mathbf{r}' - \mathbf{r} \|$ only. Small values of this distance matter most. It is singular for $\mathbf{r}' = \mathbf{r}$. The convolution of $\mathbf{M}(\mathbf{r})$ and ${\cal K}(\mathbf{r},\mathbf{r}')$ can be written as 

$$
{\cal K}(\mathbf{r},\mathbf{r}') * \mathbf{M}(\mathbf{r}) = \frac{1}{4 \pi} \int_{\Omega} \frac{- \nabla_{\mathbf{r}} \cdot \mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, . 
$$

For future reference, we will split ${\cal D}\left[ \mathbf{M}(\mathbf{r}) \right]$ in three terms as ${\cal D}\left[ \mathbf{M}(\mathbf{r}) \right] = {\cal D}_1\left[ \mathbf{M}(\mathbf{r}) \right] + {\cal D}_2 \left[ \mathbf{M}(\mathbf{r}) \right] + {\cal B} \left[ \mathbf{M}(\mathbf{r}) \right]$, where ${\cal D}_1$ is merely a scaling  

$$
{\cal D}_1\left[ \mathbf{M}(\mathbf{r}) \right] = \frac{1}{\chi_{mag}} \mathbf{M}(\mathbf{r})
$$

where ${\cal D}_2$ is the integro-differential part of the equation

$$
{\cal D}_2\left[ \mathbf{M}(\mathbf{r}) \right] = - 
\nabla_{\mathbf{r}} \, {\cal K}(\mathbf{r},\mathbf{r}') * \mathbf{M}(\mathbf{r}) \, .
$$

and where ${\cal B}$ is the boundary integral term

$$
{\cal B} \left[ \mathbf{M}(\mathbf{r}) \right] = \nabla_{\mathbf{r}} \, \frac{1}{4 \pi}  \int_{\partial \Omega} \frac{ \mathbf{M}_{\partial}(\mathbf{r}') \cdot \mathbf{n}(\mathbf{r'}) }{\|\mathbf{r}' - \mathbf{r} \|} \, dS' \, . 
$$

We also write the integro-differential equation for $\mathbf{M}(\mathbf{r})$ in explicit form. The first of the three scalar equations is 

$$
\frac{1}{\chi_{mag}} M_x(\mathbf{r}) +  
\frac{\partial}{\partial x} \, \frac{1}{4 \pi}
\int_{\Omega} \frac{- \nabla_{\mathbf{r}} \cdot \mathbf{M}(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' + 
\frac{\partial}{\partial x} \frac{1}{4 \pi}  \int_{\partial \Omega} \frac{ \mathbf{M}_{\partial}(\mathbf{r}') \cdot \mathbf{n}(\mathbf{r'}) }{\|\mathbf{r}' - \mathbf{r} \|} \, dS' = H_{ext,x} \, , 
$$

and simarly for $M_y(\mathbf{r})$ and $M_z(\mathbf{r})$.   

### Section 2.3: Mesh Generation and Shape Functions 

We wish to compute the magnetization vector $\mathbf{M}(\mathbf{r})$ inside the metallic plate $\Omega$. We therefore generate a mesh on $\Omega$. 

<b>Mesh Generation</b> We denote the mesh by $\Omega^h$. We assume that the mesh consists of tetrahedral elements only. We denote the number of nodes, edges, facets and elements of $\Omega^h$ by $N_n$, $N_{ed}$, $N_f$ and $N_e$, respectively. Assume the nodes and the elements to be denoted by ${\mathbf x}_i$ for $1 \leq i \leq N_n$ and by $P_{\alpha}$ for $1 \leq \alpha \leq N_e$, respectively. The union of all elements $P_{\alpha}$ forms the entire domain of computation. That is, we have that $\cup P_{\alpha} | 1 \leq \alpha \leq N_e = \Omega$. We will use this elementary fact in the computation of the matrix $\underline{\underline{A}}$ and the right-hand side vector $\mathbf{b}$.    

We assume that the information to decompose an element $P_{\alpha}$ into a set of facets, a facet into a set of edges and an edge into a set of points to be available (representation of the mesh $\Omega^h$ as a directed a-cyclic graph (DAG) with corresponding operations to find parent nodes and child nodes. See e.g. [wikipedia entry on polygon mesh](https://en.wikipedia.org/wiki/Polygon_mesh).

<b>Shape Functions</b> The mesh $\Omega^h$ allows to define a set of nodal P1 linear Lagrange basis functions $\{ \phi_i(\mathbf{r}) | 1 \leq i \leq N _n \}$ (see classical finite elements method). The functions $\phi_i(\mathbf{r})$ are scalar functions. To be able to cast a vector-valued partial differential equation in weak or variational form, we define $3 \, N_n$ vector-valued shape functions $\boldsymbol{\phi}_{3\,i-2} = \left( \phi_i, 0,0\right)$ (weighting function for the PDE for $M_x({\mathbf r})$), 
$\boldsymbol{\phi}_{3i-1} = \left( 0, \phi_i, 0\right)$ (weighting function for $M_y({\mathbf r})$) and 
$\boldsymbol{\phi}_{3i} = \left( 0, 0, \phi_i \right)$ (weighting function for $M_z({\mathbf r})$ for $1 \leq i \leq N_n$.  

More information on the mesh generation is provided in the [notebook](./mom_3d_plate_gmsh.ipynb). 

### Section 3.3: Galerkin Method    

Here we describe the [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method) that allows to convert the system of integral partial differential-equations for $\mathbf{M}(\mathbf{r})$ into a weak or variational formulation. This variational formulation will allow a spatial discretization. See also example of the weak formulation of the Poisson equation as [wiki on weak formulation](https://en.wikipedia.org/wiki/Weak_formulation)). (Requires example of weak form of a c oupled system of differential equations.) 

Assume $g(\mathbf{r})$ and $h(\mathbf{r})$ to be scalar functions on $\Omega$ such that 

$$
g(\mathbf{r}) = \int_{\Omega} \frac{h(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, . 
$$

By decomposing $\Omega$ into a set of volumes $P_{\beta}$, we find that

$$
g(\mathbf{r}) = \int_{\Omega} \frac{h(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' = \sum_{\beta=1}^{N_e}  \int_{P_{\beta}} \frac{h(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, . 
$$

Then by decomposing $\Omega$ again into a set of volumes $P_{\alpha}$, the volume average of $g(\mathbf{r})$ can be written as  

$$
\begin{eqnarray}
\int_{\Omega} g(\mathbf{r}) \, d\Omega = 
\sum_{\alpha=1}^{N_e}  \int_{P_{\alpha}} g(\mathbf{r}) \, d\Omega 
& = & \sum_{\alpha,\beta=1}^{N_e}  \int_{P_{\alpha}} \int_{P_{\beta}} \frac{h(\mathbf{r}')}{\|\mathbf{r}' - \mathbf{r} \|} \, d\Omega' \, d\Omega \, . 
\end{eqnarray}
$$  

(Extend to weighted averaging of $g(\mathbf{r})$, Fourier analysius, expansion in sets of orthogonal functions.) 

**Expansion of Magnetization Components** The numerical approximation of the magnetization components can be expanded in as 

\begin{eqnarray}
M_x(\mathbf{r}) = \sum_{i=1}^{N_n} c_i \, \phi_i(\mathbf{r}) \\ 
M_y(\mathbf{r}) = \sum_{i=1}^{N_n} d_i \, \phi_i(\mathbf{r}) \\
M_z(\mathbf{r}) = \sum_{i=1}^{N_n} e_i \, \phi_i(\mathbf{r}) 
\end{eqnarray}

Note that components $M_x(\mathbf{r})$, $M_y(\mathbf{r})$ and $M_z(\mathbf{r})$ are expanded in the same basis with proper expansion coefficients. In case of linear Lagrange shape functions on tetrahedra, the element $P_{\alpha}$ has 4 local degrees of freediom per components, thus in total $3*4=12$ local degrees of freedom. Mapping from local to global degrees of freedom. The expansion coefficients $c_i$, $d_i$ and $e_i$ can stacked into a global vector $\mathbf{u}$ of dimension $3 \, N_n$ where $\mathbf{u}^T = \begin{pmatrix} \mathbf{c}^T \mathbf{d}^T \mathbf{e}^T \end{pmatrix}$. This global vector is the solution of the Galerkin-MoM linear system 

$$
\underline{\underline{A}} \, \mathbf{u} = \mathbf{b}
$$

where the matrix $\underline{\underline{A}} = \underline{\underline{A}}^{(1)} + \underline{\underline{A}}^{(2)}$ is a dense $3 \, N_n$-by-$3 \, N_n$ matrix consisting of 
mass matrix part $\underline{\underline{A}}^{(1)}$ and a stiffness part $\underline{\underline{A}}^{(2)}$. The mass matrix is the spatial discretization of the differential operator ${\cal D}_1 [\mathbf{M}(\mathbf{r})]$ and is very similar to a mass matrix in a classical finite element formulation. The stiffness matrix is the spatial discretization of the differential operator ${\cal D}_2 [\mathbf{M}(\mathbf{r})]$. Assembling this matrix is the main challenge of this project. 

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

resulting in $3 \, N_n$ equations for the $3 \, N_n$ components of the vector of expansion coefficients $\mathbf{u}$. The mass matrix and load vector can be assembled by a loop over elements as in classical Galerlin FEM. The stiffness matrix requires dedicated attention due to the integral term.

### Section 4.3: Element-by-element Assembly of Stiffness Matrix

**Definition of Shape Functions on a Single Tetrahedral Element** The tetrahedral mesh element $P_{\alpha}$ has $4$ nodes. Assume these nodes to be denote ${\mathbf x}_1$, ${\mathbf x}_2$, ${\mathbf x}_3$ and ${\mathbf x}_4$ (in local numbering on this single element). The linear Lagrange shape functions on this element can be expressed as $\phi_k({\mathbf r}) = a_k x + b_k y + c_k y + d_k$ for $1 \leq k \leq 4$, where $a_k$, $b_k$, $c_k$ and $d_k$ are coefficients. These coefficients are choosen such that the shape functions satisfy the constraint that $\phi_k({\mathbf x}_{\ell}) = \delta_{k\ell}$ for $1 \leq k,\ell \leq 4$. Linear system for these coefficients. 

**Divergence of Vector Shape Function** The divergence of the vector-valued shape functions can be expressed as 

$$
\nabla_{\mathbf{r}} \cdot \boldsymbol{\phi}_k = \frac{\partial}{\partial x} (a_k x + b_k y + c_k y + d_k) = a_k \text{ if } \mod(k,3) = 1 \\
\nabla_{\mathbf{r}} \cdot \boldsymbol{\phi}_k = \frac{\partial}{\partial y} (a_k x + b_k y + c_k y + d_k)= b_k \text{ if } \mod(k,3) = 2 \\
\nabla_{\mathbf{r}} \cdot \boldsymbol{\phi}_k = \frac{\partial}{\partial z} (a_k x + b_k y + c_k y + d_k)= c_k \text{ if } \mod(k,3) = 0
$$

**Bilinear Form** By writing $d\Omega = \sum_{\alpha} dP_{\alpha}$ and $d\Omega' = \sum_{\beta} dP_{\beta}$, we arrive at 

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

where (row determined by the expansion of the magnetization, columns determined by the test function) (assuming tetrahedra and $N_e^{\alpha}=4$). 

$$
\underline{\underline{A}}_{loc,\alpha\beta}^{(2)} = 
\begin{pmatrix} 
{\mathbf a}^{\alpha} ({\mathbf a}^{\alpha\beta})^T & {\mathbf a}^{\alpha} ({\mathbf b}^{\alpha\beta})^T & {\mathbf a}^{\alpha} ({\mathbf c}^{\alpha\beta})^T \\ 
{\mathbf b}^{\alpha} ({\mathbf a}^{\alpha\beta})^T & {\mathbf b}^{\alpha} ({\mathbf b}^{\alpha\beta})^T & {\mathbf b}^{\alpha} ({\mathbf c}^{\alpha\beta})^T \\
{\mathbf c}^{\alpha} ({\mathbf a}^{\alpha\beta})^T & {\mathbf c}^{\alpha} ({\mathbf b}^{\alpha\beta})^T & {\mathbf c}^{\alpha} ({\mathbf c}^{\alpha\beta})^T 
\end{pmatrix}
$$

More explicitly we have that 

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

## Section 4: Possible Project Roadmaps  

### Section 1.4: Problem Formulation 

<b>Physical Problem</b> Study problem formulation. Study how the magnetization vector $\mathbf{M}(\mathbf{r})$ change with  of the plate, dimensions of the plate, and the direction of the external magnetic field. Use to this end literature or a reference code.

### Section 2.4: Solution Method  

Explore weak formulation, weighted residual method, Galerkin approximation and linear system formulation for a sequence of problems of increasing complexity.  

1. Poisson equation;
2. one-dimensional Fredholm integro-differential equation. A model problem can be constructed using Sympy;
3. see grad-div equation without kernel;

### Section 3.4: Assembly of the Mass Matrix and Load Vector 

1. loop over elements in the mesh to assemble the mass matrix and the lopad vector. See notebook [mom_hcubature](./mom_hcubature.ipynb); 

### Section 4.4: Numerical Assembly of the Stiffness Matrix 

1. use of quadrature implemented in [hcubature.jl](https://github.com/JuliaMath/HCubature.jl) for integration in 1D (possibly singular, look into number of function evaluatiohs), 2D (reference triangle, coordinate transformation, general triangle), 3D (reference tetrahedra, coordinate transformation, general tetrahedra), 4D (by calling hcubature for 2D twice) and 6D (by calling hcubature for 2D twice). See notebook [mom_hcubature](./mom_hcubature.ipynb);
2. formulate and compute 6D integrals for tetra/tetra interaction assuming no parallel facets. $P_{\alpha}$ with nodes ${\mathbf r}_1 = (0,0,0)$, ${\mathbf r}_2 = (1,0,0)$, ${\mathbf r}_3 = (0,1,0)$ and ${\mathbf r}_4 = (0,0,1)$. 
$P_{\beta}$ with nodes ${\mathbf r}_1 = (0,0,0)$, ${\mathbf r}_2 = (1,0,0)$, ${\mathbf r}_3 = (0.5,1,0)$ and ${\mathbf r}_4 = (0.5,0.5,1)$. Verify that no two faces are parallel to each other;
3. extend previous case to parallel faces. $P_{\alpha}$ reference triangle. $P_{\beta}$ shift of $P_{\alpha}$;
4. apply Euler theorem for homogeneous function and compute resulting 4D integrals numerically; 
5. loop over elements over the mesh to compute the load vector. Assume constant external field. Provide more details (expression of the local vector per element) here. See notebook [mom_3d_plate_gmsh](./mom_3d_plate_gmsh.ipynb); 
6. loop over elements over the mesh to compute the mass matrix. Provide more details (expression of the local matrix per element) here;
7. use of hcubature to compute the stiffness matrix elements. Possibly regularization in the kernel;

### Section 5.4: Analytical Assembly of the Stiffness Matrix 

1. implement point-edge interactions (cfr. seperate notebook); 
2. implement edge-edge interactions  (cfr. seperate notebook); 
2. implement point-facet interactions  (cfr. seperate notebook);

## Section 5: The Julia Programming Language

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

### References of the Method of Moments 

(Provide more guidance on what and how to read.) 

1. R. F. Harrington, <i>Field Computation by Moment Methods</i>, Wiley-IEEE Press, 1993. 
2. Morandi, A., Fabbri, M., & Ribani, P. L. (2010). <i>A modified formulation of the volume integral equations method for 3-D magnetostatics</i>. IEEE transactions on magnetics, 46(11), 3848-3859. And similar references;
3. Le-Duc, T., Meunier, G., Chadebec, O., Guichon, J. M., & Bastos, J. P. A. (2013). <i> General integral formulation for the 3D thin shell modeling</i>. IEEE transactions on magnetics, 49(5), 1989-1992. And similar references; 

### References on the weighted residual method, the variational formulation and the finite element method 

1. [Introduction to Numerical Methods for Variational Problems](https://link.springer.com/book/10.1007/978-3-030-23788-2) by Hans Petter Langtangen and Kent-Andre Mardal. The [book](https://link.springer.com/book/10.1007/978-3-030-23788-2) is freely available; 
2. [Wolfgang Bangerth's video lectures](https://www.math.colostate.edu/~bangerth/videos.html); 
3. [wiki Finite Element Method](https://en.wikipedia.org/wiki/Finite_element_method): Section 3 for the weak form and Section 4 for the finite element discretization;  
4. [Comsol Multiphysics Finite Element Method](https://www.comsol.com/multiphysics/finite-element-method): more information and illustrations; 
5. [Comsol Multiphysics Brief Introduction to the Weak Form](https://www.comsol.com/blogs/brief-introduction-weak-form): good introduction to a theoretical concept that provides a basis for the finite element method; 
6. [Ferrite Introduction to FEM](https://ferrite-fem.github.io/Ferrite.jl/stable/topics/fe_intro/)

**Reference Tetrahedron** Assume the reference tetrahedron to have the nodes ${\mathbf x}_1 = (0,0,0)$, ${\mathbf x}_2 = (1,0,0)$, ${\mathbf x}_3 = (0,1,0)$ and ${\mathbf x}_4 = (0,0,1)$. Then the shape function are defined as 

$$
\phi_1({\mathbf r}) = 1 - x - y - z \\
\phi_2({\mathbf r}) = x \\
\phi_3({\mathbf r}) = y \\
\phi_4({\mathbf r}) = z 
$$

**Two-Tetrahedra Test Case** For future reference, we define a test case involving two tetrahedra $P_{\alpha}$ and $P_{\beta}$. Assume $P_{\alpha}$ to be the reference element. 
Assume $P_{\beta}$ to be the reference element shifted along the $x$-axis by a distance of $x_0$. Then 

$$
\phi^{\alpha}_1({\mathbf r}) = 1 - x - y - z \\
\phi^{\alpha}_2({\mathbf r}) = x \\
\phi^{\alpha}_3({\mathbf r}) = y \\
\phi^{\alpha}_4({\mathbf r}) = z 
$$

$$
\phi^{\beta}_1({\mathbf r}) = 1 - (x-x_0) - y - z \\
\phi^{\beta}_2({\mathbf r}) = x - x_0 \\
\phi^{\beta}_3({\mathbf r}) = y \\
\phi^{\beta}_4({\mathbf r}) = z 
$$

$$
{\mathbf a}^{\alpha} = \begin{pmatrix} -1 \\ 1 \\ 0 \\ 0 \end{pmatrix}
\hspace{.5cm}
{\mathbf b}^{\alpha} = \begin{pmatrix} -1 \\ 0 \\ 1 \\ 0 \end{pmatrix}
\hspace{.5cm}
{\mathbf c}^{\alpha} = \begin{pmatrix} -1 \\ 0 \\ 0 \\ 1 \end{pmatrix}
$$


To solve a Fredholm integral equation in SymPy, you cannot use a single "built-in" solver function. Instead, you represent the equation symbolically and use techniques like the Method of Undetermined Coefficients to solve for the unknown function.A classic Fredholm equation of the second kind looks like this:
$$
y(x)=f(x)+\lambda \int _{a}^{b}K(x,t)y(t)\,dt
$$
Here is a complete example solving the Fredholm integral equation:
$$
y(x)=x+\int _{0}^{1}x\cdot t\cdot y(t)\,dt
$$
Step-by-Step Python Code

Since the kernel $K(x, t) = x \cdot t$ is separable, the integral $\int_{0}^{1} t \cdot y(t) \, dt$ evaluates to a constant, let's call it \(C\). This allows us to solve it algebraically.

When you run this script below, SymPy performs the following calculations: It evaluates the definite integral $\int_{0}^{1} t(t + C t) \, dt = \int_{0}^{1} (1+C)t^2 \, dt = \frac{1+C}{3}$ It sets up the linear equation $C = \frac{1+C}{3}$ and solves for $C$, yielding $C = \frac{1}{2}$. It substitutes $C$ back into $y(x) = x + Cx$ to yield the exact analytical solution:$\mathbf{y(x)=}\frac{\mathbf{3x}}{\mathbf{2}}$. 


```julia
from sympy import symbols, Integral, solve, Eq

# 1. Define symbols
x, t, C = symbols('x t C')

# 2. Assume the form of y(x) based on the structure of the equation
# Since y(x) = x + x * Integral(t * y(t)), y(x) must be of the form: x + C * x
y_x = x + C * x

# 3. Substitute this assumed form into the integral part to find C
# The constant C is defined as the integral of t * y(t) from 0 to 1
y_t = y_x.subs(x, t)  # Replace x with t for the integration variable
integral_expr = Integral(t * y_t, (t, 0, 1))

# 4. Set up the compatibility equation: C = Integral(t * y(t))
compat_eq = Eq(C, integral_expr.doit())

# 5. Solve for the unknown constant C
C_value = solve(compat_eq, C)[0]

# 6. Substitute C back into our assumed y(x) to get the final solution
final_solution = y_x.subs(C, C_value)

print(f"The constant C is: {C_value}")
print(f"The solution y(x) is: {final_solution}")

```
