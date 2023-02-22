# Finite Element Methods for Electrical Engineering Applications (EE4375)

## Section 1:/ Introduction

This course consists of three blocks. 

### First Block 
In the first block we discuss the finite difference method to solve the Poisson equation in one and two spatial dimensions. We restrict ourselves to the interval and the square discretized by a uniform mesh. The Poisson equation models a diffusion process and typically appears in the computation of electrostatic and magneto-static fields. Our motivation is to introduce initial problem formulations and simplified numerical solution methods. 

### Second Block 
In the second block we discuss the variational formulation of the Poisson equation in one spatial dimension. This variational formulation allows to subsequently solve the problem by a finite element method using non-uniform meshes. We show that the versatility of the finite element method in solving problems that are more complex than in the first block. 

### Third Block 
In the third and last block we discuss the finite element solution of the Poisson equation in two spatial dimension. We discuss the mesh generation using triangle and the construction of the discrete problem using a loop over all elements. We illustrate the method in the computation of magnetostatic fields in transformers and electrical machines.  

### Assessment 
Each block will be concluded with a graded homework assignment. The three homework assignments have to be defended during the oral exam. 

### Learning Activities 
The course consists of two hours of lectures and four hours of computer lab sessions each week. We will use the programing language julia (see julialang.org) and Jupiter notebooks to develop the assignments. We will use the gmsh mesh generation software. You are kindly requested to install Julia, the Jupiter notebook system and the Julia interface to gmsh (Gmsh.jl) on your laptop. You are most welcome to give a look at the references below. 

### References for First Block
1. [wiki on the Laplace equation](https://en.wikipedia.org/wiki/Laplace's_equation)
2. [wiki on the Poisson equation](https://en.wikipedia.org/wiki/Poisson%27s_equation): including an example on electrostatics
3. [wiki on boundary value problem](https://en.wikipedia.org/wiki/Boundary_value_problem)
4. [wiki on electrostatics](https://en.wikipedia.org/wiki/Electrostatics)
5. [wiki on magnetostatics](https://en.wikipedia.org/wiki/Magnetostatics)
6. [wiki on discrete Poisson matrix](https://en.wikipedia.org/wiki/Discrete_Poisson_equation)  
7. [wiki on finite difference method](https://en.wikipedia.org/wiki/Finite_difference_method) 

### References for Second Block
1. [Introduction to Numerical Methods for Variational Problems](https://link.springer.com/book/10.1007/978-3-030-23788-2) by Hans Petter Langtangen and Kent-Andre Mardal. The [book](https://link.springer.com/book/10.1007/978-3-030-23788-2) is freely available; 
2. [Wolfgang Bangerth's video lectures](https://www.math.colostate.edu/~bangerth/videos.html); 
3. [wiki Finite Element Method](https://en.wikipedia.org/wiki/Finite_element_method): Section 3 for the weak form and Section 4 for the finite element discretization;  
4. [Comsol Multiphysics Finite Element Method](https://www.comsol.com/multiphysics/finite-element-method): more information and illustrations; 
5. [Comsol Multiphysics Brief Introduction to the Weak Form](https://www.comsol.com/blogs/brief-introduction-weak-form): good introduction to a theoretical concept that provides a basis for the finite element method; 
6. [Sphinx Finite Element Method](http://hplgit.github.io/INF5620/doc/pub/sphinx-fem/): reference for implementation;

## Section 2:/ Lectures 
The course slides are available in the slides directory. 

## Section 3:/ Lab Sessions
1. [first lab session](lab-sessions/1-lab-session.ipynb) and [small GMSH manual](extended-lab-sessions/gmsh/Mesh-Generation-using-Gmsh.ipynb)
2. [second lab session](lab-sessions/2-lab-session.ipynb)
3. [third lab session](lab-sessions/3-lab-session.ipynb)
4. [fourth lab session](lab-sessions/4-lab-session.ipynb)
5. [fifth lab session](lab-sessions/5-lab-session.ipynb)
6. [sixth lab session](lab-sessions/6-lab-session.ipynb)
7. [seventh lab session](lab-sessions/7-lab-session.ipynb)


```julia

```
