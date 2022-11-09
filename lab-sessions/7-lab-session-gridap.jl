{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "670e494f",
   "metadata": {},
   "source": [
    "# Scalar diffusion and convection-diffusion problem on channel geometry with uniform mesh using Gridap  "
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a5bba0aa",
   "metadata": {},
   "source": [
    "## Import libraries "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 104,
   "id": "ed7f6390",
   "metadata": {},
   "outputs": [],
   "source": [
    "using Gridap"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "f1937cbb",
   "metadata": {},
   "source": [
    "## Section 1:/ Diffusion Problem "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 105,
   "id": "cbece9d8",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "CartesianDiscreteModel()"
      ]
     },
     "execution_count": 105,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "domain = (0,3,0,1)\n",
    "partition = (30,10)\n",
    "model = CartesianDiscreteModel(domain,partition)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 113,
   "id": "e932c353",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "TrialFESpace()"
      ]
     },
     "execution_count": 113,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# Define Lagrangian reference element\n",
    "order = 1\n",
    "reffe = ReferenceFE(lagrangian,Float64,order)\n",
    "Vh = FESpace(model, reffe; conformity=:H1, dirichlet_tags = [7,8])\n",
    "Uh = TrialFESpace(Vh, [0, 1]) # Left boundary: u = 0, right boundary: u = 1"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 114,
   "id": "0dcfa5cc",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "Measure()"
      ]
     },
     "execution_count": 114,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "degree = 2\n",
    "Ω = Triangulation(model)\n",
    "dΩ = Measure(Ω,degree)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 115,
   "id": "31094c6d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "b (generic function with 1 method)"
      ]
     },
     "execution_count": 115,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# working version \n",
    "ff(x) = 0\n",
    "a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ\n",
    "b(v) = ∫( v*ff )*dΩ"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 116,
   "id": "fccd04bd",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "AffineFEOperator()"
      ]
     },
     "execution_count": 116,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "op = AffineFEOperator(a,b,Uh,Vh)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 117,
   "id": "0a95809a",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "LinearFESolver()"
      ]
     },
     "execution_count": 117,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "ls = LUSolver()\n",
    "solver = LinearFESolver(ls)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 118,
   "id": "bed6fad5",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "SingleFieldFEFunction():\n",
       " num_cells: 300\n",
       " DomainStyle: ReferenceDomain()\n",
       " Triangulation: BodyFittedTriangulation()\n",
       " Triangulation id: 1733168475160219817"
      ]
     },
     "execution_count": 118,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "uh = solve(solver,op)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 119,
   "id": "8a75152f",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "([\"results.vtu\"],)"
      ]
     },
     "execution_count": 119,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "writevtk(Ω,\"results\",cellfields=[\"uh\"=>uh])"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "ac082bdb",
   "metadata": {},
   "source": [
    "## Section 2:/ Convection-Diffusion Problem "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 215,
   "id": "528d1d2d",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "CartesianDiscreteModel()"
      ]
     },
     "execution_count": 215,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "domain = (0,3,0,1)\n",
    "partition = (30,10)\n",
    "model = CartesianDiscreteModel(domain,partition)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 216,
   "id": "df81863e",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "TrialFESpace()"
      ]
     },
     "execution_count": 216,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# Define Lagrangian reference element\n",
    "order = 1\n",
    "reffe = ReferenceFE(lagrangian,Float64,order)\n",
    "Vh = FESpace(model, reffe; conformity=:H1, dirichlet_tags = [7,8])\n",
    "Uh = TrialFESpace(Vh, [0, 1]) # Left boundary: u = 0, right boundary: u = 1"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 217,
   "id": "1bf9bf87",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "Measure()"
      ]
     },
     "execution_count": 217,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "degree = 1\n",
    "Ω = Triangulation(model)\n",
    "dΩ = Measure(Ω,degree)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 218,
   "id": "19b7a3c7",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "b (generic function with 1 method)"
      ]
     },
     "execution_count": 218,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# working version \n",
    "ff(x) = 0\n",
    "velocity(x) = VectorValue( 1.0, 0.0 )\n",
    "a(u,v) = ∫( 0.1*∇(v)⋅∇(u) + (velocity ⋅ ∇(u))*v)*dΩ\n",
    "b(v) = ∫( v*ff )*dΩ"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 219,
   "id": "728948b9",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "AffineFEOperator()"
      ]
     },
     "execution_count": 219,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "op = AffineFEOperator(a,b,Uh,Vh)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 220,
   "id": "5bb1b54c",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "LinearFESolver()"
      ]
     },
     "execution_count": 220,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "ls = LUSolver()\n",
    "solver = LinearFESolver(ls)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 221,
   "id": "7471a56b",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "SingleFieldFEFunction():\n",
       " num_cells: 300\n",
       " DomainStyle: ReferenceDomain()\n",
       " Triangulation: BodyFittedTriangulation()\n",
       " Triangulation id: 221304387296365017"
      ]
     },
     "execution_count": 221,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "uh = solve(solver,op)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 222,
   "id": "3c4b7e6f",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "([\"results.vtu\"],)"
      ]
     },
     "execution_count": 222,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "writevtk(Ω,\"results\",cellfields=[\"uh\"=>uh])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "3d478320",
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "54785f07",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Julia 1.8.1",
   "language": "julia",
   "name": "julia-1.8"
  },
  "language_info": {
   "file_extension": ".jl",
   "mimetype": "application/julia",
   "name": "julia",
   "version": "1.8.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
