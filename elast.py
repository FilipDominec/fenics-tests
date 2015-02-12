# -*- coding: utf-8 -*-
from dolfin import *
from mshr import *

# Build mesh
L = 3.
geometry =  Rectangle(Point(0.,0.), Point(L, 1.))-Circle(Point(1.,.45), .2, 30)
mesh = generate_mesh(geometry, 100)
#TODO mesh = Mesh('lego_beam.xml')

# Construct facet markers
bndry = FacetFunction("size_t", mesh)
for f in facets(mesh):
     mp = f.midpoint()
     if near(mp[0], 0.0): bndry[f] = 1  # front side
     elif near(mp[0], L): bndry[f] = 2  # rear side



# Build function spaces (Taylor-Hood??)
U = VectorFunctionSpace(mesh, "Lagrange", 2) ## -- displacement
V = VectorFunctionSpace(mesh, "Lagrange", 2) ## -- velocity
P = FunctionSpace(mesh, "Lagrange", 1)	     ## -- pressure
W = MixedFunctionSpace([U, V, P])

bc_force = DirichletBC(W.sub(1),Constant((1., 0.)), bndry, 1)
bc_force = DirichletBC(W.sub(1),Constant((1., 0.)), bndry, 1)

## TODO build the system of equations
## Equations 
# ∂u/∂t = v			∂v/∂t = div P			J = 1
# ∫∂u/∂t φ dx  =  
# ...


	
plot(bndry, interactive=True)



## TODO gt=Constant(0.), and later: gt.assign(100*t)
## * difference between the 0 and " registers is that 0 is only populated with yanked text
## * qstat; qsub -I -l nodes=1:ppn=1 -q -queue1 -X



