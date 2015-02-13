#!/bin/python
from dolfin import *
import mshr

geom = mshr.Rectangle(Point(0.,0.), Point(1., 3.)) #+ mshr.Rectangle(Point(0.,0.), Point(2., 1.))
mesh = mshr.generate_mesh(geom, 30)

bndry = FacetFunction('size_t', mesh)
for facet in facets(mesh):
    mp = facet.midpoint()
    bndry[facet] = 1 if (near(mp[0], 0.) or near(mp[0], 1.) or near(mp[1], 0.) or near(mp[1], 3.)) else 0
#plot(bndry, interactive=True)

U = FunctionSpace(mesh, 'Lagrange', 1)
u, v = TrialFunction(U), TestFunction(U)

F = u*v*dx + inner(grad(u), grad(v))*dx + Expression('1+0*x[0]*x[1]')*v*dx
w = Function(U)
solve(lhs(F) == rhs(F), w, DirichletBC(U, 0., bndry, 1))
plot(w, interactive=True, axes=True)


