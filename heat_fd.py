#!/bin/python
# -*- coding: utf-8 -*-
"""
Weak formulation of the heat equation problem

a(u,v) is known as a bilinear form and L(v) as a linear form
dentify the terms with the unknown u and collect them in a(u,v), and similarly collect all terms with only known functions in L(v)
"""
import dolfin as df
mesh = df.UnitSquareMesh(16, 16, 'crossed')
V = df.FunctionSpace(mesh, 'Lagrange', 3)

#  def boundary(x, on_boundary):
#      return on_boundary
#  bc = DirichletBC(V, 0.0, boundary)

## Create boundary markers
boundary_parts = df.FacetFunction('size_t', mesh)
left   = df.AutoSubDomain(lambda x: df.near(x[0], 0.0))
right  = df.AutoSubDomain(lambda x: df.near(x[0], 1.0))
bottom = df.AutoSubDomain(lambda x: df.near(x[1], 0.0)) # or df.near(x[1], 1.0)
left  .mark(boundary_parts, 1)
right .mark(boundary_parts, 2)
bottom.mark(boundary_parts, 2)

u = df.Function(V)
v = df.TestFunction(V)

## for standard Poisson problem:    -Δu = f     ---->   -∫ Δu v dV  =  ∫∇u ∇v dV  -  (∫v ∂u/∂n ds)  =  ∫ f v dV
#  f = df.Expression("1")
#  a = df.inner(df.grad(u), df.grad(v)) * df.dx 
#  L = f * v * df.dx 

## for non-linear problem:           u_t - ∇(u² ∇u) = f  ---->  ∫u_t v dV  +  ∫u² ∇u ∇v dV  =  ∫ f v dV
#u_t = 0.1 K = 0.1
f = df.Expression("1")
## TODO assign , project
a = (df.inner(df.grad(u), df.grad(v))  - f * v)   *  df.dx

u = df.Function(V)
df.solve(a == 0, u, df.DirichletBC(V, 0.0, bottom))
df.plot(u, interactive=True, axes=True)





## for problem with convection:
# b = df.Expression(("-(x[1] - 0.5)", "x[0] - 0.5"))
# dsN = df.Measure("ds", subdomain_id=1, subdomain_data=boundary_parts)
# +  df.inner(b, df.grad(v))*v



