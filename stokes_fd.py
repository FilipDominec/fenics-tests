
from dolfin import *
import mshr

# Define domain
center = Point(0.3, 0.2)
radius = 0.05
L = .8
W = 0.41
geometry =  mshr.Rectangle(Point(0.0,0.0), Point(L, W)) \
           -mshr.Circle(Point(.3,.2), radius, 30) -mshr.Circle(Point(.25,.3), radius, 30)

# Build mesh
mesh = mshr.generate_mesh(geometry, 100)

# Construct facet markers
bndry = FacetFunction("size_t", mesh)
for f in facets(mesh):
     mp = f.midpoint()
     if near(mp[0], 0.0): bndry[f] = 1  # inflow
     elif near(mp[0], L): bndry[f] = 2  # outflow
     elif near(mp[1], 0.0) or near(mp[1], W): bndry[f] = 3  # walls
     elif mp.distance(center) <= radius:      bndry[f] = 3  # cylinder
plot(bndry, interactive=True)

# define dsN = Measure("ds", subdomain_id=1, subdomain_data=boundary_parts) TODO

# Build function spaces (Taylor-Hood)
## ..-Babushka rules of stability require L2-space for motion, and L1-space for pressure
V = VectorFunctionSpace(mesh, "Lagrange", 2) 
P = FunctionSpace(mesh, "Lagrange", 1)
W = MixedFunctionSpace([V, P])


noslip = Constant((0, 0))
## here we can access the velocity as `W.sub(0)' (not directly as V)
bc_in    = DirichletBC(W.sub(0), Constant((1.0,0)), bndry, 1)
#bc_out   = DirichletBC(W.sub(0), Constant((1.0,0)), bndry, 2)
bc_walls = DirichletBC(W.sub(0), Constant((0.0,0)), bndry, 3)

# Define unknown and test function(s)
(v_, p_) = TestFunctions(W)  ## version with -s at the end returns a tuple
(v , p)  = TrialFunctions(W)

def a(u,v): return inner(grad(u), grad(v))*dx
def b(p,v): return p*div(v)*dx
def L(v):   return inner(f, v)*dx

## == Stokes pro
#f = Expression("0")
#F = a(v,v_) + b(p,v_) + b(p_,v) - (inner(f, p_)*dx)  
f = Expression(("0","0"))
F = a(v,v_) + b(p,v_) + b(p_,v) - (inner(f, v_)*dx) #TODO

# w = Function(W)
# solve(lhs(F)==rhs(F), w, bcs)
w = Function(W)
solve(lhs(F)==rhs(F), w, [bc_in, bc_walls])
(v,p) = w.split()

plot(p, interactive=True, axes=True)
plot(v, interactive=True, axes=True)

