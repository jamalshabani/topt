from firedrake import *

# Import "gmesh" mesh
mesh = Mesh("cantileverBeam.msh")

Id = Identity(mesh.geometric_dimension()) # Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1)

# Create initial design
###### Begin Initial Design #####
rho = Function(V, name = "Material density")
rho_i = Function(V, name = "Material density")
x, y = SpatialCoordinate(mesh)
rho.interpolate(Constant(0.4))
###### End Initial Design #####

# Total volume of the domain |omega|
omega = assemble(Function(V).interpolate(1.0) * dx)

# Define the constant parameters used in the problem
alpha = 1.0e-2 # Perimeter weight
delta = 1.0e-3 
epsilon = 5.0e-5

alpha_d_e = alpha / epsilon
alpha_m_e = alpha * epsilon

# Downward traction force on the right corner
f = Constant((0, -1))

# Young's modulus of the beam and poisson ratio
E = 1.0
nu = 0.3 #nu poisson ratio

mu = E/(2 * (1 + nu))
_lambda = (E * nu)/((1 + nu) * (1 - 2 * nu))

# Define h(x)=x^2
def h(rho):
    return delta * (1 - rho)**2 + rho**2

# Define double-well potential W(a) function i.e W(x) = x^2(1 - x)^2
def W(rho):
    return pow(rho, 2) * pow((1 - rho), 2)

# Define stress and strain tensors
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

def sigma(u, Id):
    return _lambda * tr(epsilon(u)) * Id + 2 * mu * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV, name = "Displacement")

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
func1 = inner(f, u) * ds(8)
func2 = alpha_d_e * W(rho) * dx
func3 = alpha_m_e * inner(grad(rho), grad(rho)) * dx

J = func1 + func2 + func3

# Define the weak form for forward PDE
a_forward = h(rho) * inner(sigma(u, Id), epsilon(v)) * dx
L_forward = inner(f, v) * ds(8)
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_legrange = h(rho) * inner(sigma(u, Id), epsilon(u)) * dx
L_legrange = inner(f, u) * ds(8)
R_legrange = a_legrange - L_legrange
L = J - R_legrange

beam = VTKFile('output/beam.pvd')

# implement line search

def projectGradientDescent():

    dJdrho = Function(V, name = "Gradient w.r.t rho")
    volume = assemble(rho * dx) * 3

    # Solve forward PDE
    solve(R_fwd == 0, u, bcs = bcs)

    dJdrho.interpolate(assemble(derivative(L, rho)).riesz_representation(riesz_map="l2"))
    #dJdrho.interpolate(dJdrho - assemble(dJdrho * dx)/omega)
    rho.interpolate(rho - 5.0 * dJdrho)
    return rho




if __name__ == "__main__":
    for i in range(5000):
        print(i)
        rho = projectGradientDescent()

        if i%10 == 0:
            beam.write(rho)

    