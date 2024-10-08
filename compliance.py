# HOW TO RUN PROGRAM CODE:
# python3 minimumComplianceTAO.py -tao_monitor -tax_max_it 100
# Option: -tao_monotor prints the Function value at each iteration
# Option: -tao_max_it specify maximum number of iterations. Here we set to 100
# where P_{\varepsilon}(\rho) is the Modica-Mortola perimeter term
# alpha = 10^{-2} is the penalization parameter
# f = Constant((0, -1)) is a downward boundary force
# \rho is the design variable i.e \rho = 0 means "void" or "hole" and \rho = 1 means "solid" or "material"
# \Omega is the design domain. It is a cantilever beam with length 1 and width 1/3 clamped on its left side \Gamma_D

# This problem is self-adjoint. Lucky for us, we do not need to solve the adjoint PDE to compute sensitivity

def parse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-tao_monitor', '--tao_monitor', action = 'store_true', help = 'TAO monitor')
    parser.add_argument('-l', '--lagrange', type = float, default = 5.0, help = 'Lagrange multiplier for structural material')
    parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
    parser.add_argument('-tao_gatol', '--tao_gatol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is less than this')
    parser.add_argument('-tao_grtol', '--tao_grtol', type = float, default = 1.0e-7, help = 'Stop if relative norm of gradient is less than this')
    parser.add_argument('-tao_gttol', '--tao_gttol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is reduced by this factor')
    parser.add_argument('-v', '--volume', type = float, default = 0.4, help = 'Volume percentage for structural material')
    parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of Modica-Mortola')
    parser.add_argument('-e', '--epsilon', type = float, default = 5.0e-3, help = 'Phase-field regularization parameter')
    parser.add_argument('-o', '--output', type = str, default = 'output', help = 'Output folder')
    parser.add_argument('-m', '--mesh', type = str, default = 'cantileverBeam.msh', help = 'Dimensions of meshed beam')
    options = parser.parse_args()
    return options

options = parse()

from firedrake import *
from petsc4py import PETSc

# Import "gmesh" mesh
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) # Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1)

# Create initial design
###### Begin Initial Design #####
rho = Function(V, name = "Material density")
rho_i = Function(V, name = "Material density")
x, y = SpatialCoordinate(mesh)
rho.interpolate(Constant(options.volume))
###### End Initial Design #####

# Total volume of the domain |omega|
rho_i.interpolate(Constant(1.0))
omega = assemble(rho_i * dx)


# Define the constant parameters used in the problem
alpha = options.kappa # Perimeter weight
lagrange = options.lagrange # Lagrange multiplier for Volume constraint
delta = 1.0e-3 
epsilon = options.epsilon

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
func4 = lagrange * rho * dx  # Volume constraint

J = func1 + func2 + func3 + func4

# Define the weak form for forward PDE
a_forward = h(rho) * inner(sigma(u, Id), epsilon(v)) * dx
L_forward = inner(f, v) * ds(8)
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_legrange = h(rho) * inner(sigma(u, Id), epsilon(u)) * dx
L_legrange = inner(f, u) * ds(8)
R_legrange = a_legrange - L_legrange
L = J - R_legrange

beam = VTKFile(options.output + '/beam.pvd')

def FormObjectiveGradient(tao, x, G):

	i = tao.getIterationNumber()
	if (i%5) == 0:
		# Save output files after each 5 iterations
		rho_i.interpolate(rho)
		beam.write(rho_i, u, time = i)

	with rho.dat.vec as rho_vec:
		rho_vec.set(0.0)
		rho_vec.axpy(1.0, x)

	volume = assemble(rho * dx)/omega
	print("The volume fraction is {}".format(volume))
	print(" ")

	# Solve forward PDE
	solve(R_fwd == 0, u, bcs = bcs)

	dJdrho = assemble(derivative(L, rho))
	with dJdrho.dat.vec as dJdrho_vec:
		G.set(0.0)
		G.axpy(1.0, dJdrho_vec)

	f_val = assemble(J)
	return f_val

# Setting lower and upper bounds
lb = Function(V)
ub = Function(V)
lb.interpolate(Constant(0.0))
ub.interpolate(Constant(1.0))

with lb.dat.vec as lb_vec:
	rho_lb = lb_vec

with ub.dat.vec as ub_vec:
	rho_ub = ub_vec

# Setting TAO solver
tao = PETSc.TAO().create(PETSc.COMM_SELF)
tao.setType('bncg')
tao.setObjectiveGradient(FormObjectiveGradient, None)
tao.setVariableBounds(rho_lb, rho_ub)
tao.setFromOptions()

# Initial design guess
with rho.dat.vec as rho_vec:
	x = rho_vec.copy()

# Solve the optimization problem
tao.solve(x)
tao.destroy()
