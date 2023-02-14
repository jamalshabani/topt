def parse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-tao_bncg_type', '--tao_bncg_type', type = str, default = 'ssml_bfgs', help = 'BNCG algorithm type')
    parser.add_argument('-tao_bncg_alpha', '--tao_bncg_alpha', type = float, default = 0.5, help = 'Scalar preconditioning')
    parser.add_argument('-tao_monitor', '--tao_monitor', action = 'store_true', help = 'TAO monitor')
    parser.add_argument('-ls', '--lagrange_s', type = float, default = 5.0, help = 'Lagrange multiplier for structural material')
    parser.add_argument('-lr', '--lagrange_r', type = float, default = 0.5, help = 'Lagrange multiplier for responsive material')
    parser.add_argument('-tao_converged_reason', '--tao_converged_reason', action = 'store_true', help = 'TAO convergence reason')
    parser.add_argument('-tao_ls_type', '--tao_ls_type', type = str, default = 'more-thuente', help = "TAO line search")
    parser.add_argument('-tao_view', '--tao_view', action = 'store_true', help = "View convergence details")
    parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
    parser.add_argument('-tao_gatol', '--tao_gatol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is less than this')
    parser.add_argument('-tao_grtol', '--tao_grtol', type = float, default = 1.0e-7, help = 'Stop if relative norm of gradient is less than this')
    parser.add_argument('-tao_gttol', '--tao_gttol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is reduced by this factor')
    parser.add_argument('-vs', '--volume_s', type = float, default = 0.4, help = 'Volume percentage for structural material')
    parser.add_argument('-vr', '--volume_r', type = float, default = 0.4, help = 'Volume percentage for responsive material')
    parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of Modica-Mortola')
    parser.add_argument('-e', '--epsilon', type = float, default = 5.0e-3, help = 'Phase-field regularization parameter')
    parser.add_argument('-o', '--output', type = str, default = 'output1', help = 'Output folder')
    parser.add_argument('-m', '--mesh', type = str, default = '1_to_6_mesh', help = 'Dimensions of meshed beam')
    parser.add_argument('-es', '--esmodulus', type = float, default = 0.1, help = 'Elastic Modulus for structural material')
    parser.add_argument('-er', '--ermodulus', type = float, default = 1.0, help = 'Elastic Modulus for responsive material')
    parser.add_argument('-p', '--power_p', type = float, default = 2.0, help = 'Power for elasticity interpolation')
    parser.add_argument('-q', '--power_q', type = float, default = 2.0, help = 'Power for multiple-well function')
    options = parser.parse_args()
    return options

options = parse()

from firedrake import *
from petsc4py import PETSc
import time
import numpy as np

start = time.time()

# Import gmesh
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
# Try using DG0 function space

V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

# Create initial design
###### Begin Initial Design #####
mesh_coordinates = mesh.coordinates.dat.data[:]
M = len(mesh_coordinates)

rho =  Function(VV, name = "Design variable")
rho_i = Function(V, name = "Material density")
rho2 = Function(V, name = "Structural material")  # Structural material 1(Blue)
rho3 = Function(V, name = "Responsive material")  # Responsive material 2(Red)

x, y = SpatialCoordinate(mesh)
rho2 = interpolate(Constant(options.volume_s), V)
rho3 = interpolate(Constant(options.volume_r), V)

rho2 = 0.5 + 0.5 * sin(4*pi*x) * sin(8*pi*y)
rho3 = 0.3 + 0.3 * cos(4*pi*x) * cos(8*pi*y)

rho2 = interpolate(rho2, V)
rho3 = interpolate(rho3, V)


rho = as_vector([rho2, rho3])
rho = interpolate(rho, VV)
###### End Initial Design #####

# Define the constant parameter used in the problem
kappa = Constant(options.kappa)
lagrange_r = Constant(options.lagrange_r)
lagrange_s = Constant(options.lagrange_s)

# Total volume of the domain |omega|
omega = assemble(interpolate(Constant(1.0), V) * dx)

delta = Constant(1.0e-3)
epsilon = Constant(options.epsilon)
kappa_d_e = Constant(kappa / epsilon)
kappa_m_e = Constant(kappa * epsilon)

# Define the boundary/traction force
f = Constant((0, -1))

# Young's modulus of the beam and poisson ratio
E_v = Constant(delta)
E_s = Constant(options.esmodulus)
E_r = Constant(options.ermodulus)
nu = Constant(0.3) #nu poisson ratio

mu_v = E_v/(2 * (1 + nu))
lambda_v = (E_v * nu)/((1 + nu) * (1 - 2 * nu))

mu_s = E_s/(2 * (1 + nu))
lambda_s = (E_s * nu)/((1 + nu) * (1 - 2 * nu))

mu_r = E_r/(2 * (1 + nu))
lambda_r = (E_r * nu)/((1 + nu) * (1 - 2 * nu))

def v_v(rho):
    return 1 - rho.sub(0) - rho.sub(1)

def v_s(rho):
    return rho.sub(0)

def v_r(rho):
    return rho.sub(1)

# Define h_v(rho)=rho_v^(p)
def h_v(rho):
    return pow((1 - rho.sub(0) - rho.sub(1)), options.power_p)

# Define h_s(rho)=rho_s^(p)
def h_s(rho):
    return pow(rho.sub(0), options.power_p)

# Define h_r(rho)=rho_r^(p)
def h_r(rho):
    return pow(rho.sub(1), options.power_p)


# Define the double-well potential function
# W(x, y) = (x + y)^q * (1 - x)^q * (1 - y)^q
def W(rho):
    return pow((rho.sub(0) + rho.sub(1)), options.power_q) * pow((1 - rho.sub(0)), options.power_q) * pow((1 - rho.sub(1)), options.power_q)

# Define strain tensor epsilon(u)
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

# Define the stress tensor sigma_v(u) for void
def sigma_v(u, Id):
    return lambda_v * tr(epsilon(u)) * Id + 2 * mu_v * epsilon(u)

# Define the stress tensor sigma_s(u) for structural material
def sigma_s(u, Id):
    return lambda_s * tr(epsilon(u)) * Id + 2 * mu_s * epsilon(u)

# Define the stress tensor sigma_r(u) for responsive material
def sigma_r(u, Id):
    return lambda_r * tr(epsilon(u)) * Id + 2 * mu_r * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV, name = "Displacement")
p = Function(VV, name = "Adjoint variable")

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
J = inner(f, u) * ds(8)
func1 = kappa_d_e * W(rho) * dx

func2_sub1 = inner(grad(v_v(rho)), grad(v_v(rho))) * dx
func2_sub2 = inner(grad(v_s(rho)), grad(v_s(rho))) * dx
func2_sub3 = inner(grad(v_r(rho)), grad(v_r(rho))) * dx

func2 = kappa_m_e * (func2_sub1 + func2_sub2 + func2_sub3)

func3 = lagrange_s * v_s(rho) * dx
func4 = lagrange_r * v_r(rho) * dx

# Objective function + Modica-Mortola functional + Volume constraint
P = func1 + func2 + func3 + func4
JJ = J + P

# Define the weak form for forward PDE
a_forward_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(v)) * dx
a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
a_forward = a_forward_v + a_forward_s + a_forward_r

L_forward = inner(f, v) * ds(8)
R_fwd = a_forward - L_forward

# Define the Lagrangian
# The problem is self-adjoint so we replace langrange multiplier(p) with u
a_lagrange_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(u)) * dx
a_lagrange_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(u)) * dx
a_lagrange_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(u)) * dx
a_lagrange   = a_lagrange_v + a_lagrange_s + a_lagrange_r

L_lagrange = inner(f, u) * ds(8)
R_lagrange = a_lagrange - L_lagrange
L = JJ - R_lagrange

# Beam .pvd file for saving designs
beam = File(options.output + '/beam.pvd')
dJdrho2 = Function(V)
dJdrho3 = Function(V)
dJdrho2_project = Function(V)
dJdrho3_project = Function(V)

N = M * 2
index_2 = []
index_3 = []

for i in range(N):
    if (i%2) == 0:
        index_2.append(i)
    if (i%2) == 1:
        index_3.append(i)

def FormObjectiveGradient(tao, x, G):

    i = tao.getIterationNumber()
    if (i%5) == 0:
        rho_i.interpolate(rho.sub(1) - rho.sub(0))
        rho2.interpolate(rho.sub(0))
        rho3.interpolate(rho.sub(1))
        beam.write(rho_i, rho2, rho3, u, time = i)

    with rho.dat.vec as rho_vec:
        rho_vec.set(0.0)
        rho_vec.axpy(1.0, x)

    # Solve forward PDE
    solve(R_fwd == 0, u, bcs = bcs)

    # Evaluate the objective function
    # objective_value = assemble(J)
    # print("The value of objective function is {}".format(objective_value))

    # Print volume fraction of structural material
    volume_s = assemble(v_s(rho) * dx)/omega
    print("The volume fraction(Vs) is {}".format(volume_s))

    # Print volume fraction of responsive material
    volume_r = assemble(v_r(rho) * dx)/omega
    print("The volume fraction(Vr) is {}".format(volume_r))
    print(" ")

    # Compute gradiet w.r.t rho2 and rho3
    dJdrho2.interpolate(assemble(derivative(L, rho.sub(0))))
    dJdrho3.interpolate(assemble(derivative(L, rho.sub(1))))
 
    dJdrho2_array = dJdrho2.vector().array()
    dJdrho3_array = dJdrho3.vector().array()

    G.setValues(index_2, dJdrho2_array)
    G.setValues(index_3, dJdrho3_array)

    # print(G.view())

    f_val = assemble(J)
    return f_val

# Setting lower and upper bounds
lb = as_vector((0, 0))
ub = as_vector((1, 1))
lb = interpolate(lb, VV)
ub = interpolate(ub, VV)

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

end = time.time()
print("\nExecution time (in seconds):", (end - start))
