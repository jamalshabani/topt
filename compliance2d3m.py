def parse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--iterations', type = int, default = 100, help = 'Number of iterations')
    parser.add_argument('-tol', '--tolerance', type = float, default = 1.0e-5, help = 'Convergence tolerance')
    parser.add_argument('-e2', '--e2modulus', type = float, default = 0.01, help = 'Elastic Modulus for material 2')
    parser.add_argument('-e3', '--e3modulus', type = float, default = 1.0, help = 'Elastic Modulus for material 3')
    parser.add_argument('-v2', '--volume2', type = float, default = 0.3, help = 'Target volume for structural material')
    parser.add_argument('-v3', '--volume3', type = float, default = 0.3, help = 'Target volume for responsive material')
    parser.add_argument('-pncg', '--pncg', type = str, default = 'fr', help = 'Nonlinear conjugate gradient method type')
    parser.add_argument('-v', '--volume', type = float, default = 0.4, help = 'Volume fraction occupied by solid')
    parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of the perimeter')
    parser.add_argument('-e', '--epsilon', type = float, default = 2.5e-3, help = 'Phase-field regularization parameter')
    parser.add_argument('-o', '--output', type = str, default = 'output', help = 'Output folder')
    parser.add_argument('-m', '--mesh', type = str, default = 'cantileverBeam.msh', help = 'Dimensions of meshed beam')
    options = parser.parse_args()
    return options

options = parse()

# Get nonlinear conjugate gradient method type
pncg_type = options.pncg

from firedrake import *

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

rho2 = Function(V, name = "Material 1")  # Structural material 1(Blue)
rho3 = Function(V, name = "Material 2")  # Responsive material 2(Red)

dJdrho2 = Function(V, name = "Gradient w.r.t rho2")
dJdrho3 = Function(V, name = "Gradient w.r.t rho3")

projdJdrho2 = Function(V, name = "Projected gradient2")
projdJdrho3 = Function(V, name = "Projected gradient3")

prevdJdrho2 = Function(V, name = "Previous gradient2")
prevdJdrho3 = Function(V, name = "Previous gradient3")

deltadJdrho2 = Function(V, name = "Delta gradient2")
deltadJdrho3 = Function(V, name = "Delta gradient3")

rho2.interpolate(Constant(options.volume2))
rho3.interpolate(Constant(options.volume3))
###### End Initial Design #####

# Total volume of the domain |omega|
omega = assemble(Function(V).interpolate(1.0) * dx)

# Define the constant parameters used in the problem
kappa = options.kappa # Perimeter weight
delta = 1.0e-6 
epsilon = options.epsilon

kappa_d_e = kappa / epsilon
kappa_m_e = kappa * epsilon
iterations = options.iterations

# Downward Load
f = Constant((0, -1))

# Youngs modulii of the beam and poisson ratio
E1 = Constant(delta)
E2 = Constant(options.e2modulus)
E3 = Constant(options.e3modulus)
nu = Constant(0.3) #nu poisson ratio

mu1 = E1/(2 * (1 + nu))
lambda1 = (E1 * nu)/((1 + nu) * (1 - 2 * nu))

mu2 = E2/(2 * (1 + nu))
lambda2 = (E2 * nu)/((1 + nu) * (1 - 2 * nu))

mu3 = E3/(2 * (1 + nu))
lambda3 = (E3 * nu)/((1 + nu) * (1 - 2 * nu))

# Define h(x)=x^2
def h(x):
	return pow(x, 2)

# Define the triple-well potential function
# W(x, y) = (1 - x - y)^2 * (x + y)^2 + (1 - x)^2 * x^2 + (1 - y)^2 * y^2
def W(rho):
	material1 = pow((1 - rho2 - rho3), 2) * pow((rho2 + rho3), 2)
	material2 = pow((1 - rho2), 2) * pow(rho2, 2)
	material3 = pow((1 - rho3), 2) * pow(rho3, 2)
	return material1 + material2 + material3


# Define stress and strain tensors
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

# Define the stress tensor sigma1(u) for void
def sigma1(u, Id):
	return lambda1 * tr(epsilon(u)) * Id + 2 * mu1 * epsilon(u)

# Define the stress tensor sigma2(u) for structural material
def sigma2(u, Id):
	return lambda2 * tr(epsilon(u)) * Id + 2 * mu2 * epsilon(u)

# Define the stress tensor sigma3(u) for responsive material
def sigma3(u, Id):
	return lambda3 * tr(epsilon(u)) * Id + 2 * mu3 * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV, name = "Displacement")

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the objective function
func1 = inner(f, u) * ds(8)
func2 = kappa_d_e * W(rho) * dx

func3_sub1 = inner(grad(1 - rho2 - rho3), grad(1 - rho2 - rho3)) * dx
func3_sub2 = inner(grad(rho2), grad(rho2)) * dx
func3_sub3 = inner(grad(rho3), grad(rho3)) * dx

func3 = kappa_m_e * (func3_sub1 + func3_sub2 + func3_sub3)

J = func1 + func2 + func3

# Define the weak form for forward PDE
a_forward1 = h(1 - rho2 - rho3) * inner(sigma1(u, Id), epsilon(v)) * dx
a_forward2 = h(rho2) * inner(sigma2(u, Id), epsilon(v)) * dx
a_forward3 = h(rho3) * inner(sigma3(u, Id), epsilon(v)) * dx
a_forward = a_forward1 + a_forward2 + a_forward3

L_forward = inner(f, v) * ds(8)
R_fwd = a_forward - L_forward

# Define the Lagrangian
a_lagrange1 = h(1 - rho2 - rho3) * inner(sigma1(u, Id), epsilon(u)) * dx
a_lagrange2 = h(rho2) * inner(sigma2(u, Id), epsilon(u)) * dx
a_lagrange3 = h(rho3) * inner(sigma3(u, Id), epsilon(u)) * dx
a_lagrange   = a_lagrange1 + a_lagrange2 + a_lagrange3

L_lagrange = inner(f, u) * ds(8)

R_lagrange = a_lagrange - L_lagrange
L = J - R_lagrange

beam = VTKFile(options.output + '/beam.pvd')

# TO DO!
# 1. Add 3D support

def projectedNonlinearConjugateGradient(type):
    stepSize = 100
    c = 0.95
    beta = 0.5

    # Solve forward PDE
    solve(R_fwd == 0, u, bcs = bcs)
    currentObjValue = assemble(J)

    # Compute gradients w.r.t to designs
    dJdrho2.interpolate(assemble(derivative(L, rho2)).riesz_representation(riesz_map = "l2"))
    dJdrho3.interpolate(assemble(derivative(L, rho3)).riesz_representation(riesz_map = "l2"))

    # Do gradient projection into appropriate spaces for volume constraint
    projdJdrho2.interpolate(dJdrho2 - assemble(dJdrho2 * dx)/omega)
    projdJdrho3.interpolate(dJdrho3 - assemble(dJdrho3 * dx)/omega)

    if type == 'gd':
        alpha2 = 0
        alpha3 = 0

    elif type == 'fr':
        if assemble(inner(prevdJdrho2, prevdJdrho2) * dx) == 0.0:
            alpha2 = 0
            alpha3 = 0
        else:
            alpha2 = assemble(inner(projdJdrho2, projdJdrho2) * dx) / assemble(inner(prevdJdrho2, prevdJdrho2) * dx)
            alpha3 = assemble(inner(projdJdrho3, projdJdrho3) * dx) / assemble(inner(prevdJdrho3, prevdJdrho3) * dx)

    prevdJdrho2.interpolate(projdJdrho2)
    prevdJdrho3.interpolate(projdJdrho3)
    
    projdJdrho2.interpolate(projdJdrho2 + alpha2 * projdJdrho2)
    projdJdrho3.interpolate(projdJdrho3 + alpha3 * projdJdrho3)
    # Update design
    rho2.interpolate(rho2 - stepSize * projdJdrho2)
    rho3.interpolate(rho3 - stepSize * projdJdrho3)
    
    # Solve forward PDE
    solve(R_fwd == 0, u, bcs = bcs)
    nextObjValue = assemble(J)

    if nextObjValue > currentObjValue - c * stepSize * (assemble(inner(projdJdrho2, projdJdrho2) * dx) + assemble(inner(projdJdrho3, projdJdrho3) * dx)):
        stepSize = stepSize * beta

        # Update design
        rho2.interpolate(rho2 - stepSize * projdJdrho2)
        rho3.interpolate(rho3 - stepSize * projdJdrho2)

        # Solve forward PDE
        solve(R_fwd == 0, u, bcs = bcs)
        nextObjValue = assemble(J)
    else:
        # Update design
        rho2.interpolate(rho2 - stepSize * projdJdrho2)
        rho3.interpolate(rho3 - stepSize * projdJdrho3)

    objResidual = abs(nextObjValue - currentObjValue)

    # Compute volume fractions
    volume2 = assemble(rho2 * dx)/omega
    volume3 = assemble(rho3 * dx)/omega

    return rho2, rho3, u, volume2, volume3, currentObjValue, dJdrho2, dJdrho3, objResidual
    
converged = False

if __name__ == "__main__":

    i = 0
    while i < iterations and not converged:

        rho2, rho3, u, volume2, volume3, objValue, dJdrho2, dJdrho3, objResidual = projectedNonlinearConjugateGradient(pncg_type)

        # Convergence check
        if i == 0:
            dJdrhoNorm0 = assemble(inner(dJdrho2, dJdrho2) * dx)

        dJdrhoNorm = assemble(inner(dJdrho2, dJdrho2) * dx)
        residual = dJdrhoNorm/dJdrhoNorm0

        if residual < options.tolerance or objResidual < options.tolerance:
            converged = True

        rho_i.interpolate(rho3 - rho2)

        print("Iteration.: {0} Obj.: {1:.5f} Vol2.: {2:.2f} Vol3.: {3:.2f} |Delta F|.: {4:.6f}".format(i + 1, objValue, volume2, volume3, objResidual))

        # Save designs for processing using Paraview or VisIt
        if i%10 == 0:
            rho2.interpolate(conditional(gt(rho2, Constant(1.0)), Constant(1.0), rho2))
            rho2.interpolate(conditional(lt(rho2, Constant(0.0)), Constant(0.0), rho2))

            rho3.interpolate(conditional(gt(rho3, Constant(1.0)), Constant(1.0), rho3))
            rho3.interpolate(conditional(lt(rho3, Constant(0.0)), Constant(0.0), rho3))

            rho_i.interpolate(conditional(gt(rho_i, Constant(1.0)), Constant(1.0), rho_i))
            rho_i.interpolate(conditional(lt(rho_i, Constant(-1.0)), Constant(-1.0), rho_i))
            beam.write(rho2, rho3, rho_i, u, time = i)

        i = i + 1
            
    if converged:
        print("Converged after {} iterations......".format(i))
    else:
        print("Maximum iterations reached......")
