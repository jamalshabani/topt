def parse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--iterations', type = int, default = 100, help = 'Number of iterations')
    parser.add_argument('-tol', '--tolerance', type = float, default = 1.0e-6, help = 'Convergence tolerance')
    parser.add_argument('-pncg', '--pncg', type = str, default = 'gd', help = 'Nonlinear conjugate gradient method type')
    parser.add_argument('-v', '--volume', type = float, default = 0.4, help = 'Volume fraction occupied by solid')
    parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of the perimeter')
    parser.add_argument('-e', '--epsilon', type = float, default = 5.0e-5, help = 'Phase-field regularization parameter')
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
dJdrho = Function(V, name = "Gradient w.r.t rho")
projdJdrho = Function(V, name = "Projected gradient")
prevdJdrho = Function(V, name = "Previous gradient")
currdJdrho = Function(V, name = "Current gradient")
deltadJdrho = Function(V, name = "Delta gradient")
x, y = SpatialCoordinate(mesh)
rho.interpolate(Constant(options.volume))
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
func2 = kappa_d_e * W(rho) * dx
func3 = kappa_m_e * inner(grad(rho), grad(rho)) * dx

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

    # Compute gradients w.r.t to design 
    dJdrho.interpolate(assemble(derivative(L, rho)).riesz_representation(riesz_map = "l2"))

    # Do gradient projection into appropriate spaces for volume constraint
    projdJdrho.interpolate(dJdrho - assemble(dJdrho * dx)/omega)

    if type == 'gd':
        alpha = 0

    elif type == 'fr':
        if assemble(inner(prevdJdrho, prevdJdrho) * dx) == 0.0:
            alpha = 0
        else:
            alpha = min(1, assemble(inner(projdJdrho, projdJdrho) * dx) / assemble(inner(prevdJdrho, prevdJdrho) * dx))

    elif type == 'pr':
        if assemble(inner(prevdJdrho, prevdJdrho) * dx) == 0.0:
            alpha = 0
        else:
            deltadJdrho.interpolate(projdJdrho - prevdJdrho)
            alpha = min(1, assemble(inner(deltadJdrho, projdJdrho) * dx) / assemble(inner(prevdJdrho, prevdJdrho) * dx))

    elif type == 'hs':
        if assemble(inner(prevdJdrho, prevdJdrho) * dx) == 0.0:
            alpha = 0
        else:
            deltadJdrho.interpolate(projdJdrho - prevdJdrho)
            alpha = min(1, assemble(inner(deltadJdrho, projdJdrho) * dx) / assemble(inner(prevdJdrho, deltadJdrho) * dx))
    
    elif type == 'dy':
        if assemble(inner(prevdJdrho, prevdJdrho) * dx) == 0.0:
            alpha = 0
        else:
            deltadJdrho.interpolate(projdJdrho - prevdJdrho)
            alpha = min(1, assemble(inner(projdJdrho, projdJdrho) * dx) / assemble(inner(prevdJdrho, deltadJdrho) * dx))
    
    prevdJdrho.interpolate(projdJdrho)
    print(alpha)
    
    
    projdJdrho.interpolate(projdJdrho + alpha * projdJdrho)
    # Update design
    rho.interpolate(rho - stepSize * projdJdrho)
    
    # Solve forward PDE
    solve(R_fwd == 0, u, bcs = bcs)
    nextObjValue = assemble(J)

    if nextObjValue > currentObjValue - c * stepSize * assemble(inner(projdJdrho, projdJdrho) * dx):
        stepSize = stepSize * beta

        # Update design
        rho.interpolate(rho - stepSize * projdJdrho)

        # Solve forward PDE
        solve(R_fwd == 0, u, bcs = bcs)
        nextObjValue = assemble(J)
    else:
        # Update design
        rho.interpolate(rho - stepSize * projdJdrho)

    objResidual = abs(nextObjValue - currentObjValue)

    # Compute volume fractions
    volume = assemble(rho * dx)/omega

    return rho, u, volume, currentObjValue, dJdrho, objResidual
    
converged = False

if __name__ == "__main__":

    i = 0
    while i < iterations and not converged:

        rho, u, volume, objValue, dJdrho, objResidual = projectedNonlinearConjugateGradient(pncg_type)

        # Convergence check
        if i == 0:
            dJdrhoNorm0 = assemble(inner(dJdrho, dJdrho) * dx)

        dJdrhoNorm = assemble(inner(dJdrho, dJdrho) * dx)
        residual = dJdrhoNorm/dJdrhoNorm0

        if residual < options.tolerance or objResidual < options.tolerance:
            converged = True

        rho_i.interpolate(rho)

        print("Iteration.: {0} Obj.: {1:.5f} Vol.: {2:.2f} Residual.: {3:.6f} |Delta F|.: {4:.6f}".format(i + 1, objValue, volume, residual, objResidual))

        # Save designs for processing using Paraview or VisIt
        if i%10 == 0:
            rho_i.interpolate(conditional(gt(rho_i, Constant(1.0)), Constant(1.0), rho_i))
            rho_i.interpolate(conditional(lt(rho_i, Constant(0.0)), Constant(0.0), rho_i))
            beam.write(rho_i, u, time = i)

        i = i + 1
            
    if converged:
        print("Converged after {} iterations......".format(i))
    else:
        print("Maximum iterations reached......")
    