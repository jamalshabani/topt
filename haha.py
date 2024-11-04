def parse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--iterations', type = int, default = 2000 , help = 'Number of iterations')
    parser.add_argument('-tol', '--tolerance', type = float, default = 1.0e-6, help = 'Convergence tolerance')
    parser.add_argument('-k', '--kappa', type = float, default = 1.0e-1, help = 'Weight of the perimeter')
    parser.add_argument('-e', '--epsilon', type = float, default = 2.5e-3, help = 'Phase-field regularization parameter')
    parser.add_argument('-o', '--output', type = str, default = 'beam', help = 'Output folder')
    parser.add_argument('-m', '--mesh', type = str, default = 'cantileverBeam.msh', help = 'Mesh file in .msh format')
    options = parser.parse_args()
    return options
options = parse()

from firedrake import *

# Import "gmesh" mesh
mesh = Mesh(options.mesh)

Id = Identity(mesh.geometric_dimension()) # Identity tensor

E = [1.0, 1.0e-6] #Young's modulus for each materials with index -1 as void
volFractions = [0.4, 0.6] #Volume fractions with index -1 as void
nMaterials = len(volFractions)  #Number of materials
nu = 0.3 #nu poisson ratio
pDimension = 2

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
D = MixedFunctionSpace([V for i in range(nMaterials)])
U = VectorFunctionSpace(mesh, 'CG', 1, dim = pDimension)

# Create initial design guess
###### Begin Initial Design#####
rhoList = []
muList = []
lambdaList = []
rho = Function(D, name = "Material densities")
for i in range(nMaterials):
    (rho.sub(i)).rename("Material {0}".format(i))
    (rho.sub(i)).interpolate(Constant(volFractions[i]))
    muList.append(E[i]/(2 * (1 + nu)))
    lambdaList.append((E[i] * nu)/((1 + nu) * (1 - 2 * nu)))
###### End Initial Design#####

# Total volume of the domain |omega|
omega = assemble(Function(V).interpolate(1.0) * dx)

# Define the constant parameters used in the problem
kappa = options.kappa # Perimeter weight
epsilon = options.epsilon

kappa_d_e = kappa / epsilon
kappa_m_e = kappa * epsilon
iterations = options.iterations

# Downward Load
f = Constant((0, -1))

def rhoSub(index, rho):
    if index < nMaterials - 1:
        return rho.sub(index)
    else:
        (rho.sub(index)).interpolate(1 - sum(rho.subfunctions[ : nMaterials - 1]))
        return rho.sub(index)

# Define the multi-well potential
def W(rho):
    wList = []
    for i in range(nMaterials):
        wList.append(pow(rhoSub(i, rho), 2) * pow((1 - rhoSub(i, rho)), 2))
    return sum(wList)

# Define strain tensor epsilon(u)
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

# Define the stress tensor sigma(u)
def sigmaSimp(rho, u, Id):
    simpList = []
    sigmaList = []
    for i in range(nMaterials):
        simpList.append(pow(rhoSub(i, rho), 2))
        sigmaList.append(lambdaList[i] * tr(epsilon(u)) * Id + 2 * muList[i] * epsilon(u))
    return simpList, sigmaList


# Define test function and beam displacement
v = TestFunction(U)
u = Function(U, name = "Displacement")
p = Function(U, name = "Adjoint variable")

# The left side of the beam is clamped
bcs = DirichletBC(U, Constant((0, 0)), 7)

# Define the objective function
func1 = inner(f, u) * ds(8)
func2 = kappa_d_e * W(rho) * dx
func3 = kappa_m_e * inner(grad(rho), grad(rho)) * dx

J = func1 + func2 + func3

# Define the weak form for forward PDE
def bilinearLinearForms(rho, u, v, Id):
    bilinearList = []
    simpList, sigmaList = sigmaSimp(rho, u, Id)
    for i in range(nMaterials):
        bilinearList.append(simpList[i] * inner(sigmaList[i], epsilon(v)) * dx)
    linearForm = inner(f, v) * ds(8)
    return sum(bilinearList) - linearForm

# Define the Lagrangian
L = J - bilinearLinearForms(rho, u, u, Id)

beam = VTKFile(options.output + '/beam.pvd')

dJdrho = Function(D, name = "Gradient w.r.t rho")
projdJdrho = Function(D, name = "Projected gradient")
def projectedNonlinearConjugateGradient():
    stepSize = 100

    iteration = 0
    while iteration < options.iterations:
        # Solve forward PDE
        solve(bilinearLinearForms(rho, u, v, Id) == 0, u, bcs = bcs)
        objValue = assemble(J)
    # Update designs
        for i in range(nMaterials - 1):
            (dJdrho.sub(i)).interpolate(assemble(derivative(L, rho.sub(i))).riesz_representation(riesz_map = "l2"))
            (projdJdrho.sub(i)).interpolate(dJdrho.sub(i) - assemble(dJdrho.sub(i) * dx)/omega)
            (rho.sub(i)).interpolate(rho.sub(i) - stepSize * projdJdrho.sub(i))

        (rho.sub(nMaterials - 1)).interpolate(1 - sum(rho.subfunctions[ : nMaterials - 1]))
        iteration = iteration + 1

        volumeFractions = []
        for i in range(nMaterials):
            volumeFractions.append(assemble(rho.sub(i) * dx)/omega)
        print("Iteration.: {0} Obj.: {1:.5f} Vol.: {2:.5f} ".format(iteration, objValue, volumeFractions[0]))

        if iteration % 10 == 0:
            for i in range(nMaterials):
                (rho.sub(i)).interpolate(conditional(gt(rho.sub(i), Constant(1.0)), Constant(1.0), rho.sub(i)))
                (rho.sub(i)).interpolate(conditional(lt(rho.sub(i), Constant(0.0)), Constant(0.0), rho.sub(i)))
            beam.write(*rho.subfunctions, u, time = iteration)
    return rho, u

if __name__ == "__main__":

    rho, u = projectedNonlinearConjugateGradient()
