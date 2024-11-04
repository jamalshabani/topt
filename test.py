def parse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--iterations', type = int, default = 5000, help = 'Number of iterations')
    parser.add_argument('-tol', '--tolerance', type = float, default = 1.0e-6, help = 'Convergence tolerance')
    parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of the perimeter')
    parser.add_argument('-e', '--epsilon', type = float, default = 2.5e-3, help = 'Phase-field regularization parameter')
    parser.add_argument('-o', '--output', type = str, default = 'output', help = 'Output folder')
    parser.add_argument('-m', '--mesh', type = str, default = 'cantileverBeam.msh', help = 'Mesh file in .msh format')
    options = parser.parse_args()
    return options
options = parse()

from firedrake import *

# Import "gmesh" mesh
mesh = Mesh(options.mesh)

Id = Identity(mesh.geometric_dimension()) # Identity tensor

E = [1.0, 1.0, 1.0e-6] #Young's modulus for each materials with index -1 as void
volFractions = [0.4, 0.2, 0.6] #Volume fractions with index -1 as void
nMaterials = len(volFractions)  #Number of materials
nu = 0.3 #nu poisson ratio
pDimension = 2

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
D = VectorFunctionSpace(mesh, 'CG', 1, dim = nMaterials)
U = VectorFunctionSpace(mesh, 'CG', 1, dim = pDimension)

# Create initial design guess
###### Begin Initial Design#####
rhoList = []
muList = []
lambaList = []
for i in range(nMaterials):
    dFunction = Function(V, name = "Material {0}".format(i)).interpolate(Constant(volFractions[i]))
    rhoList.append(dFunction)

    muList.append(E[i]/(2 * (1 + nu)))
    lambaList.append((E[i] * nu)/((1 + nu) * (1 - 2 * nu)))

rho = Function(D).interpolate(as_vector([rhoList[i] for i in range(nMaterials)]))

a, b = rho.subfunctions[0], rho.subfunctions[1]
print(a, b)
###### End Initial Design#####

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

print("Hahahaah")