from fenics import *
import numpy as np
import matplotlib.pyplot as plt


def solver_bcs(kappa, f, boundary_conditions, Nx, Ny,
               degree=1,
               subdomains=[],
               linear_solver='Krylov',
               abs_tol=1E-5,
               rel_tol=1E-3,
               max_iter=1000):
    # Create mesh and define function space
    mesh = UnitSquareMesh(Nx, Ny)
    V = FunctionSpace(mesh, 'P', degree)

    # Check if we have subdomains
    if subdomains:
        if not isinstance(kappa, (list, tuple, np.ndarray)):
            raise TypeError(
                'kappa must be array if we have sudomains, not %s'
                % type(kappa))
        #materials = CellFunction('size_t', mesh)
        materials = MeshFunction('size_t', mesh, 0)
        materials.set_all(0)
        for m, subdomain in enumerate(subdomains[1:], 1):
            subdomain.mark(materials, m)

        kappa_values = kappa
        V0 = FunctionSpace(mesh, 'DG', 0)
        kappa  = Function(V0)
        help = np.asarray(materials.array(), dtype=np.int32)
        kappa.vector()[:] = np.choose(help, kappa_values)
    else:
        if not isinstance(kappa, (Expression, Constant)):
            raise TypeError(
                'kappa is type %s, must be Expression or Constant'
                % type(kappa))

    # Define boundary subdomains
    tol = 1e-14

    class BoundaryX0(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0, tol)

    class BoundaryX1(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 1, tol)

    class BoundaryY0(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], 0, tol)

    class BoundaryY1(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], 1, tol)

    # Mark boundaries
    #boundary_markers = FacetFunction('size_t', mesh)
    boundary_markers = MeshFunction('size_t', mesh, 1)
    boundary_markers.set_all(9999)
    bx0 = BoundaryX0()
    bx1 = BoundaryX1()
    by0 = BoundaryY0()
    by1 = BoundaryY1()
    bx0.mark(boundary_markers, 0)
    bx1.mark(boundary_markers, 1)
    by0.mark(boundary_markers, 2)
    by1.mark(boundary_markers, 3)

    # Redefine boundary integration measure
    ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

    # Collect Dirichlet conditions
    bcs = []
    for i in boundary_conditions:
        if 'Dirichlet' in boundary_conditions[i]:
            bc = DirichletBC(V, boundary_conditions[i]['Dirichlet'],
                             boundary_markers, i)
            bcs.append(bc)

    debug1 = False
    if debug1:

        # Print all vertices that belong to the boundary parts
        for x in mesh.coordinates():
            if bx0.inside(x, True): print('%s is on x = 0' % x)
            if bx1.inside(x, True): print('%s is on x = 1' % x)
            if by0.inside(x, True): print('%s is on y = 0' % x)
            if by1.inside(x, True): print('%s is on y = 1' % x)

        # Print the Dirichlet conditions
        print('Number of Dirichlet conditions:', len(bcs))
        if V.ufl_element().degree() == 1:  # P1 elements
            d2v = dof_to_vertex_map(V)
            coor = mesh.coordinates()
            for i, bc in enumerate(bcs):
                print('Dirichlet condition %d' % i)
                boundary_values = bc.get_boundary_values()
                for dof in boundary_values:
                    print('   dof %2d: u = %g' % (dof, boundary_values[dof]))
                    if V.ufl_element().degree() == 1:
                        print('    at point %s' %
                              (str(tuple(coor[d2v[dof]].tolist()))))

    # Define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)

    # Collect Neumann integrals
    integrals_N = []
    for i in boundary_conditions:
        if 'Neumann' in boundary_conditions[i]:
            if boundary_conditions[i]['Neumann'] != 0:
                g = boundary_conditions[i]['Neumann']
                integrals_N.append(g*v*ds(i))

    # Collect Robin integrals
    integrals_R_a = []
    integrals_R_L = []
    for i in boundary_conditions:
        if 'Robin' in boundary_conditions[i]:
            r, s = boundary_conditions[i]['Robin']
            integrals_R_a.append(r*u*v*ds(i))
            integrals_R_L.append(r*s*v*ds(i))

    # Simpler Robin integrals
    integrals_R = []
    for i in boundary_conditions:
        if 'Robin' in boundary_conditions[i]:
            r, s = boundary_conditions[i]['Robin']
            integrals_R.append(r*(u - s)*v*ds(i))

    # Sum integrals to define variational problem
    a = kappa*dot(grad(u), grad(v))*dx + sum(integrals_R_a)
    L = f*v*dx - sum(integrals_N) + sum(integrals_R_L)

    # Simpler variational problem
    F = kappa*dot(grad(u), grad(v))*dx + \
        sum(integrals_R) - f*v*dx + sum(integrals_N)
    a, L = lhs(F), rhs(F)

    # Set linear solver parameters
    prm = LinearVariationalSolver.default_parameters()
    if linear_solver == 'Krylov':
        prm["linear_solver"] = 'gmres'
        prm["preconditioner"] = 'ilu'
        prm["krylov_solver.absolute_tolerance"] = abs_tol
        prm["krylov_solver.relative_tolerance"] = rel_tol
        prm["krylov_solver.maximum_iterations"] = max_iter
    else:
        prm["linear_solver"] = 'lu'

    # Compute solution
    u = Function(V)
    solve(a == L, u, bcs, solver_parameters=prm)

    return u


def demo_bcs():
    "Compute and plot solution using a combination of boundary conditions"

    # Define manufactured solution in sympy and derive f, g, etc.
    import sympy as sym
    x, y = sym.symbols('x[0], x[1]')            # needed by UFL
    u = 1 + x**2 + 2*y**2                       # exact solution
    u_e = u                                     # exact solution
    u_00 = u.subs(x, 0)                         # restrict to x = 0
    u_01 = u.subs(x, 1)                         # restrict to x = 1
    f = -sym.diff(u, x, 2) - sym.diff(u, y, 2)  # -Laplace(u)
    f = sym.simplify(f)                         # simplify f
    g = -sym.diff(u, y).subs(y, 1)              # compute g = -du/dn
    r = 1000                                    # Robin data, arbitrary
    s = u                                       # Robin data, u = s

    # Collect variables
    variables = [u_e, u_00, u_01, f, g, r, s]

    # Turn into C/C++ code strings
    variables = [sym.printing.ccode(var) for var in variables]

    # Turn into FEniCS Expressions
    variables = [Expression(var, degree=2) for var in variables]

    # Extract variables
    u_e, u_00, u_01, f, g, r, s = variables

    # Define boundary conditions
    boundary_conditions = {0: {'Dirichlet': u_00},   # x = 0
                           1: {'Dirichlet': u_01},   # x = 1
                           2: {'Robin':     (r, s)}, # y = 0
                           3: {'Neumann':   g}}      # y = 1

    # Compute solution
    kappa = Constant(1)
    Nx = Ny = 8
    u = solver_bcs(kappa, f, boundary_conditions, Nx, Ny,
                   degree=1, linear_solver='direct')

    # Compute maximum error at vertices
    mesh = u.function_space().mesh()
    vertex_values_u_e = u_e.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = np.max(np.abs(vertex_values_u_e -
                              vertex_values_u))
    print('error_max =', error_max)
    print(errornorm(u_e, u))

    # Save and plot solution
    #vtkfile = File('poisson_extended/solution_bcs.pvd')
    #vtkfile << u
    plot(u)

demo_bcs()

plt.show()