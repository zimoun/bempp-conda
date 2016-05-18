# coding: utf8

import bempp.api

bempp.api.test()

print(bempp.api.__path__)
print(bempp.api._config.version)
print(bempp.api.__version__)
# coding: utf8

print('running laplace.py')

# # Solving a Laplace problem with Dirichlet boundary conditions

# ### Background

# In this tutorial we will solve a simple Laplace problem inside the unit sphere $\Omega$ with Dirichlet boundary conditions. The PDE is given by
#
# $$
# \Delta u = 0
# $$
#
# in $\Omega$ with boundary conditions
#
# $$
# u = g
# $$
# on the boundary $\Gamma$ of $\Omega$. The boundary data is a source $\hat{u}$ located at the point $(.9,0,0)$.
# $$
# \hat{u}(\mathbf x)=\frac{1}{4\pi\sqrt{(x-.9)^2+y^2+z^2}}.
# $$

# For this example we will use an direct integral equation of the first kind. Let
# $$
# g(\mathbf x,\mathbf y) = \frac{1}{4\pi |\mathbf x-\mathbf y|}
# $$
# the Green's function in three dimensions with $|\mathbf x|^2=x^2+y^2+z^2$. Then from Green's representation theorem it follows that every function $u$ harmonic in $\Omega$ satisfies
#
# $$
# u(\mathbf x) = \int_{\Gamma} g(\mathbf x,\mathbf y)\frac{\partial u(\mathbf y)}{\partial n(\mathbf{y})}ds(\mathbf y)-\int_{\Gamma}\frac{\partial g(\mathbf x,\mathbf y)}{\partial n(\mathbf{y})}u(\mathbf y)ds(\mathbf y),~\mathbf x\in\Omega.
# $$
#
# Taking the limit $\mathbf x\rightarrow \Gamma$ we obtain the boundary integral equation
#
# $$
# \left[V\frac{\partial u}{\partial n}\right](\mathbf x)=\frac12 u(\mathbf x)+\left[Ku\right](\mathbf x),~\mathbf x\in\Gamma.
# $$
#
# Here, $V$ and $K$ are the single and double-layer potential boundary operators defined by
#
# $$
# \begin{align}
# \left[V\phi\right](\mathbf x)&=\int_{\Gamma}g(\mathbf x,\mathbf y)\phi(\mathbf y)ds(y)\\
# \left[K\phi\right](\mathbf x)&=\int_{\Gamma}\frac{\partial g(\mathbf x,\mathbf y)}{\partial n(\mathbf{y})}\phi(\mathbf y)ds(\mathbf y)
# \end{align}
# $$
#
# for $x\in\Gamma$.
#

# ### Implementation

# In the following we demonstrate how to solve this problem with BEM++. We first define the known Dirichlet boundary data. In this example we will use a Python function for it. Other ways are possible (such as a vector of coefficients at the nodes of a mesh).

# In[1]:

import bempp.api
import numpy as np

def dirichlet_data(x, n, domain_index, result):
    result[0] = 1./(4 * np.pi * ((x[0] - .9)**2 + x[1]**2 + x[2]**2)**(0.5))


# A valid Python function to define a BEM++ GridFunction takes the inputs `x`,`n`,`domain_index` and `result`. `x` is a three dimensional coordinate vector. `n` is the normal direction. The `domain_index` allows to identify different parts of a physical mesh in order to specify different functions on different subdomains. `result` is a Numpy array that will store the result of the function call. For scalar problems it just has one component `result[0]`.

# We now define a mesh or grid in BEM++ notation. Normally one reads a grid from a file. BEM++ supports import and export to Gmsh with other data formats to follow soon. However, for this problem we do not need a complicated mesh but will rather use the built-in function `sphere` that defines a simple spherical grid.

# In[2]:

grid = bempp.api.shapes.sphere(h=0.1)


# We now define the spaces. For this example we will use two spaces, the space of continuous, piecewise linear functions and the space of piecewise constant functions. The space of piecewise constant functions has the right smoothness for the unknown Neumann data. We will use continuous, piecewise linear functions to represent the known Dirichlet data.

# In[3]:

piecewise_const_space = bempp.api.function_space(grid, "DP", 0) # A disccontinuous polynomial space of order 0
piecewise_lin_space = bempp.api.function_space(grid, "P", 1)    # A continuous piecewise polynomial space of order 1


# We can now define the operators. We need the identity operator, and the single-layer, respectively double-layer, boundary operator.  The general calling convention for an operator is
#
#     op = factory_function(domain_space,range_space,dual_to_range_space,...)
#
# Typically, for a Galerkin discretisation only the domain space and the dual space (or test space) are needed. BEM++ also requires a notion of the range of the operator. This makes it possible to define operator algebras in BEM++ that can be used almost as if the operators are continuous objects.

# In[4]:

identity = bempp.api.operators.boundary.sparse.identity(
    piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
dlp = bempp.api.operators.boundary.laplace.double_layer(
    piecewise_lin_space, piecewise_lin_space, piecewise_const_space)
slp = bempp.api.operators.boundary.laplace.single_layer(
    piecewise_const_space, piecewise_lin_space, piecewise_const_space)


# We now define the GridFunction object on the sphere grid that represents the Dirichlet data.

# In[5]:

dirichlet_fun = bempp.api.GridFunction(piecewise_lin_space, fun=dirichlet_data)


# The below code will assemble the identity and double-layer boundary operator and evaluate the right-hand side of the boundary integral equation. This is an exact analogue of the underlying mathematical formulation. Depending on the grid size this command can take a bit since here the actual operators are assembled. The left-hand side only consists of the single-layer potential operator in this example. This is here not yet assembled as it is not yet needed. In BEM++ operators are only assembled once they are needed.

# In[6]:

rhs = (.5*identity+dlp)*dirichlet_fun
lhs = slp


# The following code solves the boundary integral equation iteratively using Conjugate Gradients. BEM++ offers a CG and GMRES algorithm. Internally these are just simple interfaces to the corresponding SciPy functions with the difference that the BEM++ variants accept BEM++ operators and GridFunctions as objects instead of just operators and vectors.

# In[7]:

neumann_fun, info = bempp.api.linalg.cg(slp, rhs, tol=1E-3)


# We now want to provide a simple plot of the solution in the (x,y) plane for z=0. First we need to define points at which to plot the solution.

# In[8]:

n_grid_points = 150
plot_grid = np.mgrid[-1:1:n_grid_points*1j,-1:1:n_grid_points*1j]
points = np.vstack((plot_grid[0].ravel(),plot_grid[1].ravel(),np.zeros(plot_grid[0].size)))


# The variable `points` now contains in its columns the coordinates of the evaluation points. We can now use Green's representation theorem to evaluate the solution on these points. Note in particular the last line of the following code. It is a direct implementation of Green's representation theorem.

# In[9]:

slp_pot = bempp.api.operators.potential.laplace.single_layer(piecewise_const_space,points)
dlp_pot = bempp.api.operators.potential.laplace.double_layer(piecewise_lin_space,points)
u_evaluated = slp_pot*neumann_fun-dlp_pot*dirichlet_fun


# We now want to create a nice plot from the computed data. We only plot a slice through $z=0$. For a full three dimensional visualization BEM++ allows to export data to Gmsh. Since the solution decays quickly we will use a logarithmic plot.

# In[10]:

# The next command ensures that plots are shown within the IPython notebook
# get_ipython().magic('matplotlib inline')

# Filter out solution values that are associated with points outside the unit circle.
u_evaluated = u_evaluated.reshape((n_grid_points,n_grid_points))
radius = np.sqrt(plot_grid[0]**2+plot_grid[1]**2)
u_evaluated[radius>1] = np.nan

print('  -> bempp computations done.')

# Plot the image
try:
    import matplotlib
    matplotlib.rcParams['figure.figsize'] = (5.0, 4.0) # Adjust the figure size in IPython

    from matplotlib import pyplot as plt

    fig = plt.figure()
    plt.imshow(np.log(np.abs(u_evaluated.T)), extent=(-1,1,-1,1),origin='lower')
    plt.title('Computed solution')
    #plt.show(block=False)
    fig.savefig('laplace.png',  bbox_inches='tight')

except:
    print('something wrong with Matplotlib')

print('\n passed. DONE \n')

# coding: utf-8

print('running scattering.py')

# # Scattering from a sphere using a combined direct formulation

# ### Background

# In this tutorial we will solve the problem of scattering from the unit sphere $\Omega$ using a combined integral formulation and an incident wave:
#
# $$
# u^{\text{inc}}(\mathbf x) = e^{i k x},
# $$
#
# where $\mathbf x = (x, y, z)^t$.
#
# The PDE is given by the Helmholtz equation:
#
# $$
# \Delta u + k^2 u = 0, \quad \text{ in } \mathbb{R}^3 \backslash \Omega,
# $$
#
# where $u=u_s+u_{inc}$ is the total acoustic field and $u_{s}$ satisfies the Sommerfeld radiation condition
#
# $$
# \frac{\partial u_s}{\partial r}-iku_s=o(r^{-1})
# $$
#
# for $r:=\frac{\mathbf x}{|x|}\rightarrow\infty$.
#
# From Green's representation formula one can derive that
#
# $$
# u(\mathbf x) = u_{inc}-\int_{\Gamma}g(\mathbf x,\mathbf y)u_n(\mathbf y)ds(y).
# $$
#
# Here, $g(\mathbf x, \mathbf y)$ is the acoustic Green's function given as
#
# $$
# g(\mathbf x, \mathbf y):=\frac{e^{i k \|\mathbf{x}-\mathbf{y}\|}}{4 \pi \|\mathbf{x}-\mathbf{y}\|}.
# $$
#
# The problem has therefore been reduced to computing the normal derivative $u_n$ on the boundary $\Gamma$. This is achieved through the following boundary integral equation formulation.
#
# $$
# (\frac12I + D_k' - i \eta S_k) u_n = \frac{\partial u^{\text{inc}}}{\partial \nu}(x) - i \eta u^{\text{inc}}(x), \quad x \in \Gamma.
# $$
#
# where $I,D_k'$ and $S_k$ are respectively the identity operator, the adjoint double layer boundary operator and the single layer boundary operator. More details of the derivation of this formulation and its properties can be found in the article [*S. N. Chandler-Wilde, I. G. Graham, S. Langdon and E. A. Spence, Numerical-asymptotic boundary integral methods in high frequency acoustic scattering*](http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=8539370&fileId=S0962492912000037).
#

# ### Implementation

# First we import the bempp module and numpy.

# In[1]:

import bempp.api
import numpy as np


# We define the wavenumber

# In[2]:

k = 15.


# The rhs of the combined formulation is defined as follows.

# In[3]:

def combined_data(x, n, domain_index, result):
    result[0] = 1j * k * np.exp(1j * k * x[0]) * (n[0]-1)


# The following command creates a sphere mesh.

# In[4]:

grid = bempp.api.shapes.regular_sphere(5)


# As basis functions we use piecewise constant functions over the elements of the mesh. The corresponding space is initialized as follows.

# In[5]:

piecewise_const_space = bempp.api.function_space(grid, "DP", 0)


# We now initialize the boundary operators.
# A boundary operator always takes at least three space arguments: a domain space, a range space and the test space (dual to the range). In this example we only work on the space $L^2(\Gamma)$ and we can choose all spaces to be identical.

# In[6]:

identity = bempp.api.operators.boundary.sparse.identity(
    piecewise_const_space, piecewise_const_space, piecewise_const_space)
adlp = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(
    piecewise_const_space, piecewise_const_space, piecewise_const_space,k)
slp = bempp.api.operators.boundary.helmholtz.single_layer(
    piecewise_const_space, piecewise_const_space, piecewise_const_space,k)


# Standard arithmetic operators can be used to create linear combinations of boundary operators.

# In[7]:

lhs = .5*identity + adlp - 1j * k * slp


# We now discretize the right-hand side. Here, this means converting the Python callable `data_fun` into a GridFunction defined on the elements of the grid.

# In[8]:

grid_fun = bempp.api.GridFunction(piecewise_const_space, fun=combined_data)


# We can now use GMRES to solve the problem.

# In[9]:

from bempp.api.linalg import gmres
neumann_fun,info = gmres(lhs, grid_fun, tol=1E-5)


# `gmres` returns a grid function `neumann_fun` and an integer `info`. When everything works fine info is equal to 0.

# At this stage, we have the surface solution of the integral equation. Now we will evaluate the solution in the domain of interest. We define the evaluation points as follows.

# In[10]:

Nx = 200
Ny = 200
xmin, xmax, ymin, ymax = [-3, 3, -3, 3]
plot_grid = np.mgrid[xmin:xmax:Nx * 1j, ymin:ymax:Ny * 1j]
points = np.vstack((plot_grid[0].ravel(), plot_grid[1].ravel(), np.zeros(plot_grid[0].size)))
u_evaluated = np.zeros(points.shape[1], dtype=np.complex128)
u_evaluated[:] = np.nan


# Then we create a single layer potential operator and use it to evaluate the solution at the evaluation points.

# In[11]:

x, y, z =points
idx = np.sqrt(x**2 + y**2) > 1.0


# The variable idx allows to compute only points located outside the unit circle of the plane. We use a single layer potential operator to evaluate the solution at the observation points.

# In[12]:

from bempp.api.operators.potential import helmholtz as helmholtz_potential
slp_pot=helmholtz_potential.single_layer(piecewise_const_space,points[:,idx],k)
res = np.real(np.exp(1j *k * points[0,idx]) - slp_pot.evaluate(neumann_fun))
u_evaluated[idx] = res.flat


# We can now easily plot a slice of the domain solution.

# In[13]:

u_evaluated=u_evaluated.reshape((Nx, Ny))

print('  -> bempp computations done.')

# get_ipython().magic('matplotlib inline')
# Plot the image
try:
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize =(10, 8))
    plt.imshow(np.real(u_evaluated.T),extent=[-3, 3, -3, 3])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar()
    plt.title("Scattering from the unit sphere, solution in plane z=0")
    #plt.show(block=False)
    fig.savefig('scattering.png',  bbox_inches='tight')
except:
    print('something wrong with Matplotlib')

print('\n passed. DONE \n')
# coding: utf8

print('running maxwell.py')

# # Electromagnetic scattering from a screen

# ### Background

# In this tutorial we consider the scattering of an electromagnetic wave from a perfectly conducting screen $\Gamma:=[-2,2]\times[-1,1]\times0$. The time-harmonic Maxwell equation for the electric field $E$ reduces to
#
# $$
# \nabla\times\nabla\times \mathbf{E} -k^2 \mathbf{E} = 0
# $$
#
# in $\mathbb{R}^3\backslash\Gamma$, where $k:=2\pi/\lambda$ is the wavenumber and $\lambda$ is the wavelength. The electric field $\mathbf{E}$ is the sum of the incident field $\mathbf{E}^{(i)}$ and the scattered field $\mathbf{E}^{(s)}$. Here, we assume that the incident field is given by
#
# $$
# \mathbf{E}^{(i)}:=\begin{bmatrix} e^{ikz} & 0 & 0 \end{bmatrix},
# $$
# which is a wave travelling in the $z$ direction and polarised in the $x$ direction. On the screen the tangential component $\mathbf{E}_t:=n\times \mathbf{E}$ must be zero. Towards infinity we impose the Silver-Muller radiation condition
#
# $$
# \lim_{|\mathbf{x}|\rightarrow\infty} |\mathbf{x}|\left(\nabla\times \mathbf{E}^{(s)}\times\hat{\mathbf{x}}-ik\mathbf{E}^{(s)}\right) = 0,
# $$
# where $\hat{\mathbf{x}}=\mathbf{x}/|\mathbf{x}|$.
#
# The scattered wave $\mathbf{E}^{(s)}$ has the representation
#
# $$
# \mathbf{E}^{(s)} = -\Psi_{\mathbf{SL}}(\mathbf{\lambda}),
# $$
#
# where $\lambda$ is the jump of the Neumann trace of the scattered field $\mathbf{E}^{(s)}$ across the screen. The Maxwell electric field potential operator $\Psi_{\mathbf{SL}}$ is defined as
#
# $$
# \Psi_{\mathbf{SL}}(\mathbf{v}):=ik\int_{\Gamma}G(\mathbf{x},\mathbf{y})\mathbf{v}(\mathbf{y})d\Gamma(\mathbf{y})-
# \frac{1}{ik}\nabla_{\mathbf{x}}\int_{\Gamma}G(\mathbf{x},\mathbf{y})(\nabla_{\Gamma}\cdot\mathbf{v})(\mathbf{y})d\Gamma(\mathbf{y})
# $$
# with $G(\mathbf{x},\mathbf{y}):=\frac{e^{ik|\mathbf{x}-\mathbf{y}|}}{4\pi|\mathbf{x}-\mathbf{y}|}$.
#
# The associated boundary operator is denoted by $\mathcal{S}$. It is obtained as average from the tangential traces of the electric field potential operator from both sides of the screen. The boundary integral equation is now
#
# $$
# \mathcal{S}\mathbf{\lambda} = n\times E^{(i)}.
# $$
# The $-$ sign is missing in comparison to the representation formula since we want to satisfy the boundary conditions for the negative incident wave so that the tangential trace of the total field is zero on the screen.
#
# More details about the mathematical background can be found in the wonderful overview paper by [Buffa and Hiptmair (2003)](http://www.sam.math.ethz.ch/~hiptmair/Courses/CEM/BUH03.pdf).

# ### Implementation

# We start with the usual imports.

# In[1]:

import bempp.api
import numpy as np


# To avoid preconditioning issues we assemble the operators in dense mode so that we can solve via LU decomposition later on.

# In[2]:

bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'


# We need to define the wavenumber of the problem.

# In[3]:

wavelength = .5
k = 2 * np.pi / wavelength
#k = 1
#wavelength = 2 * np.pi / k

# The function `incident_field` defines the incident wave. The function `tangential_trace` computes the tangential trace of the incident field.

# In[4]:

def incident_field(x):
    return np.array([np.exp(1j * k * x[2]), 0. * x[2], 0. * x[2]])

def tangential_trace(x, n, domain_index, result):
    result[:] = np.cross(incident_field(x), n, axis=0)


# We define a structured grid in the x-y plane.

# In[5]:

grid = bempp.api.structured_grid((-2, -1), (2, 1), (60, 30))
# grid = bempp.api.shapes.cube(length=2, origin=(0, 0, 0), h=0.1)


# A suitable space for Maxwell problems is the Raviart-Thomas space of order 0. Higher-order Raviart-Thomas spaces are currently not available in BEM++.

# In[6]:

space = bempp.api.function_space(grid, "RT", 0)


# In the following we define the Maxwell single-layer boundary operator and the Maxwell identity operator. For Maxwell problems the suitable inner product is not the standard $L^2$ inner product but the anti-symmetric pairing $\langle \mathbf{\mu},\mathbf{\eta}\rangle_{\mathbf{\tau},\mathbf{\Gamma}}:=\int_{\Gamma}\overline{\mathbf{\mu}(\mathbf{x})}\cdot (\mathbf{\eta}(\mathbf{x})\times n(\mathbf{x}))d\Gamma(\mathbf{x})$. Hence, the standard identity is not suitable and a specific Maxwell identity operator is defined in BEM++, which implements this pairing.

# In[7]:

slp = bempp.api.operators.boundary.maxwell.electric_field(space, space, space, k)
identity = bempp.api.operators.boundary.sparse.maxwell_identity(space, space, space)


# The following command creates a grid function from the tangential trace.

# In[8]:

trace_fun = bempp.api.GridFunction(space, fun=tangential_trace)


# We have to multiply the GridFunction `trace_data` with the Maxwell mass matrix `id` in order to obtain the anti-symmetric dual pairing of the basis functions with the incident wave.

# In[9]:

rhs = identity * trace_fun


# We use a direct solver to solve the system. It is not a large problem and Krylov methods without preconditioning converge poorly for the electric field integral equation.

# In[10]:

from bempp.api.linalg import lu
lambda_data = lu(slp, rhs)


# Now that the solution $\mathbf{\lambda}$ is computed we want to plot the total field. First, we define a grid of points in the x-z plane.

# In[11]:

# Create a grid of points

nx = 151
nz = 151
extent = 5
x, y, z = np.mgrid[-extent:extent:nx * 1j, 0:0:1j, -extent:extent:nz * 1j]
points = np.vstack((x.ravel(), y.ravel(), z.ravel()))


# We now initialize the single-layer potential operator. By default potentials are assembled using accelerated H-Matrix techniques. Here, we set the accuracy of the assembly to 1E-2. This is quite low but sufficient for plotting accuracy.

# In[12]:

bempp.api.global_parameters.hmat.eps = 1E-2
slp_pot = bempp.api.operators.potential.maxwell.electric_field(space, points, k)


# The following commands now compute the total field by first computing the scattered field from the representation formula and then summing into it the incident field.

# In[13]:

scattered_field_data = -slp_pot * lambda_data
incident_field_data = incident_field(points)
field_data = scattered_field_data + incident_field_data


# In electromagnetic scattering it is often useful to visualize the squared electric field density. This value is computed below.

# In[14]:

squared_field_density = np.real(np.sum(field_data * field_data.conj(), axis=0))


# Finally, we can plot everything using a simple Matplotlib plot.

print('  -> bempp computations done.')

# In[15]:

# get_ipython().magic('matplotlib inline')
try:
    import matplotlib
    from matplotlib import pyplot as plt
    # Adjust the figure size in IPython
    matplotlib.rcParams['figure.figsize'] = (10.0, 8.0)

    fig = plt.figure()
    plt.imshow(squared_field_density.reshape((nx,nz)).T,
               cmap='coolwarm', origin='lower',
               extent=[-extent, extent, -extent,extent])
    plt.colorbar()
    plt.title("Squared Electric Field Density")
    #plt.show(block=False)
    fig.savefig('maxwell.png',  bbox_inches='tight')
except:
    print('something wrong with Matplotlib')

print('\n passed. DONE \n')
