from fenics import *
import numpy as np 
import time 


from utils import smoothmax
from utils import unity_partition_1d
from utils import unity_partition_2d
from utils import two_mountains
from utils import obs_1d


def solver_op(example, n, T=1, num_steps=100):
    
    dt = T/num_steps
    
    if example=='example_1':
        list_u=np.empty((0, n+1))    
        mesh = IntervalMesh(n, -1, 1)
    elif example=='example_2':
        list_u=np.empty((0, (n+1)*(n+1)))    
        mesh = RectangleMesh(Point(-1, -1), Point(1,1), n, n)

    Vh = FunctionSpace(mesh, "CG", 1)

    u = Function(Vh)    
    u0 = Constant(0)
    
    v = TestFunction(Vh)
    if example=='example_1':
        psi = interpolate(bump2(degree=2), Vh)

    elif example=='example_2':
        psi = interpolate(two_mountains(degree=2),Vh)
        
    un = interpolate(u0, Vh)
    
    f = Constant(-3)
    
    eps = Constant(pow(10, -8))
    F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (dt/eps)*inner(smoothmax(-u+psi), v)*dx - un*v*dx - dt*f*v*dx
    bc = DirichletBC(Vh, Constant(0), "on_boundary")

    t = 0
    vtksol = File("output/popsol.pvd")
    vtkobs = File("output/popobs.pvd")
    #vtkel = File("output/el.pvd")
    vtkobs << psi
    
    xyz = Vh.tabulate_dof_coordinates()
    x = xyz[:,0]
    
    if example=='example_1':
        y=0
    elif example=='example_2':
        y = xyz[:,1]
    for n in range(num_steps):
        t+=dt
        solve(F==0, u, bcs=bc)
        vtksol << (u, t)
        un.assign(u)         
        list_u = np.append(list_u, [u.vector().get_local()],axis=0)

    return list_u, psi.vector().get_local(), x, y

def sol_one_realization(batches, example, n=100, T=1, num_steps=100):
    
    if example=='example_1':
        list_w=np.empty((0, n+1))
    elif example=='example_2':
        list_w=np.empty((0, (n+1)*(n+1)))

    
    
    dt = T/num_steps
    
    if example=='example_1':
        mesh = IntervalMesh(n, -1, 1)
    elif example=='example_2':
        mesh = RectangleMesh(Point(-1, -1), Point(1,1), n, n)

    Vh = FunctionSpace(mesh, "CG", 1)
    w = Function(Vh)
    

    u0 = Constant(0)
    
    v = TestFunction(Vh)
   
    #psi = Constant(0.)
    if example=='example_1':
        psi = interpolate(bump2(degree=2), Vh)
        chi = unity_partition_1d(n=0, element=Vh.ufl_element())
        pi = 2

    elif example=='example_2':
        psi = interpolate(two_mountains(degree=2),Vh)
        chi = unity_partition_2d(n=0, element=Vh.ufl_element())
        pi = 4

#    chi = unity_partition_2d(n=0, element=Vh.ufl_element())
   
    
    wn = interpolate(u0, Vh)
    
    f = Constant(-3)
    
    eps = Constant(pow(10, -8))    
    F_R = w*v*dx + dt*dot(chi*pi*grad(w), grad(v))*dx - (dt/eps)*inner(smoothmax(-w+psi), v)*dx - wn*v*dx - dt*chi*pi*f*v*dx

    bc = DirichletBC(Vh, Constant(0), "on_boundary")

    t = 0
    vtksol_ran = File("output/popsol_ran.pvd")

    vtkobs = File("output/popobs.pvd")
    #vtkel = File("output/el.pvd")
    vtkobs << psi
    
    xyz = Vh.tabulate_dof_coordinates()
    x = xyz[:,0]
    if example=='example_1':
        y=0
    elif example=='example_2':
        y = xyz[:,1]

    for n in range(num_steps):
        t+=dt
        solve(F_R==0, w, bcs=bc)
        vtksol_ran << (w, t)
        wn.assign(w) 
        
        chi.n = batches[n]

        list_w = np.append(list_w, [w.vector().get_local()],axis=0)
    return list_w, psi.vector().get_local(), x, y


def solver_avg_mb_op(bat,nx,nt,T,example,itera):

    if example=='example_1':
        avr=np.zeros((nt,nx+1))
    elif example=='example_2':
        avr=np.zeros((nt,(nx+1)*(nx+1)))
    
    avr_time=0    
    for i in range(itera):
        print(f'Realization number: {i}')
        start_time = time.time()
        wn_list, obstacle,X,Y = sol_one_realization(bat[i],example,nx, T, nt)
        avr+=wn_list
        avr_time+=time.time() - start_time
    return avr/itera, avr_time/itera

