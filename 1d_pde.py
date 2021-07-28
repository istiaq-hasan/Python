import numpy                    #here we load numpy library that provides a bunch of useful matrix operations akin to MATLAB
from matplotlib import pyplot   #2D plotting library that we will use to plot our results

# lineSingle = '------------------------------------------------'

# print("Solving 1D Diffusion Equation using Finite Difference Method\n")

#setting the grid

nx = 500                     #grid points
# dx  = 2 / (nx - 1)          #grid spacing
dx = 0.1
nt = 5000                     #number of timesteps
nu = 0.9                    #viscosity
cfl = 0.4                   
# dt = cfl*dx**2/nu           #based on von neumaan stability analysis
dt = 0.001

x = numpy.arange(0, (nx + 1) * dx, dx)
t = numpy.arange(0, (nt + 1) * dt, dt)


#innitial condition

# print(lineSingle)
# print("Computing Innitial Solution...")


u = numpy.ones(len(x))
u[int(nx / 4):int(3 * nx / 4)] = 2      #Square Wave Profile

# u = numpy.random.rand(nx)           #random value

# print("Printing Innitial Solution...")
# print(lineSingle)

# print(u)

pyplot.plot(x, u, label='Initial Solution')

#discritization

# print(lineSingle)
# print("Calculating Numerical Solution......")
# print(lineSingle)

un = numpy.ones(len(x))  

for n in range(nt+1):               #time marching
    un = u.copy()       
    for i in range(1, nx-1):        #Space marching
        
        u[i] = un[i] + dt * (un[i] * (1 - un[i])) + dt/dx**2 *(un[i+1] - 2*un[i] + un[i-1]) #Central Differnece Scheme
        
        u[0] = un[0] + dt * (un[0] * (1 - un[0])) + dt/dx**2 * 2 * (un[1] - un[0])
        u[nx - 1] = un[nx - 1] + dt * (un[nx - 1] * (1 - un[nx - 1])) + dt/dx**2 * 2 * (un[nx - 2] - un[nx - 1])
    if n > 0 and n%1000 == 0:
        pyplot.plot(x, u, label='Numerical Solution at time ' + str(t[n]))
        

# print(lineSingle)
# print("Printing Numerical Solution......")
# print(lineSingle)
        
# # print(u)

# print(lineSingle)
# print("Plotting Innitial & Numerical Solution")
# print(lineSingle)

pyplot.plot(x, u, label='Numerical Solution' + str(t[-1]))
pyplot.title('1D Diffusion Convecction')
pyplot.xlabel('Grid Space')
pyplot.ylabel('Velocity')
pyplot.grid()
pyplot.legend()
pyplot.savefig('data.png')
pyplot.show()