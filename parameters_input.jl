##################################################################################
#                                     W A R N I N G                              #
#                 Thy must not exceed the recurrence time T_rec                  #
#                                                                                #

Nx = 512            # Nodes on space
Nv = 2048           # Nodes on velocity space
dt = 0.01           # Time step
ntime = 100*1800    # Iteration number for time-advance

L = 5*pi                       # Length of box
mode = 1                        # Perturbation mode
kmode = 2pi*mode/L              # Perturbation wavenumber
eps = 5.e-2

vmin = 0
vmax = 6

insToSave = 36                   # Ammount of nonzero instants to save

simulationName = "k0_4-m1-1800wp"      # Makes folder with that name to save current data


# Distribution function at time t = 0

function f_init(x, v, eps, kmode)
    # Set the initial conditions of the simulation. (Initial distribution function)
    # Inputs, x: positions vector, v: velocity vector, 
    #         eps: perturbation amplitude, kmode: perturbation wavenumber.
    
    Nx = length(x)

    f = zeros(Nx,Nv)
    
    f0 = (exp.(-(v.^2)./2)')./sqrt(2pi)       # Valentini et al. 2004

    f = ( 1 + eps*cos.(x*kmode) )*f0
    return f
end

#                                                                                #
#                                                                                #
##################################################################################

# Calculate derivated quantities
dx = L/Nx
x = [i*dx for i in 0:Nx]
dx = L/(Nx-1)

vrange = vmax - vmin
dv = vrange/Nv
v = [vmin + i*dv for i in 0:Nv]
dv = vrange/(Nv-1)

# end
println("Parameters loaded successfully.")
