using JLD          # To save data

# Import functions and parameters
include("include/function_master.jl")
include("parameters_input.jl")

# Create folder to store current simulation
savePath = "data/"*simulationName
run(`mkdir -p $savePath`)
run(`cp parameters_input.jl $savePath/parameters_input.jl`)
# Copy plotting script to simulation folder
run(`cp include/plot.jl $savePath/plot.jl`)
println("Folder $savePath is ready.")

# Define global needed quantities
Ek = zeros(Nx, ntime +1);
EElect = zeros(ntime +1);
EKin = copy(EElect);
f = zeros(insToSave+1, Nx, Nv);

# Initialize perturbed distribution function at t=0
f0 = f_init(x, v, eps, kmode)

# Ensure periodic boundary conditions and take distribution to fourier space
f0[:,1:3] = f0[:,Nv-3:Nv] = 0.0;   # kills the distribution at the boundaries
f0 = fourierspace2d(f0);          # fxv --> fku


# Solve electrostatic 1D Vlasov equation
println("Entering main loop...")
tic()
fku, Ek, EElect, EKin = CHIN_A_INT(dt, ntime, x, dv, v, f0, insToSave);
toc()

# Bring quantities to real space
for i in 1:insToSave+1;
    f[i,:,:] = realspace2d(fku[i,:end,:end]);
end
E = realspace1d(Ek);

# Save data
println("Now saving calculated quantities...")
save(savePath*"/data.jld","f", f, "E", E, "EElect", EElect, "EKin", EKin)
