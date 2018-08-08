using JLD
using PyPlot

include("../../include/functions.jl")
include("parameters_input.jl")

# Loading calculated quantities
ft, EField, EKin = load("data.jld","ft", "EField", "EKin")

run(`mkdir -p plot`)
cd("plot")

taxis = [i*dt for i=0:ntime]
itsPerSave = floor(Int32, dt*ntime/insToSave)
ipw = floor(Int32, 1/dt)  # Ammount of time steps to make up 1 electron plasma period.
wPerPoint = 1 # Plot every wPerPoint plasma periods

ipw = ipw*wPerPoint # No not modify

# Other quantities of interest
EEField = ElectrostaticEnergy(dx, EField); # Total electrostatic energy of the electric field
EPow = EPowerPerMode(EField); # Power of EField for each mode. EPow[time:space]


######################## Plots ########################
figur=0

# Power spectrum of EField
figure(figur)
figur=figur+1

maxmode = Nx
#max_2 = Int(maxmode/2)
mode_zero = 1
modesToPlot = 4
for mode = mode_zero+modesToPlot:-1:mode_zero+1
    plot(taxis[1:ipw:end], log10.(EPow[1:ipw:end,mode]), linestyle="-", marker="", linewidth=0.4, label=string("mode=",mode-mode_zero));
    legend(loc="best");
    ylabel(string("log(EField Energy) per mode"));
    xlabel("Time, \$\\omega_{pe}^{-1}\$");
end
grid()
savefig("EPowPerMode.pdf")


# Electric field energy
figure(figur)
figur=figur+1

plot(taxis[1:ipw:end], log10.(EEField[1:ipw:end]), linestyle="-", marker="");
title("Potential electric energy")
grid()
ylabel("log(Field Energy)");
xlabel("Time, \$\\omega_{pe}^{-1}\$");
savefig("EEField.pdf")

# Kinetic energy
figure(figur)
figur=figur+1

plot(taxis[1:ipw:end],EKin[1:ipw:end], linestyle="-", marker="", color="red");
title("Kinetic energy")
xlabel("Time, \$\\omega_{pe}^{-1}\$")
ylabel("Kinetic energy")
grid()
savefig("EKin.pdf")

# Total energy (EEField + EKin)
figure(figur)
figur=figur+1

Etotal = EEField + EKin;
emax = 1.00001*maximum(Etotal[ipw:end])
emin = 0.99999*minimum(Etotal[ipw:end])

plot(taxis[ipw:ipw:end],Etotal[5:ipw:end], linestyle="-", marker="", color="red");
grid()
title("Total energy")
axis([taxis[ipw], taxis[end], emin, emax]); # getting rid of first plasma period
xlabel("Total energy, \$ E_p + E_k \$")
xlabel("Time, \$\\omega_{pe}^{-1}\$")
savefig("ETot.pdf")


# Level curves (filled)

for inst in 1:(insToSave+1)

    figure(figur)
    figur=figur+1
    
    contourf(x, v, ft[inst,:,:]', cmap=ColorMap("jet"),54);#, vmin = 0, vmax = maximum(ft[inst,:,:]));
    colorbar()
    xlabel("Position")
    ylabel("Velocity")
    title("Distribution function in x and v, \$ t="*string((inst-1)*itsPerSave)*" \\omega_{pe}^{-1} \$")
    grid()
    axis([0,L,vmin,vmax])
    savefig("levelcurves_f_"*string((inst-1)*itsPerSave)*".pdf")
    close()
end


# Mean dist. function vs v
figure(figur)
figur=figur+1

for inst=1:length(ft[:,end,end])
    plot(v, sum(ft[inst,:,:],1)'./Nx, linestyle="-", marker="", linewidth=0.4, label="f_"*string((inst-1)*itsPerSave));
end
title("Mean distribution function")
xlabel("Velocity")
grid()
legend(loc="best");
savefig("f_vs_v.pdf")

#
figure(figur)
figur=figur+1

contourf(taxis[1:ipw:end], fftshift(kx), fftshift(EPow,2)[1:ipw:end,:]', cmap=ColorMap("jet"), 200);
colorbar()
title("E power in t and k")
grid()
xlabel("time, t")
ylabel("Wavenumber, k")
axis([taxis[1],taxis[end], -1, 1])
savefig("EPow_contour.pdf")

#
println("All graphics have been plotted.")
