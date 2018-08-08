# Chin A integrator for vlasov_1d

function CHIN_A_INT( dt::Float64, ntime::Int64,  x::Array{Float64}, dv::Float64,
                     v::Array{Float64}, fku::Array{Complex{Float64}, 2}, insToSave::Int64);

    ##  NOTE: _ fku = fku(k,u), where u is v-wavenumber
    ##        _ Ek_out = Ek_out(t,k)
    ##        _ rhok = rhok(k)
    ##        _ E = E(x)
    ##        _ rho = rho(x)
    
    Nx = size(fku,1);
    Nv = size(fku,2);
    vrange = v[end]-v[1];

    L = x[end] - x[1];
    k = wavevector(x);
    u = wavevector(v);
    du  = 2pi/(vrange*(Nv-1));     
    
    f_out = zeros(Complex{Float64}, insToSave, Nx, Nv)
    f_out[1,:,:] = fku;
    
    EElect = zeros(Float64, ntime+1);
    EKinet = similar(EElect);

    Ek_out = zeros(Complex{Float64}, ntime+1 ,Nx);

    ## Initial Electric Field and Energies    
    rhok, Ek = poisson1d(k, fk[1:end,1]);
    
    EElect[1] = EEfield(L, Ek);
    EKinet[1] = EKin(fku, du);
    
    ## Electric field output
    Ek_out[1,:] = Ek;
    E = realspace1d(Ek, ip1E);

    # For poisson1d only. Invertes rhok and Ek
    # exhaustive in-place planning pays off in time for several uses
    ip1E = plan_ifft!(E, flags = FFTW.EXHAUSTIVE);
    
    ## Chin A integrator coefficients
    coef1 = dt/2;
    coef2 = 2*dt/3;
    coef3 = dt/6;
    coef4 = dt*dt/24;
    
    Xshift = exp.(-1im*coef1*k*v'); # X-space advance according to A integrator
    Vshift = exp.(1im*coef3*E*u');  # V-space advance according to A integrator

    ## Filter for de-aliasing    
    umax = maximum(u);                         # Maximum velocity wavevector
    filterv = exp.(-36.*(abs.(u')/umax).^36);  # Ta-dá, the filter

    
#########  MAIN LOOP STARTS HERE  ########

    #Indicate progress every itdiv iterations in file "progress"
    strToFile("progress","started $(now()) \n","w")
    itdiv = (ntime+1) / instToSave
    stime = time()
    inst = 1

    for nt = 2:ntime+1
### First v-advection step
        fku = fku.*Vshift;
        
### First x-advection step
        fku = fku.*Xshift;
        
### Middle v-advection step
        ## Electric field 
        rhok, Ek = poisson1d(k, fku[1:end,1]);
        dEdx = realspace1d(rhok,ip1E);
        E = realspace1d(Ek,ip1E);

        G = coef4*E.*dEdx;
        Vshift = exp.(1im*coef2*(E + G)*u');
        fku = fku.*Vshift;
        fku = fku.*filterv;
        
### Last x-advection step
        fku = fku.*Xshift;
        
### Last v-advection step plus the first one
        ## Electric field 
        rhok, Ek = poisson1d(k, fk[1:end,1]);
        E = realspace1d(E,ip1E);
        
        ## The v-shift
        Vshift = exp.(1im*coef3*E*u');
        fku = fku.*Vshift;
        
## Electric field & energies. Final iteration value
        Ekin[nt] = EKinetic(fku, du);  # Kinetic energy from distribution output
        EElect[nt] = EEfield(L, Ek);      # Electric field energy output
        Ek_out[nt,:] = Ek;                # Electric field output

## Indicates progress every 1/8 of the total iterations
        if (mod(nt, itdiv) <1 )
            inst = inst + 1
            f_out[inst,:,:] = fku
            strToFile("progress","$(floor(100.0*nt/ntime)) % accomplished in $(time()-stime) secs.\n","a")
        end

    end
    # return quantities in the fourier space
    return f_out, Ek_out, EElect, Ekin;
    
end
