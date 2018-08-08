# All functions related to Fourier-calculations
# -> This includes transforms and creation of Fourier-space-quantities

function wavevector(var::Array{Float64})
    var_range = var[end]-var[1];
    N_var = length(var);

    k_nod = collect(1:N_var)  - N_var/2 -1;
    return 2pi*fftshift(k_nod)/var_range;    
end


function realspace1d(A::Array{Complex{Float64}}, ipf)
    return real(ipf*A)
end


function realspace1d(A::Array{Complex{Float64}}) # No Plan overloading
    return real( ifft(A) )
end


function fourierspace2d(A::Array{Float64}, p1, p2)
    # A: a two dimensional real quantity
    # p1, p2: FFT Plans for A in dimensions 1 and 2
    # Returns the real-space form of a 2d quantity taking its iFFT
    #    Plans as argument
    return p2*(p1*A);
end


function fourierspace2d(A::Array{Float64}) # No Plan Overloading
    # A: a two dimensional real quantity
    # Returns the real-space form of a 2d quantity taking its iFFT
    #    Plans as argument
    return fft(fft(A, 1), 2);
end


function realspace2d(A::Array{Complex{Float64}}, ip1, ip2);
    # A: a two dimmensional Fourier quantity
    # ip1, ip2: iFFT Plans for A in dimensions 1 and 2
    # Returns the real-space form of a 2d quantity taking its iFFT
    #    Plans as argument
    A = ip2*(ip1*A)
    return real(A)
end


function realspace2d(A::Array{Complex{Float64}}); # No Plan overloading
    # A: a two dimmensional Fourier quantity
    # Returns the real-space form of a 2d quantity without plans
     A = ifft( ifft(A, 1), 2 )
    return real(A)
end
