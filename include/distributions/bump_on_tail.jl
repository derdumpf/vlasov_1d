########################################################################################
######################     BUMP-ON-TAIL PARTICLE DISTRIBUTION    #######################
########################################################################################

function bump_on_tail(v::Float64[], rhoc::Float64, rhob::Float64, vtc::Float64, vtb::Float64, V0c::Float64,V0b::Float64)

    f0 = Float64[];   
    fnc = rhoc/(sqrt(2pi)*vtc); # Normalization factor for core species
    fnb = rhob/(sqrt(pi)*vtb);  # Normalization factor for beam species
 
    f0 = fnc.*exp.(-((v - V0c)./vtc).^2./2.)' + fnb.*exp.(-((v - V0b)./vtb).^2./2.)';
    
    return f0;
end
