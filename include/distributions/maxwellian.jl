########################################################################################
######################      MAXWELLIAN PARTICLE DISTRIBUTION     #######################
########################################################################################

function maxwellian(v::Float64[], rho::Float64, vt::Float64, V0::Float64)  

    f0 = Float64[];   
    fn = rho/(sqrt(2pi)*vt); 
    f0 = fnc.*exp.(-((v - V0c)./vtc).^2./2.)';

    return f0;
end