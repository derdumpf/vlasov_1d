# FIELD RELATED FUNCTIONS

function poisson1d( kx::Array{Float64}, fku_0::Array{Complex{Float64}} );
    #   NOTE: _ fk = fk(kx,kv)
    #         _ Ek = Ek(kx)
    #         _ rhok = rhok(kx)
        
    rhok = fku_0;             # rho_k
    rhok[1] = 0.0;            # remove (ion) background
  
    Ek = rhok./(-1im*kx);      
    Ek[1] = 0.0;              # no field due to mean density

    return rhok, Ek;
end


function EEfield( L::Float64, Ek::Array{Complex{Float64}} ); 
    #   NOTE: _ Ek = Ek(kx)
    #         _ E_field_n (just for the present instant)
    
    E_field_n = 2pi*sum(Ek.*Ek)/(2.*L);
    
    return E_field_n;
end

function EKinetic(fku::Array{Complex{Float64}}, du::Float64);
    #  NOTE: _fku = fku(k,u)
    # 
    
    fku_row = fku[1,1:end];
    Kin = ( fku_row[3] - 16*fku_row[2] + 30*fku_row[1] - 16*fku_row[end] + fku_row[end-1] )/(12*du*du)

    return Kin;
end
