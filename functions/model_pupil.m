function [Phi, T] = model_pupil(dia_pupil, Z_modes, Z_mag, T_coefs)

    [~,~,R,pupil]=create_coord(dia_pupil, 2/dia_pupil ,'exact');
    
    %modelling aberrations and lens transmision
    Phi = sum(ZernikeCalc(Z_modes, Z_mag',pupil,'NOLL'),3); %aberration file "coef.mat" must be loaded 
   
    %additional apodization in the objective pupil      
    T = T_coefs(1).*pupil - T_coefs(2).*R.^2 - T_coefs(3).*R.^4 - T_coefs(4).*R.^6; 
    T = T /max(T(:)); 
    T = T .* pupil; 
end
