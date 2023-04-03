function [Phi, T] = model_pupil(dia_pupil, Z_modes, Z_phase, Z_amp)

    [~,~,R,pupil]=create_coord(dia_pupil, 2/dia_pupil ,'exact');
    
    %modelling aberrations and lens transmision
    Phi = sum(ZernikeCalc(Z_modes, Z_phase',pupil,'NOLL'),3); %aberration file "coef.mat" must be loaded 
   
    %additional apodization in the objective pupil      
    T = Z_amp(1).*pupil - Z_amp(2).*R.^2 - Z_amp(3).*R.^4 - Z_amp(4).*R.^6; 
    T = T/max(T(:)); 

