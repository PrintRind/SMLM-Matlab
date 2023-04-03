%calculates the 3D PSF of a TIRF microscope

function PSF = calculate_PSF(Nk, lambda_0, NA, RI_fluid, RI_imm, z_max, dz, z_defocus, Nx, ux, Z_modes, Z_phase, Z_amp)

    %-------test parameters
%     Nk = 128;       %pupil diameter
%     Nx = 17;        
%     RI_fluid = 1.33; 
%     RI_imm = 1.518; 
%     lambda_0 = 670e-9;
%     Z_modes = [5, 6, 7, 8];
%     Z_phase = zeros(1,length(Z_modes));
%     Z_amp = [1, 0, 0, 0]; 
%     NA = 1.49;
%     ux = 100e-9; 
%     z_defocus = -400e-9;
%     dz = 10e-9;     %z-increment of PSF calculation
%     z_max = 200e-9; %max. z-value (should be within TIRF range)
    %----------------------------

    d_layer = 0; %intermed. layer tickness
    RI_layer = RI_fluid; %intermediate layer RI
    RI = [RI_fluid, RI_layer, RI_imm]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]
    k0 = 2*pi/lambda_0; 
    z_vec = 0:dz:z_max;   
    Nz = length(z_vec); 

    %dipole magnitudes
    mu_x = 1e-9; 
    mu_y = 1e-9; 
    mu_z = 1e-9; 
    f_obj = 1e-3; 

    [aberr, T] = model_pupil(Nk, Z_modes, Z_phase, Z_amp); %modelling pupil phase and transmission function from Z_phase/Z_amp
    
    %calculating BFP-fields for all dipole orientations
    PSF = zeros(Nx,Nx,Nz);
    E_tot = zeros(Nz); 

    %calculating defocus function
    uk = (2*k0*NA)/Nk; 
    [~,~,Kr,~] = create_coord(Nk,uk,'exact');
    Kz = sqrt(k0^2*RI_imm^2 - Kr.^2); 
    mask = T .* exp(1i*aberr + 1i*Kz*z_defocus);
    pupil_full = Kr <= k0*NA; 
    pupil_UAF = Kr <= k0*RI_fluid; 

    for m=1:length(z_vec)
        %dipole=[0,0]; %[theta, phi], e.g. [0,0] for z-dipole, [pi/2,0] for x-dipole
        [Ex_Pz, Ey_Pz] = fun_dipole_imaging(Nk, lambda_0 , NA, RI,[0,0],d_layer,z_vec(m),f_obj,mu_z,T); %z-dipole
        [Ex_Px, Ey_Px] = fun_dipole_imaging(Nk ,lambda_0 , NA, RI,[pi/2, 0],d_layer,z_vec(m),f_obj,mu_x,T); %x-dipole
        [Ex_Py, Ey_Py] = fun_dipole_imaging(Nk, lambda_0 , NA, RI,[pi/2, pi/2],d_layer,z_vec(m),f_obj,mu_y,T); %y-dipole 
    
        I_xx=abs(czt2(Ex_Px .* mask, uk, ux, Nx)).^2;
        I_yx=abs(czt2(Ey_Px .* mask, uk, ux, Nx)).^2;
        I_xy=abs(czt2(Ex_Py .* mask, uk, ux, Nx)).^2;
        I_yy=abs(czt2(Ey_Py .* mask, uk, ux, Nx)).^2;    
        I_xz=abs(czt2(Ex_Pz .* mask, uk, ux, Nx)).^2;
        I_yz=abs(czt2(Ey_Pz .* mask, uk, ux, Nx)).^2;   

        E_tot(m) = sum(sum(abs(Ex_Px).^2)) + sum(sum(abs(Ex_Py).^2)) + sum(sum(abs(Ex_Pz).^2)) + sum(sum(abs(Ey_Px).^2)) + sum(sum(abs(Ey_Py).^2)) + sum(sum(abs(Ey_Pz).^2));
        PSF(:,:,m)=(I_xx+I_yx+I_xy+I_yy+I_xz+I_yz)/E_tot(m); %normalization to total energy in BFP; we thus assume that the SAME energy comes through the objective lens for every z-distance (note that the energy would otherwise be different in the simulation for varying z-distances)
    end




