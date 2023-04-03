%estimate initial parameters for precise localization
function v0 = estimate_v0(I_photons)

    [N_img, Nx, Ny] = size(I_photons);
    
    frame = zeros(Nx,Ny); 
    frame(:,1) = 1; 
    frame(:,end) = 1; 
    frame(1,:) = 1; 
    frame(end,:) = 1; 
    
    v0 = zeros(N_img, 5);
    
    for m = 1:N_img
        tmp = squeeze(I_photons(m,:,:));
        bg = mean(tmp(frame==1)); 
        signal = sum(tmp(:) - bg);
        v0(m,:) = [0, 0 , 10, signal, bg];
    end


end
