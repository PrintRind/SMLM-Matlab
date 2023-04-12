%estimate initial parameters for precise localization
function v0 = estimate_v0(I_photons)

    [Ny, Nx, N_img] = size(I_photons);
    
    frame = zeros(Ny,Nx); 
    frame(:,1) = 1; 
    frame(:,end) = 1; 
    frame(1,:) = 1; 
    frame(end,:) = 1; 
    
    v0 = zeros(N_img, 5);
    
    
    for m = 1:N_img
        img = I_photons(:,:,m);
        tmp = squeeze(img);
        bg = mean(tmp(frame==1)); 
        signal = sum(tmp(:) - bg);

        [~, idx] = max(img(:));
        [y_max, x_max] = ind2sub([Ny,Nx], idx);

        v0(m,:) = [x_max - round(Nx/2), y_max - round(Ny/2) , 10, signal, bg];
    end


end
