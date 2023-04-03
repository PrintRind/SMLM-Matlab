function [I_molecule, del_idx] = crop_molecule_images(image_stack, preloc_positions, frame_no, PSF)
%crops single molecule images from the raw data stack
%outputs: stack of cropped single molecule images and a vector of indices
%that should be deleted (molecules too close too boundaries)

    Nx = PSF.N_image - 2; %side length of a molecule image; make it a bit smaller than the model to prevent fitting errors
    N_img = size(preloc_positions,1);
    [N_rows, N_cols, N_frames] = size(image_stack);
    I_molecule = zeros(Nx,Nx,N_img); 
    v = (1:Nx) - mean(1:Nx); %cropping index range
    
    x = preloc_positions(:,1)/PSF.ux/1e9;  %column center position in pixel
    y = preloc_positions(:,2)/PSF.ux/1e9;  %row position in pixel
    
    idx = 0; 
    del_idx = []; 
    for m = 1:N_frames %loop over frames
        image = image_stack(:,:,m); %select frame
        
        for n = 1:sum(frame_no == m) %loop over molecules in this frame
            idx = idx + 1; 
            
            row_range =  round(y(idx) + v); 
            col_range =  round(x(idx) + v); 
            
            if 1<=min(row_range) && N_rows>=max(row_range) && 1<=min(col_range) && N_cols>=max(col_range)
                I_molecule(:,:,idx) = image(row_range, col_range); 
                
            else
                I_molecule(:,:,idx) = 0; 
                del_idx = [del_idx; idx]; 
            end
            %imagesc(I_molecule(:,:,idx)); pause(0.1); 
        end
    end

end

