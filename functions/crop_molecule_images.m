function [I_molecule, TS] = crop_molecule_images(image_stack, TS, PSF, fw)
%crops single molecule images from the raw data stack
%outputs: stack of cropped single molecule images and a vector of indices
%that should be deleted (molecules too close too boundaries)
%fw...border width
%TS...thunderstorm table
%d_min..minimum distance between two neighbouring molecules

 
    frame_no = table2array(TS(:, 'frame')); %vector indicating the frame number containing the molecule

    Nx = PSF.N_image - 2*fw; %side length of a molecule image; make it a bit smaller than the model to prevent fitting errors
    N_img = size(frame_no,1);
    [N_rows, N_cols, N_frames] = size(image_stack);

    %v = (1:Nx) - mean(1:Nx); %cropping index range
    v = (1:Nx) - ceil(Nx/2);

    x = table2array(TS(:, 'x_nm_'))/PSF.ux/1e9;  %column center position in pixels
    y = table2array(TS(:, 'y_nm_'))/PSF.ux/1e9;  %row position in pixels

    N_img = size(x,1); %number of images left
    I_molecule = zeros(Nx,Nx,N_img); %init. stack of molecules images

    idx = 0; 
    del_idx = []; 
    for m = 1:N_frames %loop over frames
        image = image_stack(:,:,m); %select frame
        
        for n = 1:sum(frame_no == m) %loop over molecules in this frame
            idx = idx + 1; 
            
            row_range =  ceil(y(idx) + v); 
            col_range =  ceil(x(idx) + v); 
            
            if 1<=min(row_range) && N_rows>=max(row_range) && 1<=min(col_range) && N_cols>=max(col_range)
                I_molecule(:,:,idx) = image(row_range, col_range); 
                
            else %delete molecule if it is too close to the image borders
                I_molecule(:,:,idx) = 0; 
                del_idx = [del_idx; idx]; 
            end
            %imagesc(I_molecule(:,:,idx)); pause(0.1); 
        end
    end
    
    %deleting faulty entries
    I_molecule(:,:,del_idx) = [];
    TS(del_idx,:) = []; 
    disp([num2str(length(del_idx)), ' detected molecules deleted because of too high proximity to the camera image borders.']);

end

