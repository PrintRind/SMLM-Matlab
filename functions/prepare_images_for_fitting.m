%prepare experimental images for fitting
function Image = prepare_images_for_fitting(I_cts, camera)
    Image = double(I_cts - camera.baseline); %subtract baseline
    Image = Image * camera.amp/camera.QE; %go from counts to photo-electrons; correct for the QE --> photons
end