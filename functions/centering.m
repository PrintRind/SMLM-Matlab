
function [stack_cent, m_focus] = centering(stack) 
    [Nx, Ny, Nz ] = size(stack)
    [~, maxidx] = max(stack(:))
    [maxrow, maxcol, m_focus] = ind2sub(size(stack), maxidx)
 
    %centering the PSF within the frame
    stack_cent = circshift(stack, round(([Nx Ny]+1)/2)-[maxrow, maxcol]); 