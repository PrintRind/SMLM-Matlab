function J_Alex = Jacobian_Alex(f, v0)
   
    if size(v0,1) == 1; v0 = transpose(v0); end

    dx = 1e-2; 
    dy = dx; 
    dz = dx; 
    ds = 1; 

    dxv = [dx, 0, 0, 0, 0]'; 
    dyv = [0, dy, 0, 0, 0]'; 
    dzv = [0, 0, dz, 0, 0]'; 
    dsv = [0, 0, 0, ds, 0]'; 

    Jx = (f(v0 + dxv) - f(v0 - dxv))/(2*dx); 
    Jy = (f(v0 + dyv) - f(v0 - dyv))/(2*dy); 
    Jz = (f(v0 + dzv) - f(v0 - dzv))/(2*dz); 
    Js = (f(v0 + dsv) - f(v0 - dsv))/(2*ds); 
    %Js = (f(v0) - v0(5))/v0(4); 
  
    J_Alex = [Jx(:), Jy(:), Jz(:), Js(:), ones(numel(Jx),1) ]; 
    
end