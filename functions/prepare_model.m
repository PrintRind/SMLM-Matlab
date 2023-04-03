function f = prepare_model(F, X, Y, Z)           
    %F...Interpolant object of the 3D PSF
    %X,Y,Z...2D matrices defining the coordinates
    %v...parameter vector [x,y,z,sig,bg]
    f = @model;
    function f = model(v)
        %this function is the one minimized by the levenberg-marquardt
        %algorithm (see LM_poisson.m)
        f = v(5) + v(4) * F(Y - v(2), X - v(1), Z + v(3));  %define function which calculates theoretical molecule image given the parameters in vector v
    end        
end
