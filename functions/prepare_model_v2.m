function [f, df] = prepare_model_v2(F, dFx, dFy, dFz, x, y, ux)           
    %F, dFx..etc....spline objects of the 3D PSF and its spatial
    %derivatives
    %x,y,z..vectors defining the coordinates
    %v...parameter vector [x,y,z,sig,bg]
    f = @model;
    df = @d_model;

    function f = model(v)
        %this function is the one minimized by the levenberg-marquardt
        %algorithm (see LM_poisson.m)
        f = v(5) + v(4) * fnval(F, {x-v(1),y-v(2), 1 + v(3)}); 
    end

    function J = d_model(f_a, v)
        df_x = -v(4)*fnval(dFx, {x-v(1),y-v(2),1 + v(3)});
        df_y = -v(4)*fnval(dFy, {x-v(1),y-v(2),1 + v(3)});
        df_z = v(4)*fnval(dFz, {x-v(1),y-v(2),1 + v(3)});
        df_signal = (f_a - v(5))/v(4); 
        
        J = [df_x(:), df_y(:), df_z(:), df_signal(:), ones(numel(f_a),1)]; 
    end
        
end
