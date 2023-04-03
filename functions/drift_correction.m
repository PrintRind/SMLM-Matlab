function preloc_positions = drift_correction(preloc_positions, table, ux)
    %preloc_positions contains the frame number (1st column) and the x,y
    %positions (2nd, 3rd columns)
    %this function corrects the position vectors x and y, which must be given in units of
    %pixels. csv_filename is the drift table exported in Thunderstorm if you click on "save" in the drift-popup window.  
    %ux is the pixel unit in meters
    %table = readtable(csv_filename); %T is a csv-table exported from Thundertorm (Drift correction)
    x0 = table.X0;
    y0 = table.Y0;
    x1 = table.X1;
    y1 = table.Y1;

    x0(isnan(x0)) = [];
    y0(isnan(y0)) = [];
    x1(isnan(x1)) = [];
    y1(isnan(y1)) = [];
    t = linspace(min(x0), max(x0), 100); %time axis
    
    max_deg = min(length(x0)); %degree of polyfit

    p_x = polyfit(x0,y0,max_deg-1);
    p_y = polyfit(x1,y1,max_deg-1);

    x_fit = polyval(p_x, t);
    y_fit = polyval(p_y, t);

    figure(1)
    plot(t,x_fit, t, y_fit);
    hold on; 
    plot(x0,y0, 'o');
    plot(x1,y1, 'o');
    hold off;
    grid on; 
    xlabel("frame no."); 
    ylabel("drift / pixel"); 

    t_frame = preloc_positions(:,1); 
    preloc_positions(:,2) = preloc_positions(:,2) - polyval(p_x, t_frame)*ux*1e9; 
    preloc_positions(:,3) = preloc_positions(:,3) - polyval(p_y, t_frame)*ux*1e9; 

end