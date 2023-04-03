function pos = cell2array(positions)

    pos = []; 
    for m = 1:numel(positions)
        pos = cat(1, pos, positions{m}); 
    end
end
