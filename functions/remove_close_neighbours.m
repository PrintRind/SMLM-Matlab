function [X, Y] = remove_close_neighbours(X, Y, dist_thresh)
        belowTh = [];
        tmpPos = [X,Y];
        % calculate distances between points in upper triangle
        % matrix
        for ix=1:size(tmpPos,1)
            
            for iy = ix+1:size(tmpPos,1)
                if norm(tmpPos(ix,:)-tmpPos(iy,:)) < dist_thresh %are molecules too close?
                    belowTh = [belowTh; [ix,iy]]; %belowTh contains row indices of molecule pairs that are too close
                end
            end
        end
        % delete points below threshold
        if not(isempty(belowTh))
            X(unique([belowTh(:,1);belowTh(:,2)])) = [];
            Y(unique([belowTh(:,1);belowTh(:,2)])) = [];
        end
end