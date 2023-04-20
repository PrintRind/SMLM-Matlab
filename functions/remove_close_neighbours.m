function [TS_table_filt, images_filt] = remove_close_neighbours(TS_table, images, dist_thresh)
        %removing molecules that are too closely spaced 
        %TS_table...thunderstorm table ["ID", "frame", "x/nm", "y/nm"]
        %images...stack of single molecule images (already cropped) size e.g.[17, 17, 3000]
        %x...column coordinates in pixels
        %y...row coordinates in pixels
        %frame_no...vector containing the frame number where the molecule
        %appears
        %dist_thresh...minimum distance in pixels; molecules more closely
        %spaced are both removed
        
        TS = table2array(TS_table);
        frame_no = TS(:,2);
        frames_list = unique(frame_no); 
        no_frames = length(frames_list); %number of frames in the stack     
     
        TS_filt = []; 
        images_filt = []; 
        % calculate distances between points in upper triangle
        % matrix

        for m = 1:no_frames 

            belowTh = [];
            del_idx = [];
            idx = find(frame_no == frames_list(m)); %vector of indices of molecules in frame m
           
            N_mol = length(idx); %no. of molecules in this frame

            tmpPos = [TS(idx,3), TS(idx,4)]; %getting coords. of molecules in frame m
            tmpframes = frame_no(idx); 

            for iA=1:N_mol
                
                for iB = iA+1:N_mol
                    if norm(tmpPos(iA,:)-tmpPos(iB,:)) < dist_thresh %are molecules too close?
                        belowTh = [belowTh; [iA,iB]]; %belowTh contains row indices of molecule pairs that are too close
                    end
                end
            end

            % delete points below threshold
            if not(isempty(belowTh))
                del_idx = [del_idx; unique([belowTh(:,1);belowTh(:,2)])]; %row index of molecules to be deleted            
            end
            
            %deleting entries in TS_filt and images_filt
            tmp = TS(idx,:); 
            tmp(del_idx,:) = [];
            TS_filt = [TS_filt; tmp]; 

            tmp = images(:,:,idx);
            tmp(:,:,del_idx) = []; 
            images_filt = cat(3, images_filt, tmp); 
        end
        
        TS_table_filt = array2table(TS_filt); 
        TS_table_filt.Properties.VariableNames = TS_table.Properties.VariableNames; 
        no_deleted = size(TS_table,1) - size(TS_table_filt,1); 
        fraction_deleted = round(100*no_deleted/size(TS_table,1));
        disp([num2str(fraction_deleted), '% of entries (', num2str(no_deleted), ' molecules) deleted.']); 
end