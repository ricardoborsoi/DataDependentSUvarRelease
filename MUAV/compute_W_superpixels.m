function [imout] = compute_W_superpixels(imin, nr, nc, L, spSegs, numSuperpixels, mode)
% -------------------------------------------------------------------------
% Compute the result of the operator W
% 
% 
% 
% 
% 
% 
% -------------------------------------------------------------------------





% *************************************************************************
if strcmp(mode, 'W')
    
    Y_cube = imin;
%     Y3 = zeros(size(Y_cube));
    avg_superpx = zeros(1, numSuperpixels+1, L);
    
    for i=0:numSuperpixels
        [rowi, coli] = find(spSegs==i);

        for j=1:length(rowi)
            % Averages all pixels inside each superpixel
            if j == 1
                avg_superpx(1,i+1,:) = (1/length(rowi)) * Y_cube(rowi(j),coli(j),:);
            else
                avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y_cube(rowi(j),coli(j),:);
            end
        end
    
    %     % This is optional (for visualization)
    %     for j=1:length(rowi)
    %         Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
    %     end
    end
    
    avg_superpx = squeeze(avg_superpx)';
    imout = avg_superpx;
    
    
    % *********************************************************************
elseif strcmp(mode, 'Wast')
    Y3 = zeros(nr, nc, L);
%     avg_superpx = zeros(1, numSuperpixels+1, L);
%     avg_superpx = ==imin;
    
%     avg_superpx = permute(imin, [3 1 2]);
    avg_superpx = permute(imin, [3 2 1]);
    
%     size(imin)
%     size(Y3)
%     size(avg_superpx)
    
    for i=0:numSuperpixels
        [rowi, coli] = find(spSegs==i);
        
%         for j=1:length(rowi)
%             % Averages all pixels inside each superpixel
%             if j == 1
%                 avg_superpx(1,i+1,:) = (1/length(rowi)) * Y_cube(rowi(j),coli(j),:);
%             else
%                 avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y_cube(rowi(j),coli(j),:);
%             end
%         end

        % This is optional (for visualization)
        for j=1:length(rowi)
            Y3(rowi(j),coli(j),:) = avg_superpx(1, i+1, :);
        end
    end
    
    imout = Y3;

    
    
    
    % *********************************************************************
elseif strcmp(mode, 'WWast')
    
    
    Y_cube = imin;
    Y3 = zeros(size(Y_cube));
    avg_superpx = zeros(1, numSuperpixels+1, L);
    
    for i=0:numSuperpixels
        [rowi, coli] = find(spSegs==i);

        for j=1:length(rowi)
            % Averages all pixels inside each superpixel
            if j == 1
                avg_superpx(1,i+1,:) = (1/length(rowi)) * Y_cube(rowi(j),coli(j),:);
            else
                avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y_cube(rowi(j),coli(j),:);
            end
        end

        % This is optional (for visualization)
        for j=1:length(rowi)
            Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
        end
    end
    
%     avg_superpx = squeeze(avg_superpx)';
    imout = Y3;
    
    
    
    % *********************************************************************
else
    error('Unknown mode specified!!!')
end



