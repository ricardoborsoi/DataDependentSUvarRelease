function [mtx] = ilexico(data, nr, nc, opt)
% =========================================================================
% Convert from lexicographic ordering mack to matrix
% 
% INPUT:
%    data - lexicographically ordered input
%    nr   - number of rows of the desired image
%    nc   - number of columns of the desired image
%    opt  - way that the lexicographic ordering is performed, that is, 
%           row-wise or column wise. Options are 'row' and 'col'
% 
% OUTPUT:
%    mtx - image ordered as a cube
% 
% -------------------------------------------------------------------------
% Ricardo Borsoi
% 2017 10 10
% =========================================================================

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if length(size(data)) == 2
    
    [sz1, sz2] = size(data);
    
    if ~(sz1 == nr*nc || sz2 == nr*nc)
        error('Inconsistent image sizes!')
    end
    P = min(sz1,sz2);
    
    if sz1 > sz2
        data = data';
    end
    
    switch opt
        case 'row'
            % Convert row-wise lexocographic order back to cube
            mtx = permute(reshape(data', nc, nr, P),[2 1 3]);
        case 'col'
            % Convert column-wise lexocographic order back to cube
            mtx = reshape(data', nr, nc, P);
    end
    
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
elseif length(size(data)) == 3
    
    [sz1, sz2, sz3] = size(data);

    if ~(sz3 == nr*nc)
        error('Inconsistent image sizes!')
    end
    L = sz1;
    P = sz2;
    
    mtx = zeros(L,P,nr,nc);
    
    switch opt
        case 'row' % --->
            % Convert row-wise lexocographic order back to cube
            n=0;
            for i=1:nr
                for j=1:nc
                    n = n+1;
                    mtx(:,:,i,j) = data(:,:,n);
                end
            end
            
        case 'col' % \/ \/ \/
            % Convert column-wise lexocographic order back to cube
            n=0;
            for j=1:nc
                for i=1:nr
                    n = n+1;
                    mtx(:,:,i,j) = data(:,:,n);
                end
            end
    end
    
    
    % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
else
    error('Data have too many of too few dimensions!')
    
end


