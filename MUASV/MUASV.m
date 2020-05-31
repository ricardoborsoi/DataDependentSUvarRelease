function [A,psi_maps,M] = MUASV(data, A_init, psis_init, M0,lambda_m,lambda_a,lambda_psi,slic_size,slic_reg,rho_a,varargin)
% 
%ELMM_ADMM Unmix hyperspectral data using the Extended Linear Mixing Model
%(ELMM).
%
%   We find a stationary point of the following functional:
%   J(S,A,PSI) = 1/2 * sum_{k=1}^{N} (||x_k - S_k*a_k||_{2}^{2} + 
%   ||S_k - S0*psi_{k}||_{F}^{2}) + lambda_A R(A) + lambda_PSI R(PSI)
%
%   with S a collection of endmember matrices for each pixel, A the
%   abundances in each pixel and for each endmember, and PSI the scaling
%   factors advocated by the ELMM.
%
%   The abundances are subject to the usual nonnegativity and sum to one
%   constraints. The scaling factors and endmember spectral are nonnegative
%   as well.
%
%   R(A) is a spatial regularization term on the abundances. It can be 
%   either an anisotropic total variation term TV(A) applied on each 
%   material or a Tikhonov like regularization on the spatial gradients of 
%   the abundance maps. R(PSI) is a differentiable Tikhonov regularization 
%   on the spatial gradients of the scaling factor maps.
%
%   Mandatory inputs:
%   -data: m*n*L image cube, where m is the number of rows, n the number of
%   columns, and L the number of spectral bands.
%   -A_init: P*N initial abundance matrix, with P the number of endmembers
%   to consider, and N the number of pixels (N=m*n)
%   -psis_init: P*N initial scaling factor matrix
%   -S0: L*P reference endmember matrix
%   -lambda_m: regularization parameter on the ELMM tightness
%   -lambda_a: regularization parameter for the spatial regularization on
%   the abundances.
%   -lambda_psi: regularization parameter for the spatial regularization on
%   the scaling factors
%   The spatial regularization parameters can be scalars, in which case 
%   they will apply in the same way for all the terms of the concerned
%   regularizations. If they are vectors, then each term of the sum
%   (corresponding to each material) will be differently weighted by the
%   different entries of the vector.
%   - slic_size: 
%   - slic_reg: 
%   - rho_a: 
%
%   Optional inputs (arguments are to be provided in the same order as in 
%   the following list):
%   -verbose: flag for display in console. Display if true, no display
%   otherwise (default: true)
%   -maxiter_anls: maximum number of iterations for the ANLS loop (default:
%   100)
%   -epsilon_s: tolerance on the relative variation of S between two
%   consecutive iterations (default: 10^(-3))
%   -epsilon_a: tolerance on the relative variation of A between two
%   consecutive iterations (default: 10^(-3))
%   -epsilon_psi: tolerance on the relative variation of psi between two
%   consecutive iterations (default: 10^(-3))
%
%   Outputs:
%   -A: P*N abundance matrix
%   -psi_maps: P*N scaling factor matrix
%   -S: L*P*N tensor constaining all the endmember matrices for each pixel
%   -optim_struct: structure containing the values of the objective
%   function and its different terms at each iteration
%
%   An example of the use of this function on a real dataset is shown in 
%   the file demo.m.
%
%   The algorithm is presented in detail in:
%
%   L. Drumetz, M. A. Veganzones, S. Henrot, R. Phlypo, J. Chanussot and 
%   C. Jutten, "Blind Hyperspectral Unmixing Using an Extended Linear
%   Mixing Model to Address Spectral Variability," in IEEE Transactions on 
%   Image Processing, vol. 25, no. 8, pp. 3890-3905, Aug. 2016.
%
%   Last Revision: 17-November-2016.
%   Revision: 1.0
%


% ========================================================
% Get input parameters

minnargin = 10; % check number of inputs
maxnargin = 16;
narginchk(minnargin,maxnargin);

P = size(A_init,1); % number of endmembers

% find out if regularization parameters are scalar or vectors

% abundances
if length(lambda_a) == 1
    scalar_lambda_a = 1;
elseif length(lambda_a) == P
    scalar_lambda_a = 0;
    if size(lambda_a,1) == 1
        lambda_a = lambda_a';
    end
else
    error('ELMM_ADMM:wrong_size_regularization','lambda_a must be a scalar or a P-dimensional vector')
end

% scaling factors
if length(lambda_psi) == 1
    scalar_lambda_psi = 1;
elseif length(lambda_psi) == P
    scalar_lambda_psi = 0;
    if size(lambda_psi,1) == 1
        lambda_psi = lambda_psi';
    end
else
    error('ELMM_ADMM:wrong_size_regularization','lambda_psi must be a scalar or a P-dimensional vector')
end

% initialize optional parameters

maxiter_anls = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
norm_sr = '1,1';
verbose = true;
flag_use_parfor = false;

% change optional parameters if provided
for i = 1:length(varargin)
    switch i
        case 1
            verbose = varargin{i};
        case 2
            maxiter_anls = varargin{i};
        case 3
            epsilon_s = varargin{i};
        case 4
            epsilon_a = varargin{i};
        case 5
            epsilon_psi = varargin{i};
        case 6
            flag_use_parfor = varargin{i};
    end
end



% ========================================================
%% initialize variables

[m,n,L] = size(data); % dimensions of the data cube
N = m*n; % number of pixels

data_r = reshape(data,N,L)'; % store image as a matrix

% relative variations of the variables
rm = zeros(maxiter_anls,1);
ra = zeros(maxiter_anls,1);
rpsi =  zeros(maxiter_anls,1);

% initialize variables
A = A_init; % initialize the abundance matrix
M = repmat(M0,[1,1,N]); % initialize pixel-dependent endmember matrix
psi_maps = psis_init; % initialize scaling factors

M0ptM0 = diag(M0'*M0); % precomputing



% build handlers and necessary stuff

% forward first order horizontal difference operator
FDh = zeros(m,n);
FDh(end,1) = -1;
FDh(end,end) = 1;
FDh = fft2(FDh);

% forward first order vertical difference operator
FDv = zeros(m,n);
FDv(1,end) = -1;
FDv(end,end) = 1;
FDv = fft2(FDv);



% auxiliary functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define a circular convolution (the same for all bands) accepting a
% matrix  and returning a matrix

ConvC = @(X,FK)  reshape(real(ifft2(fft2(reshape(X', m,n,P)).*repmat(FK,[1,1,P]))), m*n,P)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert matrix to image
conv2im  = @(A)  reshape(A',m,n,P);
% convert image to matrix
conv2mat = @(A)  reshape(A,m*n,P)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Compute superpixels decomposition
Y2a = data;

% reorder and rescale data into 2-D array
[numRows,numCols,numSpectra] = size(Y2a);
scfact = mean(reshape(sqrt(sum(Y2a.^2,3)), numRows*numCols, 1));
Y2a = Y2a./scfact;
imgVec = reshape(Y2a, [numRows*numCols numSpectra]);

% compute superpixels
disp('Computing SLIC Superpixels...');
spSegs = vl_slic(single(Y2a), slic_size, slic_reg);
numSuperpixels = double(max(spSegs(:)))+1; 

%[nr, nc, L] = size(data);
[m, n, L] = size(data);

% Decompose Y into approximation and detail scales
Yc = compute_W_superpixels(data, m, n, L, spSegs, numSuperpixels, 'W');
Yd = data - compute_W_superpixels(data, m, n, L, spSegs, numSuperpixels, 'WWast');


% Initialize variables
Yc_mtx = Yc;
Yd_mtx = reshape(Yd, N, L)';

Ac = zeros(P, numSuperpixels);
Ad = zeros(P, N);

Ac_old = zeros(P, numSuperpixels);
Ad_old = zeros(P, N);

AcWast_old = zeros(P, N);
Ad_old = zeros(P, N);

Mc = zeros(L,P,numSuperpixels);
Md = zeros(L,P,N);



% Compute sum-to-one constraints
sumW = compute_W_superpixels(ones(m,n,1), m, n, 1, spSegs, numSuperpixels, 'W')';
sumImWWast = ones(1,N) - reshape(compute_W_superpixels(ones(m,n), m, n, 1, spSegs, numSuperpixels, 'WWast'), m*n, 1)';


% Initialize variables related to the nonnegative least squares solver
delta_asc = 1000;
P_old_fnnls_c = repmat(zeros(1,P), [numSuperpixels 1]);
Z_old_fnnls_c = repmat(1:P, [numSuperpixels 1]);

P_old_fnnls_d = repmat(zeros(1,P), [N 1]);
Z_old_fnnls_d = repmat(1:P, [N 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% optimization

for i = 1:maxiter_anls
    
    % store previous values of the variables
    M_old = M;
    psi_maps_old = psi_maps;
    A_old_anls = A;
    
    %% M_update
    
    if verbose
        fprintf('updating M...\n')
    end
        
    if flag_use_parfor
        parfor k =1:N
            M(:,:,k) = (data_r(:,k)*A(:,k)' + lambda_m*M0*diag(psi_maps(:,k))) ...   %
                /(A(:,k)*A(:,k)' + lambda_m*eye(P));
            M(:,:,k) = max(1e-6,M(:,:,k));
        end
    else
        for k =1:N
            M(:,:,k) = (data_r(:,k)*A(:,k)' + lambda_m*M0*diag(psi_maps(:,k))) ...   %
                /(A(:,k)*A(:,k)' + lambda_m*eye(P));
            M(:,:,k) = max(1e-6,M(:,:,k));
        end
    end

    
        
    if verbose
        fprintf('Done!\n')
    end
    
    %% A_update
    
    if verbose
        fprintf('updating A...\n')
    end
    
    
	if any(lambda_a) % superpixels spatial regularization
        
        % Average all endmember matrices inside each superpixel
        if flag_use_parfor
            parfor k=1:numSuperpixels
                idx = find(spSegs(:) == k-1);
                Mc(:,:,k) = mean(M(:,:,idx), 3);
            end
        else
            for k=1:numSuperpixels
                idx = find(spSegs(:) == k-1);
                Mc(:,:,k) = mean(M(:,:,idx), 3);
            end
        end
        
        

        
        % Unmix the image in the coarse domain
        if flag_use_parfor
            parfor k=1:numSuperpixels
                % WITH OK ASC
                Stemp = [Mc(:,:,k); sqrt(lambda_a*rho_a)*eye(P); delta_asc*ones(1,P)];
                % Ac(:,k) = fnnls(Stemp'*Stemp, Stemp'*[Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], 1e-8);
                % Ac(:,k) = nnlsm_blockpivot( Stemp, [Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], false, Ac_old(:,k) );
                [Ac(:,k), ~, P_old_fnnls_c(k,:), Z_old_fnnls_c(k,:)] = fnnlsb(Stemp'*Stemp, Stemp'*[Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], P_old_fnnls_c(k,:), Z_old_fnnls_c(k,:), 1e-8);
                % opt = solopt; opt.maxtime = 2000; opt.verbose = 0; opt.algo = 'PLB'; out = bbnnls(Stemp, [Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], Ac_old(:,k), opt); Ac(:,k) = out.x;
            end
            
        else
            for k=1:numSuperpixels
                % WITH OK ASC
                Stemp = [Mc(:,:,k); sqrt(lambda_a*rho_a)*eye(P); delta_asc*ones(1,P)];
                % Ac(:,k) = fnnls(Stemp'*Stemp, Stemp'*[Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], 1e-8);
                % Ac(:,k) = nnlsm_blockpivot( Stemp, [Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], false, Ac_old(:,k) );
                [Ac(:,k), ~, P_old_fnnls_c(k,:), Z_old_fnnls_c(k,:)] = fnnlsb(Stemp'*Stemp, Stemp'*[Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], P_old_fnnls_c(k,:), Z_old_fnnls_c(k,:), 1e-8);
                % opt = solopt; opt.maxtime = 2000; opt.verbose = 0; opt.algo = 'PLB'; out = bbnnls(Stemp, [Yc(:,k); zeros(P,1); delta_asc*sumW(1,k)], Ac_old(:,k), opt); Ac(:,k) = out.x;
            end
        end

        
        
        
        % Compte the conjugate superpixel transformation
        AcWast = compute_W_superpixels(Ac, m, n, P, spSegs, numSuperpixels, 'Wast');
        AcWast = reshape(AcWast, m*n, P )';
        
        % Compute sum-to-one constraint for detail scale
        sumImWWast = ones(1,N) - ones(1,P) * AcWast;
        
        

        
        
        if flag_use_parfor
            parfor k =1:N
                Stemp = [M(:,:,k); sqrt(lambda_a)*eye(P); delta_asc*ones(1,P)];
                Ytemp = Mc(:,:,spSegs(k)+1) * Ac(:,spSegs(k)+1);
                [Ad(:,k), ~, P_old_fnnls_d(k,:), Z_old_fnnls_d(k,:)] = fnnlsb(Stemp'*Stemp, Stemp'*[(Yd_mtx(:,k)+Ytemp); sqrt(lambda_a)*Ac(:,spSegs(k)+1); delta_asc], P_old_fnnls_d(k,:), Z_old_fnnls_d(k,:), 1e-8);
            end
        else
            for k =1:N
                Stemp = [M(:,:,k); sqrt(lambda_a)*eye(P); delta_asc*ones(1,P)];
                Ytemp = Mc(:,:,spSegs(k)+1) * Ac(:,spSegs(k)+1);
                [Ad(:,k), ~, P_old_fnnls_d(k,:), Z_old_fnnls_d(k,:)] = fnnlsb(Stemp'*Stemp, Stemp'*[(Yd_mtx(:,k)+Ytemp); sqrt(lambda_a)*Ac(:,spSegs(k)+1); delta_asc], P_old_fnnls_d(k,:), Z_old_fnnls_d(k,:), 1e-8);
            end
        end
        
        
        
        % Invert superpixel trasnformation (not necessary since we solve the problem directly for Ad+AcWast)
        A = Ad;


        % Store old variables
        AcWast_old = AcWast;
        Ad_old = Ad;
        Ac_old = Ac;


        
    else % FCLSU (without spatial regularization)  
        for k =1:N
            A(:,k) =  FCLSU(data_r(:,k),M(:,:,k));
        end
	end
    
    if verbose
        fprintf('Done!\n')
    end
    
    %% psi_update
    
    if verbose
        fprintf('updating psi...\n')
    end
    
    if any(lambda_psi) % with spatial regularization
        
        if scalar_lambda_psi %% update with scalar or vector lambda_m (done in the Fourier domain)
            if flag_use_parfor
                parfor p = 1:P
                    numerator = reshape(lambda_m*permute(M(:,p,:),[1 3 2])'*M0(:,p),m,n);
                    psi_maps_im =  real(ifft2(fft2(numerator)./((lambda_psi*(abs(FDh).^2+ abs(FDv).^2) + lambda_m * M0ptM0(p)))));
                    psi_maps(p,:) = psi_maps_im(:);
                end
            else
                for p = 1:P
                    numerator = reshape(lambda_m*permute(M(:,p,:),[1 3 2])'*M0(:,p),m,n);
                    psi_maps_im =  real(ifft2(fft2(numerator)./((lambda_psi*(abs(FDh).^2+ abs(FDv).^2) + lambda_m * M0ptM0(p)))));
                    psi_maps(p,:) = psi_maps_im(:);
                end
            end
        else
            if flag_use_parfor
                parfor p = 1:P
                    numerator = reshape(lambda_m*squeeze(M(:,p,:))'*M0(:,p),m,n);
                    psi_maps_im =  real(ifft2(fft2(numerator)./((lambda_psi(p)*(abs(FDh).^2+ abs(FDv).^2) + lambda_m * M0ptM0(p)))));
                    psi_maps(p,:) = psi_maps_im(:);
                end
            else
                for p = 1:P
                    numerator = reshape(lambda_m*squeeze(M(:,p,:))'*M0(:,p),m,n);
                    psi_maps_im =  real(ifft2(fft2(numerator)./((lambda_psi(p)*(abs(FDh).^2+ abs(FDv).^2) + lambda_m * M0ptM0(p)))));
                    psi_maps(p,:) = psi_maps_im(:);
                end
            end
        end
    else % without spatial regularization
        if flag_use_parfor
            for p=1:P
                psi_maps_temp = zeros(N,1);
                parfor k = 1:N
                    psi_maps_temp(k) =  (M0(:,p)'*M(:,p,k))/M0ptM0(p);
                end
                psi_maps(p,:) = psi_maps_temp;
            end
        else 
            for p=1:P
                psi_maps_temp = zeros(N,1);
                for k = 1:N
                    psi_maps_temp(k) =  (M0(:,p)'*M(:,p,k))/M0ptM0(p);
                end
                psi_maps(p,:) = psi_maps_temp;
            end
        end
    end
    
    if verbose
        fprintf('Done!\n')
    end
    

    % =================================================================
    % residuals of the ANLS loops
    
    rm_vect = zeros(N,1);
    for k=1:N
        rm_vect(k) = norm(M(:,:,k)-M_old(:,:,k),'fro')/norm(M_old(:,:,k),'fro');
    end
    
    rm(i) = mean(rm_vect);
    ra(i) = norm(A(:)-A_old_anls(:),2)/norm(A_old_anls(:),2);
    rpsi(i) = norm(psi_maps-psi_maps_old,'fro')/(norm(psi_maps_old,'fro'));
    

    % termination test
    if verbose
            fprintf('iteration %d of %d, rm = %f, ra = %f, rpsi= %f\n',i,maxiter_anls,rm(i),ra(i),rpsi(i));
    end

    if ((rm(i) < epsilon_s) && (ra(i) < epsilon_a)) && (rpsi(i)< epsilon_psi)
        break;
    end
    
end

end

