% =========================================================================
% Example with the Houston dataset, provided by Lucas Drumets and his collaborators
% 
% This code corresponds to example 3 in the following publication:
%   "A data dependent multiscale model for hyperspectral unmixing with spectral variability"
%   Ricardo Augusto Borsoi, Tales Imbiriba, José Carlos Moreira Bermudez
%   IEEE Transactions on Image Processing, v. 29, 3638-3651, 2020.
% =========================================================================


clear
close all
clc

addpath DATA
addpath MUASV
addpath(genpath('other_methods'))
addpath(genpath('utils'))



% Fix RNG just in case for reproducibility
rng(5, 'twister')


%% load data
load data_houston 

[m,n,L] = size(data);
N = m*n;
data_r = reshape(data,m*n,L);

load endmembers_houston % load initial reference endmembers
M0 = S0;
P = size(M0,2);

% identify materials
%materials{1} = 'vegetation'; 
%materials{2} = 'red roofs';
%materials{3} = 'concrete';
%materials{4} = 'asphalt';



%% Fully Constrained Least Squares Unmixing (FCLSU)
disp('FCLSU...')
tic
A_FCLSU = FCLSU(data_r',M0)';
time_fcls = toc;



%% Scaled version of the (partially) Constrained Least Squares (SCLSU)
disp('S-CLSU...')

tic
[A_SCLSU, psis_SCLSU] = SCLSU(data, M0);
time_scls = toc;

M_SCLSU = zeros(L,P,N);
for i = 1:m*n
   M_SCLSU(:,:,i) = M0*diag(psis_SCLSU(:,i)); 
end



% ----------------------------------------------------------
%% Full Extended Linear Mixing Model
disp('ELMM')

% initialization with S-CLSU
A_init = A_SCLSU;
psis_init = ones(size(A_init));

% regularization parameters (same parameters as in the paper [1]):
lambda_m   = 0.4; 
lambda_a   = 0.005;
lambda_psi = 0.001;

% optional parameters
nnorm = '1,1'; % Use a Total Variation on the abundances
verbose = true; % display

maxiter_anls = 100;
maxiter_admm = 100;
epsilon_m = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 5e-3;
epsilon_admm_rel = 5e-3;

tic
[A_ELMM, psis_ELMM, M_ELMM, optim_struct] = ELMM_ADMM(data, A_init, psis_init, M0,lambda_m,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_m,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel);
time_elmm = toc;



% ----------------------------------------------------------
%% Perturbed Linear Mixing Model
disp('PLMM')
A_init = A_SCLSU;

% regularization parameters (same parameters as in the paper [1]):
alpha_plmm = 1.4e-3;
beta_plmm  = 5e2;
gamma_plmm = 1;

tic
[A_PLMM, dM_lexico2, M_PLMM] = interface_PLMM(data, A_init, M0, alpha_plmm, beta_plmm, gamma_plmm);
time_plmm = toc;



% ----------------------------------------------------------
%% MUA-SV algorithm
disp('MUA-SV')

A_init = A_SCLSU;
psis_init = ones(size(A_init));

verbose = true; % display
lambda_m   = 0.5;
lambda_psi = 0.001;
lambda_a   = 0.001;
slic_size  = 5;
slic_reg   = 0.001;
rho_a      = 0.35;

epsilon_m    = 2e-3;
epsilon_a    = 2e-3;
epsilon_psi  = 2e-3;
maxiter_anls = 100;

tic
[A_MUASV, psis_MUASV, M_MUASV] = MUASV(data, A_init, psis_init, M0,lambda_m,lambda_a,lambda_psi,slic_size,slic_reg,rho_a,verbose,maxiter_anls,epsilon_m,epsilon_a,epsilon_psi);
time_MUASV = toc;







%% ==============================================
% abundance maps

A_FCLSU_im = reshape(A_FCLSU',m,n,P);
A_SCLSU_im = reshape(A_SCLSU',m,n,P);
A_ELMM_im  = reshape(A_ELMM',m,n,P);
A_PLMM_im  = reshape(A_PLMM',m,n,P);
A_MUASV_im = reshape(A_MUASV',m,n,P);

fh = figure;
[ha, pos] = tight_subplot(5, 4, 0.01, 0.1, 0.1);
for i=1:4,
    axes(ha(i));
    imagesc(A_FCLSU_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+4));
    imagesc(A_SCLSU_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+8));
    imagesc(A_PLMM_im(:,:,i), [min(A_PLMM_im(:)), max(A_PLMM_im(:))])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+12));
    imagesc(A_ELMM_im(:,:,i), [min(A_ELMM_im(:)), max(A_ELMM_im(:))])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+16));
    imagesc(A_MUASV_im(:,:,i), [min(A_MUASV_im(:)), max(A_MUASV_im(:))])
    set(gca,'ytick',[],'xtick',[])
end
set(fh, 'Position', [0 0 650 700])
axes(ha(1));
ylabel('FCLS','interpreter','latex')
title('Vegetation','interpreter','latex')
axes(ha(2));
title('Met. Roofs','interpreter','latex')
axes(ha(3));
title('Concrete','interpreter','latex')
axes(ha(4));
title('Asphalt','interpreter','latex')
axes(ha(5));
ylabel('SCLS','interpreter','latex')
axes(ha(9));
ylabel('PLMM','interpreter','latex')
axes(ha(13));
ylabel('ELMM','interpreter','latex')
axes(ha(17));
ylabel('MUA-SV','interpreter','latex')
colormap jet







%%

fprintf('\n\n Time FCLS: \t %f \n', time_fcls);
fprintf('Time SCLS: \t %f \n', time_scls);
fprintf('Time ELMM: \t %f \n', time_elmm);
fprintf('Time PLMM: \t %f \n', time_plmm);
fprintf('Time MUASV: \t %f \n\n', time_MUASV);



% Reconstruction errors
H_FCLSU = M0*A_FCLSU; % reconstruction for S-FCLSU
H_SCLSU = zeros(L,N); % reconstruction for S-CLSU
H_ELMM  = zeros(L,N); % reconstruction for ELMM
H_PLMM  = zeros(L,N); % reconstruction for PLMM
H_MUASV = zeros(L,N); % reconstruction for MUASV

for i=1:N
   H_ELMM(:,i)  = squeeze(M_ELMM(:,:,i))*A_ELMM(:,i); 
   H_SCLSU(:,i) = squeeze(M_SCLSU(:,:,i))*A_SCLSU(:,i);
   H_PLMM(:,i)  = squeeze(M_PLMM(:,:,i))*A_PLMM(:,i);
   H_MUASV(:,i) = squeeze(M_MUASV(:,:,i))*A_MUASV(:,i); 
end

RMSE_FCLSU = 1/L*sum((H_FCLSU'-data_r).^2,2);
RMSE_SCLSU = 1/L*sum((H_SCLSU'-data_r).^2,2);
RMSE_ELMM  = 1/L*sum((H_ELMM'-data_r).^2,2);
RMSE_PLMM  = 1/L*sum((H_PLMM'-data_r).^2,2);
RMSE_MUASV  = 1/L*sum((H_MUASV'-data_r).^2,2);

fprintf('\n\nReconstruction error FCLS: \t\t %f \n', sum(RMSE_FCLSU)/N );
fprintf('Reconstruction error SCLS: \t\t %f \n',     sum(RMSE_SCLSU)/N  );
fprintf('Reconstruction error ELMM: \t\t %f \n',     sum(RMSE_ELMM)/N  );
fprintf('Reconstruction error PLMM: \t\t %f \n',     sum(RMSE_PLMM)/N  );
fprintf('Reconstruction error MUASV: \t\t %f \n\n\n', sum(RMSE_MUASV)/N  );





