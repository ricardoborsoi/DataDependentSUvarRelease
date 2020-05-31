% =========================================================================
% Example with the Cuprite dataset
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


%% load data        
load cuprite_ref 

data = reshape(x',Lines,Columns,L);
[m,n,L] = size(data);

load Cvec_mod % This vector attempts to compensate the image for illumination variations, useful if using spectral libraries (see Iordache et al. 2011 paper.)
data_r = reshape(data,m*n,L)*diag(c);
N = m*n;

%% load endmembers
load endmembers_cuprite
P =  14;


% ===========================================================
%% Fully Constrained Least Squares Unmixing (FCLSU)
disp('FCLSU...')
tic
A_FCLSU = FCLSU(data_r',M0)';
time_fcls = toc;


% ===========================================================
%% Scaled version of the (partially) Constrained Least Squares (SCLSU)
disp('SCLSU...')

tic
[A_SCLSU, psis_SCLSU] = SCLSU(data, M0);
time_scls = toc;

M_SCLSU = zeros(L,P,m*n);
for i = 1:m*n
   M_SCLSU(:,:,i) = M0*diag(psis_SCLSU(:,i)); 
end


% ===========================================================
%% Extended Linear Mixing Model
disp('ELMM')

% initialization with S-CLSU
A_init = A_SCLSU;
psis_init = ones(size(A_init));

% regularization parameters (same parameters as in the original paper):
lambda_m   = 0.4; 
lambda_a   = 0.005;
lambda_psi = 0.005;

% optional parameters
nnorm = '2,1';
verbose = true; % display

maxiter_anls = 100;
maxiter_admm = 100;
epsilon_m = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);

tic
[A_ELMM, psis_ELMM, M_ELMM, optim_struct] = ELMM_ADMM(data, A_init, psis_init, M0,lambda_m,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_m,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel);
time_elmm = toc;



% ===========================================================
%% Perturbed Linear Mixing Model
disp('PLMM')

A_init = A_SCLSU;

% regularization parameters:
alpha_plmm = 3.1e-4;
beta_plmm  = 5e2;
gamma_plmm = 1;

tic
[A_PLMM, dM_lexico2, M_PLMM] = interface_PLMM(data, A_init, M0, alpha_plmm, beta_plmm, gamma_plmm);
time_plmm = toc;





% ===========================================================
%% MUA-SV algorithm
disp('MUA-SV')

A_init = A_SCLSU;
psis_init = ones(size(A_init));

verbose = true; % display

lambda_m   = 5;
lambda_psi = 0.01;
lambda_a   = 0.05;
slic_size  = 6;
slic_reg   = 0.001;
rho_a      = 0.01;

epsilon_m    = 2e-3;
epsilon_a    = 2e-3;
epsilon_psi  = 2e-3;
maxiter_anls = 100;

tic
[A_MUASV, psis_MUASV, M_MUASV] = MUASV(data, A_init, psis_init, M0,lambda_m,lambda_a,lambda_psi,slic_size,slic_reg,rho_a,verbose,maxiter_anls,epsilon_m,epsilon_a,epsilon_psi);
time_MUASV = toc;



% ===========================================================
%%

A_FCLSU_im = reshape(A_FCLSU',m,n,P);
A_SCLSU_im = reshape(A_SCLSU',m,n,P);
A_ELMM_im  = reshape(A_ELMM',m,n,P);
A_PLMM_im  = reshape(A_PLMM',m,n,P);
A_MUASV_im = reshape(A_MUASV',m,n,P);

P_sel = [1 2 7 12];
P = length(P_sel);

fh = figure;
[ha, pos] = tight_subplot(5, P, 0.01, 0.1, 0.1);
for i=1:P,
    axes(ha(i));
    imagesc(A_FCLSU_im(:,:,P_sel(i)), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+P));
    imagesc(A_SCLSU_im(:,:,P_sel(i)), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+2*P));
    imagesc(A_PLMM_im(:,:,P_sel(i)), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+3*P));
    imagesc(A_ELMM_im(:,:,P_sel(i)), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+4*P));
    imagesc(A_MUASV_im(:,:,P_sel(i)), [0 1])
    set(gca,'ytick',[],'xtick',[])
end
% set(fh, 'Position', [0 0 650 700])
axes(ha(1));
title('Alunite','interpreter','latex')
axes(ha(2));
title('Sphene','interpreter','latex')
axes(ha(3));
title('Buddingtonite','interpreter','latex')
axes(ha(4));
title('Muscovite','interpreter','latex')


axes(ha(0*P+1));
ylabel('FCLS','interpreter','latex')
axes(ha(1*P+1));
ylabel('SCLS','interpreter','latex')
axes(ha(2*P+1));
ylabel('PLMM','interpreter','latex')
axes(ha(3*P+1));
ylabel('ELMM','interpreter','latex')
axes(ha(4*P+1));
ylabel('MUA-SV','interpreter','latex')
colormap jet





%%

fprintf('\n\n Time FCLS: \t %f \n', time_fcls);
fprintf('Time SCLS: \t %f \n', time_scls);
fprintf('Time ELMM: \t %f \n', time_elmm);
fprintf('Time PLMM: \t %f \n', time_plmm);
fprintf('Time MUASV: \t %f \n\n', time_MUASV);


N = m*n;

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












