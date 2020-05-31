% =========================================================================
% Example with synthetic data generated according to the Hapke model
% Important: the data for this example was originally presented at the ELMM paper
%            by Lucas Drumetz and his collaborators.
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


rng(5, 'twister') 


% Noise SNR
SNR = 40;


% =========================================================================

% load data
load DATA/test_est_ex_hapke.mat

L  = 16;
P  = 3;
nr = 50;
nc = 50;
N  = nr*nc;

% Generate image
data_r = zeros(L, N);
for i=1:N
    data_r(:,i) = Mn_true(:,:,i) * Ath(:,i);
end
data_r = data_r';

% Add noise
sigma = sqrt(sum(sum((data_r).^2))/N/L/10^(SNR/10));
data_r = data_r + sigma*randn(N,L);

% Create data cube
data   = ilexico(data_r, nr, nc, 'col');





% Compute M0 --------------------------------
% Endmember initialization (VCA [1])
M0 = vca(data_r','Endmembers',P,'verbose','off');

% Sort M0 with respect to real/desired EM signatures to ease the comparison of estimated abundance maps
id = zeros(P,1);
for k = 1:P
    for l = 1:P
        s(l) = 180*acos( (Mth(:,k).')*M0(:,l) /(norm(Mth(:,k))*norm(M0(:,l))) )/pi; 
    end
    [~, id(k)] = min(s);
end
M0 = M0(:,id);


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
for i = 1:nr*nc
   M_SCLSU(:,:,i) = M0*diag(psis_SCLSU(:,i)); 
end



% ---------------------------------------------------
%% Full Extended Linear Mixing Model
disp('ELMM')

% initialization with S-CLSU
A_init = A_SCLSU;
psis_init = ones(size(A_init));

% select regularization parameters
switch SNR
    case 40
        lambda_m   = 0.1;
        lambda_a   = 0.0005;
        lambda_psi = 0.5; 
    case 30
        lambda_m   = 0.005; 
        lambda_a   = 0.01;
        lambda_psi = 0.01; 
    case 20
        lambda_m   = 0.005; 
        lambda_a   = 0.01;
        lambda_psi = 0.0005;
    otherwise
        error('ELMM parameters not defined for specified SNR! value')
end

% optional parameters
nnorm = '2,1'; % Use a Total Variation on the abundances
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



%--------------------------------------------------------
%% Perturbed Linear Mixing Model
disp('PLMM')

% select regularization parameters
switch SNR
    case 40
        alpha_plmm = 0.01; 
        beta_plmm  = 50;
        gamma_plmm = 1;
    case 30
        alpha_plmm = 0.01; 
        beta_plmm  = 50;
        gamma_plmm = 1;
    case 20
        alpha_plmm = 0.01; 
        beta_plmm  = 50;
        gamma_plmm = 1.5;
    otherwise
        error('ELMM parameters not defined for specified SNR! value')
end
        
A_init = A_SCLSU;

tic
[A_PLMM, dM_lexico2, M_PLMM] = interface_PLMM(data, A_init, M0, alpha_plmm, beta_plmm, gamma_plmm);
time_plmm = toc;



% ------------------------------------------------------------
%% MUA-SV algorithm
disp('MUA-SV')

run('utils/vlfeat-0.9.20/toolbox/vl_setup')

% select regularization parameters
switch SNR
    case 40
        lambda_m   = 5;
        lambda_psi = 0.01; 
        lambda_a   = 0.01;
        slic_size  = 2;
        slic_reg   = 0.0005;
        rho_a      = 0.001;
    case 30
        lambda_m   = 5;
        lambda_psi = 0.05; 
        lambda_a   = 0.01;
        slic_size  = 4;
        slic_reg   = 0.01;
        rho_a      = 0.001;
    case 20
        lambda_m   = 0.01;
        lambda_psi = 0.05; 
        lambda_a   = 0.005;
        slic_size  = 4;
        slic_reg   = 0.05;
        rho_a      = 0.001;
    otherwise
        error('MUASV parameters not defined for specified SNR! value')
end
                
maxiter_anls = 100;
epsilon_m   = 2e-3;
epsilon_a   = 2e-3;
epsilon_psi = 2e-3;

tic
[A_MUASV, psis_MUASV, M_MUASV] = MUASV(data, A_init, psis_init, M0,lambda_m,lambda_a,lambda_psi,slic_size,slic_reg,rho_a,verbose,maxiter_anls,epsilon_m,epsilon_a,epsilon_psi);
time_MUASV = toc;


%%
% =========================================================================
% Error metrics
% =========================================================================
clc

A_FCLSU_im = reshape(A_FCLSU',nr,nc,P);
A_SCLSU_im = reshape(A_SCLSU',nr,nc,P);
A_ELMM_im  = reshape(A_ELMM',nr,nc,P);
A_PLMM_im  = reshape(A_PLMM',nr,nc,P);
A_MUASV_im = reshape(A_MUASV',nr,nc,P);

fprintf('\n\nAbundance errors FCLS: \t\t %f \n', norm(Aimth(:)-A_FCLSU_im(:)) ^2/N );
fprintf('Abundance errors SCLS: \t\t %f \n',      norm(Aimth(:)-A_SCLSU_im(:)) ^2/N  );
fprintf('Abundance errors ELMM: \t\t %f \n',      norm(Aimth(:)-A_ELMM_im(:))  ^2/N  );
fprintf('Abundance errors PLMM: \t\t %f \n',      norm(Aimth(:)-A_PLMM_im(:))  ^2/N  );
fprintf('Abundance errors MUASV: \t %f \n\n\n', norm(Aimth(:)-A_MUASV_im(:)) ^2/N  );


% Spectral errors
M_SCLSU_lex2 = row2col_lexico_order(M_SCLSU, nr, nc);
M_ELMM_lex2  = row2col_lexico_order(M_ELMM, nr, nc);
M_PLMM_lex2  = row2col_lexico_order(M_PLMM, nr, nc);
M_MUASV_lex2 = row2col_lexico_order(M_MUASV, nr, nc);
SS_FCLS      = repmat(M0, [1 1 N]);

mse_FCLS = sum( (1/(N*L*P)) * (SS_FCLS(:)-Mn_true(:)).^2 );
mse_SCLS = sum( (1/(N*L*P)) * (M_SCLSU_lex2(:)-Mn_true(:)).^2 );
mse_ELMM = sum( (1/(N*L*P)) * (M_ELMM_lex2(:)-Mn_true(:)).^2 );
mse_PLMM = sum( (1/(N*L*P)) * (M_PLMM_lex2(:)-Mn_true(:)).^2 );
mse_MUASV = sum( (1/(N*L*P)) * (M_MUASV_lex2(:)-Mn_true(:)).^2 );

fprintf('Spectral MSE FCLS: \t\t %f \n', mse_FCLS);
fprintf('Spectral MSE SCLS: \t\t %f \n', mse_SCLS);
fprintf('Spectral MSE ELMM: \t\t %f \n', mse_ELMM);
fprintf('Spectral MSE PLMM: \t\t %f \n', mse_PLMM);
fprintf('Spectral MSE MUASV: \t\t %f \n\n', mse_MUASV);



% Compute spectral angles
sam_FCLS = 0;
sam_SCLS = 0;
sam_ELMM = 0;
sam_PLMM = 0;
sam_MUASV = 0;

for n=1:N
    for i=1:P
        sam_FCLS = sam_FCLS + computeSpectralAngle(SS_FCLS(:,i,n), Mn_true(:,i,n));
        sam_SCLS = sam_SCLS + computeSpectralAngle(M_SCLSU_lex2(:,i,n), Mn_true(:,i,n));
        sam_ELMM = sam_ELMM + computeSpectralAngle(M_ELMM_lex2(:,i,n), Mn_true(:,i,n));
        sam_PLMM = sam_PLMM + computeSpectralAngle(M_PLMM_lex2(:,i,n), Mn_true(:,i,n));
        sam_MUASV = sam_MUASV + computeSpectralAngle(M_MUASV_lex2(:,i,n), Mn_true(:,i,n));
    end
end

sam_FCLS = sam_FCLS/N;
sam_SCLS = sam_SCLS/N;
sam_ELMM = sam_ELMM/N;
sam_PLMM = sam_PLMM/N;
sam_MUASV = sam_MUASV/N;

fprintf('\n\nSpectral angles FCLS: \t\t %f \n', sam_FCLS);
fprintf('Spectral angles SCLS: \t\t %f \n', sam_SCLS);
fprintf('Spectral angles ELMM: \t\t %f \n', sam_ELMM);
fprintf('Spectral angles PLMM: \t\t %f \n', sam_PLMM);
fprintf('Spectral angles MUASV: \t\t %f \n\n', sam_MUASV);



% Reconstruction errors
H_FCLSU = M0*A_FCLSU; % reconstruction for S-FCLSU
H_SCLSU = zeros(L,nr*nc); % reconstruction for S-CLSU
H_ELMM  = zeros(L,nr*nc); % reconstruction for ELMM
H_PLMM  = zeros(L,nr*nc); % reconstruction for PLMM
H_MUASV = zeros(L,nr*nc); % reconstruction for MUASV

for i=1:nr*nc
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

RMSE_FCLSU_im = reshape(RMSE_FCLSU,nr,nc);
RMSE_SCLSU_im = reshape(RMSE_SCLSU,nr,nc);
RMSE_ELMM_im  = reshape(RMSE_ELMM,nr,nc);
RMSE_PLMM_im  = reshape(RMSE_PLMM,nr,nc);
RMSE_MUASV_im  = reshape(RMSE_MUASV,nr,nc);

fprintf('\n\nReconstruction error FCLS: \t\t %f \n', sum(RMSE_FCLSU)/N );
fprintf('Reconstruction error SCLS: \t\t %f \n',     sum(RMSE_SCLSU)/N  );
fprintf('Reconstruction error ELMM: \t\t %f \n',     sum(RMSE_ELMM)/N  );
fprintf('Reconstruction error PLMM: \t\t %f \n',     sum(RMSE_PLMM)/N  );
fprintf('Reconstruction error MUASV: \t\t %f \n\n\n', sum(RMSE_MUASV)/N  );



% variability maps
varPsis_SCLS  = ilexico(psis_SCLSU, nr, nc, 'col');
varPsis_ELMM  = ilexico(psis_ELMM, nr, nc, 'col');
varPsis_MUASV = ilexico(psis_MUASV, nr, nc, 'col');
vardM_PLMM    = ilexico(dM_lexico2, nr, nc, 'col');

vardM_PLMM2 = zeros(nr,nc,P);
for i=1:nr
    for j=1:nc
        for k=1:P
            vardM_PLMM2(i,j,k) = norm(squeeze(vardM_PLMM(:,k,i,j)));
        end
    end
end




% Display run time

fprintf('\n\n Time FCLS: \t %f \n', time_fcls);
fprintf('Time SCLS: \t %f \n', time_scls);
fprintf('Time ELMM: \t %f \n', time_elmm);
fprintf('Time PLMM: \t %f \n', time_plmm);
fprintf('Time MUASV: \t %f \n\n', time_MUASV);




%%
% Plot abundance maps
fh = figure;
[ha, pos] = tight_subplot(6, 3, 0.01, 0.1, 0.1);
for i=1:3,
    axes(ha(i));
    imagesc(Aimth(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+3));
    imagesc(A_FCLSU_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+6));
    imagesc(A_SCLSU_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+9));
    imagesc(A_PLMM_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+12));
    imagesc(A_ELMM_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+15));
    imagesc(A_MUASV_im(:,:,i), [0 1])
    set(gca,'ytick',[],'xtick',[])
end
set(fh, 'Position', [0 0 500 700])
axes(ha(1));
ylabel('True','interpreter','latex')
title('Endmember 1','interpreter','latex')
axes(ha(2));
title('Endmember 2','interpreter','latex')
axes(ha(3));
title('Endmember 3','interpreter','latex')
axes(ha(4));
ylabel('FCLS','interpreter','latex')
axes(ha(7));
ylabel('SCLS','interpreter','latex')
axes(ha(10));
ylabel('PLMM','interpreter','latex')
axes(ha(13));
ylabel('ELMM','interpreter','latex')
axes(ha(16));
ylabel('MUA-SV','interpreter','latex')
colormap jet

% print(['examples/estim_abundances_hapke_SNR' num2str(SNR) '_tght'],'-depsc')








