%% genHyBRmean
% This script for testing 6wk case study for low noise level(or variance)
% with genHyBRmean (iterative + regularization + mean estimation)
%
% Test:     nlevel      sigma       observation file
%               5%      0.0612      Z_5.mat 
%              10%      0.1224      Z_10.mat 
%              50%      0.6121      Z_50.mat 
%             100%      1.2243      Z_100.mat 
%
% Note that the Z with 2 ppm noise variance 
%         163.365%      2.0000      Z.mat 
%
% For each noise level, tested genHyBR w/ regularization parameter selected
% by 'optimal', 'dp' ,'wgcv', 'upre' and 'LSQR'
%
% Latest Update by T.Cho, Apr. 7, 2021


% Only need to change this section.
% Set script options
iter = 200;
NoiseOptions = [5, 10, 50];
SigmaOptions = [0.0612, 0.1224, 0.6121];
RegOptions = {'optimal','dp','wgcv','upre','LSQR'};

%%
nNO = length(NoiseOptions);
nRO = length(RegOptions);

% Setup path and filenames for saving outputs
outpath = 'output/';
outfile_fluxes = cell(nNO,nRO);
outfile_error = cell(nNO,nRO);
outfile_res = cell(nNO,nRO);
outfile_output = cell(nNO,nRO);

for j = 1:nRO
    for i = 1:nNO
        outfile_fluxes{i,j} = strcat(outpath,'Fluxes_meanEst_',num2str(NoiseOptions(i)),'_',RegOptions(j),'_',date,'_',num2str(iter),'.csv');
        outfile_error{i,j} = strcat(outpath,'Error_meanEst_',num2str(NoiseOptions(i)),'_',RegOptions(j),'_',date,'_',num2str(iter),'.csv');
        outfile_res{i,j} = strcat(outpath,'Res_meanEst_',num2str(NoiseOptions(i)),'_',RegOptions(j),'_',date,'_',num2str(iter),'.csv');
        outfile_output{i,j} = strcat(outpath,'Output_meanEst_',num2str(NoiseOptions(i)),'_',RegOptions(j),'_',date,'_',num2str(iter),'.mat');
    end
end

outfile_output = outfile_output(:);
outfile_fluxes = outfile_fluxes(:);
outfile_error = outfile_error(:);
outfile_res = outfile_res(:);


% Choose number of times
ntimes = 328;

% Read in the path of H matrix and s matrix
Hpath = 'H_6wk/';
spath = 's_6wk/';

% Create the X matrix
disp('Create the X matrix');
X = [] ;
    for i = 1:8
        for j = 1:3222*41*8
                if rem((fix((j-1)/3222) + 1 - i), 8) == 0
                X(j, i) = 1 ;
                else
                X(j, i) = 0 ;
                end
        end
    end	
X = sparse(X);

% Load distmat.mat
load('H_6wk/distmat.mat');

s_true = [];
for j = 1:ntimes
    load(strcat(spath,'s_',num2str(j),'.mat'));
    s_true = [s_true; s(:)]; 
end
beta_true = [2.5244 3.8394 3.8928 3.1219 -1.8360 -7.5656 -7.7407 -2.5368]';
x_true = [s_true(:); beta_true(:)];

parNoiseOptions = [];
for i = 1:nRO
    parNoiseOptions = [parNoiseOptions NoiseOptions(:)];
end
parNoiseOptions = parNoiseOptions(:);

parRegOptions = [];
for i = 1:nNO
   parRegOptions = [parRegOptions; RegOptions];
end
parRegOptions = parRegOptions(:);

%% Parallel computing
parfor k = 1 : (nNO*nRO)
    
    % Read in the observation vector here
    Zmat = load(strcat(Hpath,'Z_',num2str(parNoiseOptions(k)),'.mat'));
    Z = Zmat.Zn;
    
    % Set the covariance matrix parameters for the inversion
    theta = [ 2.000 10.046 555.420 9.854 ];
    theta(1) = Zmat.sigma_t;
    disp('Covariance matrix parameters');
    disp(theta);
    
    n = length(Z);
    
    % Create E
    % Spherical covariance model
    E                      = 1 - 1.5 .* (deltamat ./theta(3))  + 0.5 .* (deltamat.^3 ./ theta(3).^3);
    E(deltamat > theta(3)) = 0;
    
    % Create D
    % Create time distance matrix
    days = [];
    for i = 1:ntimes
        for j = 1:ntimes
            if rem(abs(i-j),8) == 0
                days(i,j) = abs(i-j)/8;
            else
                days(i,j) = 10^6;
            end
        end
    end

    % Spherical covariance model
    D = 1 - 1.5 .* (days   ./theta(4))  + 0.5 .* (days.^3   ./ theta(4).^3);
    D(days   > theta(4)) = 0;

    % Create the R matrix   
    disp('Create R');
    R = ones(n,1);
	R = spdiags(R,0,n,n);
    
    
    % Create augmented Q and H 
    disp('Create H and Q-kron');
    nbeta = size(X,2);
    gamma = 10;
    Qbeta = (1/gamma^2) * eye(nbeta);
    Qkron = kronMat(D,E);
    Qtil = QtilMat(Qbeta,Qkron,X);
    H_aug = matvecH_aug(ntimes,Hpath,nbeta);

    % Set up iterative methods
    % Low noise level usually has semi-convergence later than higher noise.
    % Giving different maximum iterations for iterative methods
    switch parNoiseOptions(k)
        case 5
            maxit = iter;
        case 10
            maxit = iter/2;
        case 50
            maxit = iter/6;
    end
    
    switch char(parRegOptions(k))
        case 'LSQR'
            disp(strcat('Run_', char(parRegOptions(k)), '_for_', num2str(parNoiseOptions(k)), '%_noise_level' ));
            input = HyBR_lsmrset('InSolv', 'none',...
                         'Reorth','on','Iter', maxit,'x_true',x_true,'nLevel',theta(1));
        otherwise
            disp(strcat('Run_', char(parRegOptions(k)), '_for_', num2str(parNoiseOptions(k)), '%_noise_level' ));
            input = HyBR_lsmrset('InSolv', 'tikhonov','RegPar',char(parRegOptions(k)),...
                         'Reorth','on','Iter', maxit,'x_true',x_true,'nLevel',theta(1));
    end
    
    % Run genHyBR for given noise and regulairzation parameter
    % selection methods
    tic;
    [recons, output] = genHyBR(H_aug, Z(:), Qtil, R, input);
    toc;

    % Write the outputs to csv file
    disp('Writing outputs to file');
    outfile_f = char(outfile_fluxes{k});
    outfile_e = char(outfile_error{k});
    outfile_r = char(outfile_res{k});

    dlmwrite(outfile_f,full(recons),',');
    dlmwrite(outfile_e,full(output.Enrm),',');
    dlmwrite(outfile_r,full(output.Rnrm),',');

    disp('Outputs written to file');
    disp(outfile_f)
    disp(outfile_e)
    disp(outfile_r)
    
    % Choose Output components to save
    output = rmfield(output,'U');
    output = rmfield(output,'V');
    output = rmfield(output,'QV');
    output = rmfield(output,'B');
    
    parsave(char(outfile_output{k}),output);
    disp(char(outfile_output{k}))
end