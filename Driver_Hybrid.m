%% Hybrid methods
% This script for estimating fluxes of case study for low noise level(or variance)
% using hybrid methods(genHyBRs, genHyBRmean)
%
% Latest Update by T.Cho, Nov. 28, 2021


%% User needs to defined the following options in this section.
% Add path

    addpath(genpath('genHyBR-master'));
    rmpath('genHyBR-master/RestoreTools/IterativeMethods');
    
% There are two case studies

    CaseStudy = '6wk';
    % CaseStudy = '1yr';

% Choose hybrid methods: {'genhybr', 'meanest'}
%   genhybr: genHyBRs in paper
%   meanest: genHyBRmean in paper    

    HybridOptions = 'meanest';
    % maximum iteration of hybrid method
    iter = 50; 
    
% Choose noise level: 5%, 10%, 50% (noise ratio to observation)
% as described in paper

    NoiseOptions = [5, 10, 50]; % one or more options can be chosen
    nNO = length(NoiseOptions);

% Choose regularization option : {'optimal', 'dp', 'LSQR'}
%   optimal: find optimal regularization parameter using exact s
%   dp: discrepancy principle
%   LSQR: no regularization parameter, zero(0)

    RegOptions = {'optimal', 'dp', 'LSQR'}; % one or more options can be chosen
    nRO = length(RegOptions);
    
% Setup path and filenames for saving outputs
    
    outpath = "Define path where output of this script will be saved"
    
    outfile_fluxes = cell(nNO,nRO);
    outfile_output = cell(nNO,nRO);
    
    for j = 1:nRO
        for i = 1:nNO
            outfile_fluxes{i,j} = strcat(outpath,'Fluxes_',num2str(NoiseOptions(i)),'_',RegOptions(j),'_',date,'.csv');
            outfile_output{i,j} = strcat(outpath,'Output_',num2str(NoiseOptions(i)),'_',RegOptions(j),'_',date,'.mat');
        end
    end
    
    outfile_fluxes = outfile_fluxes(:);
    outfile_output = outfile_output(:);
%%

% Choose number of times
   
    switch CaseStudy
        case '6wk'
            ntimes = 328;
        case '1yr'
            ntimes = 2920;
    end

% Read in the path of H matrix and s matrix
    
    Hpath = strcat(CaseStudy,'/H/');
    spath = strcat(CaseStudy,'/s/');
    
% Create the X matrix only if meanest   
    
    disp('Create the X matrix');
    switch HybridOptions
        case 'meanest'
            switch CaseStudy
                case '6wk' % Get 3-hourly pattern
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
                    
                case '1yr' % Get monthly pattern
                    m1 = 3222;
                    m = 9408240;
                    X = sparse(m,12);
                    mm = m./12;
                    X(1:mm,1) = 1;
                    X((mm+1):(2.*mm),2) = 1;
                    X((2.*mm +1):(3.*mm),3) = 1;
                    X((3.*mm +1):(4.*mm),4) = 1;
                    X((4.*mm +1):(5.*mm),5) = 1;
                    X((5.*mm +1):(6.*mm),6) = 1;
                    X((6.*mm +1):(7.*mm),7) = 1;
                    X((7.*mm +1):(8.*mm),8) = 1;
                    X((8.*mm +1):(9.*mm),9) = 1;
                    X((9.*mm +1):(10.*mm),10) = 1;
                    X((10.*mm +1):(11.*mm),11) = 1;
                    X((11.*mm +1):(12.*mm),12) = 1;
            end
    end
    
% Load distmat.mat from Hpath    

    load(strcat(Hpath,'distmat.mat'));

% Get True s
    
    s_true = [];
    for j = 1:ntimes
        load(strcat(spath,'s_',num2str(j),'.mat'));
        s_true = [s_true; s(:)]; 
    end
    
 % get True(not exact) beta for meanest. Used "Get_Opt_Beta.m"
    
    switch HybridOptions
        case 'meanest'
            switch CaseStudy
                case '6wk' % beta has 8 elements
                    beta_true = [2.6327, 3.5117, 3.4709, 2.7027, -1.7343, -6.7710, -6.6152, -1.6639];
                    s_true = [s_true(:); beta_true(:)];                    
                case '1yr'% beta has 12 elements
                    beta_true = [0.3568, 0.3240, 0.4357, 0.4217, 0.4265, 0.3501, 0.2961, 0.2327, -0.0163, -0.2017, -0.3655, -0.2170];
                    s_true = [s_true(:); beta_true(:)];                    
            end
    end
    
 % Setting for parallel computing   
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


%% Parallel Computing

    
    parfor k = 1 : (nNO*nRO)

        % Read in the observation vector here
        Zmat = load(strcat(spath,'Z_',num2str(parNoiseOptions(k)),'.mat'));
        Z = Zmat.Zn;

        % Set the covariance matrix parameters for the inversion
        switch CaseStudy
            case '6wk'
                theta = [ 2.000 10.046 555.420 9.854 ];
            case '1yr'
                theta = [ 2.000 10.046 585.680 12.366 ];
        end
        theta(1) = Zmat.sigma;
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
        
        switch CaseStudy
            case '1yr'
                switch HybridOptions
                    case 'meanest'
                        D_theta = [ sqrt(34.7) sqrt(14.84) sqrt(16.09) sqrt(8.32) sqrt(6.66) sqrt(7.47) sqrt(9.67) sqrt(9.37) sqrt(39.07) sqrt(72.63) sqrt(102.37) sqrt(70.43) ]';
                        sigmaQ  = [ D_theta(1) .* ones(8.*30,1); D_theta(2) .* ones(8.*31,1); D_theta(3) .* ones(8.*30,1); ...
                                    D_theta(4) .* ones(8.*31,1); D_theta(5) .* ones(8.*31,1); D_theta(6) .* ones(8.*28,1); ...
                                    D_theta(7) .* ones(8.*31,1); D_theta(8) .* ones(8.*30,1); D_theta(9) .* ones(8.*31,1); ...
                                    D_theta(10) .* ones(8.*30,1); D_theta(11) .* ones(8.*31,1); D_theta(12) .* ones(8.*31,1)];
                        % scaling sigmaQ
                        sigmaQ = sigmaQ/max(sigmaQ(:));  
                        % scaling D
                        D = (sigmaQ*sigmaQ') .* D;
                end
        end
        
        % Create the R matrix   
        disp('Create R');
        R = ones(n,1);
        R = spdiags(R,0,n,n);
        
        % Create augmented Q and H 
        disp('Create H and Q-kron');
        
        switch HybridOptions
            case 'genhybr'
                Q = kronMat(D,E);
                H = matvecH(ntimes,Hpath);
            case 'meanest'
                nbeta = size(X,2);
                gamma = 10;
                Qbeta = (1/gamma^2) * eye(nbeta);
                Q = kronMat(D,E);
                Q = QtilMat(Qbeta,Q,X);
                H = matvecH_aug(ntimes,Hpath,nbeta);
        end
        
        % Input for hybrid methods
        switch char(parRegOptions(k))
            case 'LSQR'
                disp(strcat('Run_', char(parRegOptions(k)), '_for_', num2str(parNoiseOptions(k)), '%_noise_level' ));
                input = HyBR_lsmrset('InSolv', 'none',...
                             'Reorth','on','Iter', iter,'x_true',s_true,'nLevel',theta(1));
            otherwise
                disp(strcat('Run_', char(parRegOptions(k)), '_for_', num2str(parNoiseOptions(k)), '%_noise_level' ));
                input = HyBR_lsmrset('InSolv', 'tikhonov','RegPar',char(parRegOptions(k)),...
                             'Reorth','on','Iter', iter,'x_true',s_true,'nLevel',theta(1));
        end
        
        % Run Main algorithm
        tic;
        [recons, output] = genHyBR(H, Z(:), Q, R, input);
        toc;
        
        % Write the outputs to csv file
        disp('Writing outputs to file');
        outfile_f = char(outfile_fluxes{k});
        dlmwrite(outfile_f,full(recons),',');
        disp('Outputs written to file');
        disp(outfile_f)
        
        % Choose Output components to remove from save
%         output = rmfield(output,'U');
%         output = rmfield(output,'V');
        output = rmfield(output,'QV');
%         output = rmfield(output,'B');

        parsave(char(outfile_output{k}),output);
        disp(char(outfile_output{k}))
    end
    
