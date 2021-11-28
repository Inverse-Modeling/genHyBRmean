%% Direct method
% This script is for estimating fluxes and of case study fow low noise level(or variance)
% using direct method

% Latest Update by T.Cho, Nov. 28, 2021

%% User needs to defined the following options in this section.
% Add path

    addpath(genpath('functions'));
    addpath(genpath('DirectMethod'));
    
% Fix reg. parameter - choose regularization parameter to be used in direct
% method

    fix_reg_param = 1;

% There are two case studies

    CaseStudy = '6wk';
    % CaseStudy = '1yr';
    
% Choose noise level: 5%, 10%, 50% (noise ratio to observation)
% as described in paper
    
    NoiseOptions = 50;
    
% Setup path and filenames for saving outputs
    
    outpath = "Define path where output of this script will be saved"
    
    outfile_fluxes = strcat(outpath,'Fluxes_',num2str(NoiseOptions),'_direct_',date,'.csv');
    outfile_uncert = strcat(outpath,'Uncert_',num2str(NoiseOptions),'_direct_',date,'.csv');
    
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
    
% Create the X matrix   

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
            
% Load distmat.mat from Hpath    

    load(strcat(Hpath,'distmat.mat'));
    
%% Main algorithm
    
    % Read in the observation vector here
    Zmat = load(strcat(spath,'Z_',num2str(NoiseOptions),'.mat'));
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
    % Take the inverse of the E matrix
	Einv  = inv(E);
    
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
        
        case '6wk'
                sigmaQ = ones(size(D,1),1) / fix_reg_param;
        case '1yr'
                D_theta = [ sqrt(34.7) sqrt(14.84) sqrt(16.09) sqrt(8.32) sqrt(6.66) sqrt(7.47) sqrt(9.67) sqrt(9.37) sqrt(39.07) sqrt(72.63) sqrt(102.37) sqrt(70.43) ]';
                sigmaQ  = [ D_theta(1) .* ones(8.*30,1); D_theta(2) .* ones(8.*31,1); D_theta(3) .* ones(8.*30,1); ...
                            D_theta(4) .* ones(8.*31,1); D_theta(5) .* ones(8.*31,1); D_theta(6) .* ones(8.*28,1); ...
                            D_theta(7) .* ones(8.*31,1); D_theta(8) .* ones(8.*30,1); D_theta(9) .* ones(8.*31,1); ...
                            D_theta(10) .* ones(8.*30,1); D_theta(11) .* ones(8.*31,1); D_theta(12) .* ones(8.*31,1)];
                % scaling sigmaQ
                sigmaQ = sigmaQ/max(sigmaQ(:));  
                sigmaQ = sigmaQ / fix_reg_param;
    end
    
    % scaling D
    D = (sigmaQ*sigmaQ') .* D;
    
    % Take the inverse of D
	Dinv      =  inv(D);
    
    % Create the R matrix   
    disp('Create R');
    R = ones(n,1);
    R = spdiags(R,0,n,n);

    % Create augmented Q and H 
    disp('Create H and Q-kron');
    
    % Define the sizes of the different matrices
    p   = size(X,2);
    n   = length(Z);
    m1  = size(E,1);
    ntimes = size(D,1);
    
%% Esimate the fluxes with direct method

% Use the standard kriging equations to estimate the fluxes

	% Compute the psi matrix
    tic;
	disp('Launch the script to calculate HQH + R');
        psi = HQHR(R, D, E, Hpath);
    disp(toc);
    
    % Compute HX
	disp('Calculate HX');
	tic;
	HX = zeros(n,p);
	for j = 1:size(D,1)
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        sel = (m1*(j-1)+1):(j*m1);
        HX = HX + H*X(sel,:);
        clear H;
    end
	disp(toc);
    
    % Compute the LHS of the inversion equations
	disp('Compute the left hand side of the inverse modeling equations');
        tic;
	LHS = [[psi HX] ; [HX' zeros(p,p)]];
        save(strcat(outpath,'psi_HX.mat'),'psi','HX','-v7.3');
	clear psi HX;
	disp(toc);
    
    % Compute the inversion weights (see Eqs. 25 and 28 in Michalak et al. 2004i)
	disp('Compute the weights');
	tic;
	RHS = [ Z ; zeros(p,1)];
	weights = LHS \ RHS;
	disp(toc);
    
    eta   = weights(1:(size(weights,1)-p));
	betas = weights((size(weights,1)-p+1):size(weights,1));
    
    % Calculate H * eta
	disp('Calculate H^T * eta');
	tic;
    Heta = [];
	for j = 1:size(D,1)
        load(strcat(Hpath,'H_',num2str(j),'.mat'));
        Heta = [Heta; H'*eta];
    end
	disp(toc);
    
    % Calculate Q * H^T * eta
	disp('Calculate Q * H^T * eta');
	tic;
	QHeta = [];
	
	for j = 1:size(Dinv,1)
        A1 = zeros(m1,1);
        for i = 1:size(D,1)
            sel = (m1*(i-1)+1):(i*m1);
%             A1 = A1 + Heta(sel,1) .* D(j,i) .* sigmaQ(i) .* sigmaQ(j);	
              A1 = A1 + Heta(sel,1) .* D(j,i);	  
        end % End of i loop
        temp = E * A1;
        QHeta = [QHeta; temp];
    end % End of j loop
	clear A1 temp;
	disp(toc);
    
    % Create the estimate of the fluxes
	disp('Estimate the fluxes');
	tic;
	shat = X*betas + QHeta;
	disp(toc);
    
    % Write the outputs to file

    disp('Writing outputs to file');
    dlmwrite(outfile_fluxes,full(shat),',');

    disp('Outputs written to file');
    disp(outfile_fluxes)

