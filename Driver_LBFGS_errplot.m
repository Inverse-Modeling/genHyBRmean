% This script is for testing case study with hybrid methods
%
% Latest Update by T.Cho, Jun. 10, 2021

%%
% Add path

    addpath(genpath('genHyBR-master'));
    addpath('lbfgs');
    
% Fix reg. parameter

    fix_reg_param = 1;    
    
% There are two case studies

    CaseStudy = '6wk';
    % CaseStudy = '1yr';

% For each case, we tested with fire and without fire.

    Fire = 'on';
    % Fire = 'off';

% maximum iteration of L-BFGS
    
    maxit = 50; 
    
% Choose noise level: 5%, 10%, 50%, 100% (noise ratio to observation)
    
    NoiseOptions = 50;
    
% Setup path and filenames for saving outputs
    
    outpath = strcat('output/',CaseStudy,'/Fire_',Fire,'/lbfgs/');
    
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
    spath = strcat(CaseStudy,'/Fire_',Fire,'/');

% Create s_true
    s_true = [];
    for j = 1:ntimes
        load(strcat(spath,'s_',num2str(j),'.mat'));
        s_true = [s_true; s(:)]; 
    end
    
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
    
    m = size(X,1);
    switch CaseStudy
        
        case '6wk'
                sigmaQ = ones(size(D,1),1) / (fix_reg_param^2);
        case '1yr'
                D_theta = [ sqrt(34.7) sqrt(14.84) sqrt(16.09) sqrt(8.32) sqrt(6.66) sqrt(7.47) sqrt(9.67) sqrt(9.37) sqrt(39.07) sqrt(72.63) sqrt(102.37) sqrt(70.43) ]';
                sigmaQ  = [ D_theta(1) .* ones(8.*30,1); D_theta(2) .* ones(8.*31,1); D_theta(3) .* ones(8.*30,1); ...
                            D_theta(4) .* ones(8.*31,1); D_theta(5) .* ones(8.*31,1); D_theta(6) .* ones(8.*28,1); ...
                            D_theta(7) .* ones(8.*31,1); D_theta(8) .* ones(8.*30,1); D_theta(9) .* ones(8.*31,1); ...
                            D_theta(10) .* ones(8.*30,1); D_theta(11) .* ones(8.*31,1); D_theta(12) .* ones(8.*31,1)];
                % scaling sigmaQ
                sigmaQ = sigmaQ/max(sigmaQ(:));  
                sigmaQ = sigmaQ / (fix_reg_param^2);
    end
    
    % scaling D
    D = (sigmaQ*sigmaQ') .* D;
    
    % Take the inverse of D
	Dinv      = inv(D);
    
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
    
%% Esimate the fluxes with L-BFGS-B algorithm

    disp('Create an initial guess for the L-BFGS algorithm');
    shat0 = zeros(m,1);

    % Pre-calculate matrix products where possible 
    % Calculate inv(Q) * X
    % We'll refer to this matrix product as varible "B" from now on
	% This step is really slow and is the time-limiting step in the inverse model
    
    disp('Pre-calculate matrix products where possible');
	
	% B = inv(Q) * X
	B = [];
    
    for j = 1:ntimes
        disp(j);
        tic;
            B1 = zeros(m1,size(X,2));
                for i = 1:size(Dinv,1)
                sel = (m1.*(i-1)+1):(i.*m1);
%                 B1 = B1 + X(sel,:) .* Dinv(j,i) .* (1./sigmaQ(i)) .* (1./sigmaQ(j));
                B1 = B1 + X(sel,:) .* Dinv(j,i);
                end % End of i loop
            temp = Einv * B1;
            B = [B; temp];
        disp(toc);
    end % End of j loop
	clear B1 temp;

	save(strcat(outpath,'nlevel_',num2str(NoiseOptions),'_B.mat'),'B');
% 	load(strcat(outpath,'nlevel_',num2str(NoiseOptions),'_B.mat'));

%%
    % Set the time counter
        disp('Time at the initiation of the L-BFGS algorithm:');
        disp(clock);
        
    % For the cost function
    f1 = @(shat) cost_gradient_fun(Z, R, X, B, Dinv, Einv, Hpath, shat);

    % Create an empty flux estimate
    shat = [];    
    
    % Set options for the L-BFGS-B algorithm
    options = struct('HessUpdate','lbfgs','GradObj','on','Display','iter','MaxIter',maxit,'GradConstr',false);

    
    outname = strcat(outpath,'fluxes_LBFGS',num2str(maxit),'.csv');
    
    % Run L-BFGS-B
    disp('Run the L-BFGS-B algorithm');

    [shat,costfun,exitflag,gradient,s_errplot] = fminlbfgs_errplot(f1,shat0,options,s_true);

    % Print out time information
    disp('Time at the end of the L-BFGS algorithm:');
    disp(clock);

    % Write the outputs to file
    disp('Writing outputs to file');
    dlmwrite(outname,full(shat),',');
	save(strcat(outpath,'fluxes_LBFGS',num2str(maxit),'_nlevel_',num2str(NoiseOptions),'.mat'),'shat');
	
    disp('Outputs written to file');
    disp(outname);	
    
    convergence = s_errplot;
    save(strcat(outpath,'convergence_LBFGS_',num2str(maxit),'_nlevel_',num2str(NoiseOptions),'.mat'),'convergence');
    