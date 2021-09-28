%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_run												%
% PURPOSE: Launch a geostatistical inverse model (GIM), implemented with Lagrange multipliers for non-negativity	%
% S. Miller, Jan. 8, 2016												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%

%---------------------------%
% Summary of this script    %
%---------------------------%

	% L-BFGS algorithm
	% Here is the link to the L-BFGS algorithm that I've implemented in this code: 
	% https://www.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer
	% Note: there are several codes available that implement L-BFGS in Matlab. 


%----------------------------------------------------------%
% OPTIONAL: Loop over different numbers of iterations.     %
%----------------------------------------------------------%

        % There are two options here.
        % 1) I could construct a parfor loop and loop over each iteration
        % 2) I could launch a shell script that iterates through 1:50.
        % I might do the latter. 

        % Then I can test convergence by iteration number.
        % Submitting 50 jobs is not that big of a deal.

        maxit = str2num(getenv('maxit'));
        disp('maxit');
        disp(num2str(maxit));


%----------------------------%
% Initial setup parameters   %
%----------------------------%
	
	%-----------------------------------------%
	% **Set the name of the output folder**   %
	%-----------------------------------------%

	% This script will save the estimated fluxes into the following folder:
	outpath = '/home-4/smill191@jhu.edu/work/smiller/GIM_test/LBFGS_year/';
	
	
	%---------------------------------------------%
	% **Set the covariance matrix parameters**    %
	%---------------------------------------------%

        % Set the covariance matrix parameters for the inversion
        theta = [ 2.000 sqrt(34.7) sqrt(14.84) sqrt(16.09) sqrt(8.32) sqrt(6.66) sqrt(7.47) ...
                sqrt(9.67) sqrt(9.37) sqrt(39.07) sqrt(72.63) sqrt(102.37) sqrt(70.43)  585.68 12.366 ];	
	
	% Display the covariance matrix parameters on screen
	disp('Covariance matrix parameters');
	disp(theta);
	

%---------------------------------%
% **Read in the observations**    %
%---------------------------------%

	disp('Read in the observation vector');

%       ntimes = 2920; 
%        load(strcat(Hpath,'H_1.mat'));
%       n      = size(H,1);
%        Z      = zeros(n,1);
%        for j = 1:ntimes;
%       disp(num2str(j));
%        load(strcat(Hpath,'H_',num2str(j),'.mat'));
%        load(strcat(Hpath,'s_',num2str(j),'.mat'));
%       Z = Z + H*s;
%        clear H;
%        end;
%       Z = Z + normrnd(0,theta(1),[length(Z),1]);
%       save(strcat(outpath,'Z.mat'),'Z');

        load('/home-4/smill191@jhu.edu/work/smiller/GIM_test/minres_year/Z.mat');


%--------------------------------%
% Set the path to the H matrix   %
%--------------------------------%

        Hpath = '~/work/smiller/data/OCO2_H_matrix_year/';


%----------------------------%
% **Create the X matrix**    %
%----------------------------%

        disp('Create the X matrix');

        tic;
        n = length(Z);

%        %-------------------------------------------------%
%        % Create a temporary X object for a 31-day month  %
%        %-------------------------------------------------%
%
%        X = [] ;
%        m1 = 3222
%        for i = 1:8
%                for j = 1:(m1.*31.*8)
%                        if rem((fix((j-1)/3222) + 1 - i), 8) == 0
%                        X(j, i) = 1 ;
%                        else
%                        X(j, i) = 0 ;
%                        end
%                end
%        end
%        Xtemp31 = sparse(X);
%
%        %-------------------------------------------------%
%        % Create a temporary X object for a 30-day month  %
%        %-------------------------------------------------%
%
%        X = [] ;
%        for i = 1:8
%                for j = 1:(m1.*30.*8)
%                        if rem((fix((j-1)/3222) + 1 - i), 8) == 0
%                        X(j, i) = 1 ;
%                        else
%                        X(j, i) = 0 ;
%                        end
%                end
%        end
%        Xtemp30 = sparse(X);
%
%
%        %-------------------------------------------------%
%        % Create a temporary X object for a 28-day month  %
%        %-------------------------------------------------%
%
%        X = [] ;
%        for i = 1:8
%                for j = 1:(m1.*28.*8)
%                        if rem((fix((j-1)/3222) + 1 - i), 8) == 0
%                        X(j, i) = 1 ;
%                        else
%                        X(j, i) = 0 ;
%                        end
%                end
%        end
%        Xtemp28 = sparse(X);
%
%        %-------------------------------------------------%
%        % Append the temporary X matrices together
%        %-------------------------------------------------%
%
%        X = blkdiag(Xtemp30,Xtemp31,Xtemp30,Xtemp31,Xtemp31,Xtemp28,Xtemp31,Xtemp30,Xtemp31,Xtemp30,Xtemp31,Xtemp31);
%
%        disp(toc);
%
	
	%-------------------------------------------------%
	% Alternate construction of X with 12 columns     %
	%-------------------------------------------------%

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


%-------------------------%
% Create the E matrix     %
%-------------------------%

	disp('Create the E and D matrices');

	% The E and D matrices are used to construct the Q covariance matrix. 
	% The E matrix describes how covariances in Q decay in space, and the D matrix describes how covariances in Q decay in time.
	% See the following paper by Vineet Yadav for more details on the E and D matrices: http://www.geosci-model-dev.net/6/583/2013/gmd-6-583-2013.html
	% The E matrix has dimensions (r x r)
	% The D matrix has dimensions (q x q)
	
	
	%----------------------------------------------%
	% **Read in the geographic distance matrix**   %
	%----------------------------------------------%
	
	% This matrix (dimensions r x r) should define the distance between each emissions grid box
	% This matrix is symmetric with zeros on the diagonals. 
	% This matrix should typically has units of km. Whatever units you choose should match the units of theta(3)

        load('/home-4/smill191@jhu.edu/work/smiller/GIM_test/distmat.mat');
	
	
	%------------%
	% Create E   %
	%------------%

	% No need to edit this section
	% Note: I use a spherical covariance model because of its sparse properties.

	% Multiply the distance matrix by the covariance model
	% Exponential covariance model
	% decaymat = exp(-1 .* deltamat ./ theta(3));

        % Spherical covariance model
        E                      = 1 - 1.5 .* (deltamat ./theta(14))  + 0.5 .* (deltamat.^3 ./ theta(14).^3);
        E(deltamat > theta(14)) = 0;

        Einv = inv(E);
	
	
%----------------------------%
% Create the D matrix        %
%----------------------------%

	% The D matrix describes the time between each time period in the inversion
	% It will have dimensions (q x q)
	% The D matrix is usually easier to set up than the E matrix, and I have included sample code below
	
	%-----------------------------------%
	% **Create time distance matrix**   %
	%-----------------------------------%

        days = [];
        for i = 1:2920
        for j = 1:2920
            if rem(abs(i-j),8) == 0
                days(i,j) = abs(i-j)/8;
            else
                days(i,j) = 10^6;
            end
        end
        end

	%------------%
	% Create D   %
	%------------%

        % Spherical covariance model
        D = 1 - 1.5 .* (days   ./theta(15))  + 0.5 .* (days.^3   ./ theta(15).^3);
        D(days   > theta(15)) = 0;

        Dinv = inv(D);

        disp(toc);


%----------------------------%
% Create the sigmaQ vector   %
%----------------------------%

        m = size(X,1);

        % sigmaQ = [ theta(2) .* ones(m1.*8.*30,1); theta(3) .* ones(m1.*8.*31,1); theta(4) .* ones(m1.*8.*30,1); ...
        %         theta(5) .* ones(m1.*8.*31,1); theta(6) .* ones(m1.*8.*31,1); theta(7) .* ones(m1.*8.*28,1); ...
        %         theta(8) .* ones(m1.*8.*31,1); theta(9) .* ones(m1.*8.*30,1); theta(10) .* ones(m1.*8.*31,1); ...
        %         theta(11) .* ones(m1.*8.*30,1); theta(12) .* ones(m1.*8.*31,1); theta(13) .* ones(m1.*8.*31,1)];

        % Create a vector that has one value of sigmaQ for each unique time period in the inverse model
        sigmaQ = [ theta(2) .* ones(8.*30,1); theta(3) .* ones(8.*31,1); theta(4) .* ones(8.*30,1); ...
                theta(5) .* ones(8.*31,1); theta(6) .* ones(8.*31,1); theta(7) .* ones(8.*28,1); ...
                theta(8) .* ones(8.*31,1); theta(9) .* ones(8.*30,1); theta(10) .* ones(8.*31,1); ...
                theta(11) .* ones(8.*30,1); theta(12) .* ones(8.*31,1); theta(13) .* ones(8.*31,1)];


%------------------------%
% Create the R matrix    %
%------------------------%
		
	disp('Create R');
	
	% No need to edit this section, unless you want a more complicated setup for R.
	R = (theta(1).^2) .* ones(n,1);
	R = spdiags(R,0,n,n);


%------------------------------------------------------%
% Create the initial guess for the L-BFGS-B algorithm  %
%------------------------------------------------------%

	disp('Create an initial guess for the L-BFGS algorithm');

	% For now, I'll use a really simple initial guess. I want to test
	% how the algorithm converges without making any assumptions about the 
	% intelligence of our initial guess.

	shat0 = zeros(m,1);


%-------------------------------------------------%
% Pre-calculate matrix products where possible    %
%-------------------------------------------------%

	% Calculate inv(Q) * X
	% We'll refer to this matrix product as varible "B" from now on
	% This step is really slow and is the time-limiting step in the inverse model

	disp('Pre-calculate matrix products where possible');
	
	% B = inv(Q) * X
	B = [];
	p = size(X,2);
	m1  = size(E,1);
	ntimes = size(D,1);

%	for j = 1:ntimes;
% disp(j);
% tic;
%	B1 = zeros(m1,size(X,2));
%		for i = 1:size(Dinv,1);
%		sel = (m1.*(i-1)+1):(i.*m1);
%		B1 = B1 + X(sel,:) .* Dinv(j,i) .* (1./sigmaQ(i)) .* (1./sigmaQ(j));
%		end; % End of i loop
%	temp = Einv * B1;
%	B = [B; temp];
% disp(toc);
%	end; % End of j loop
%	clear B1 temp;
%
%	save(strcat(outpath,'B.mat'),'B');
	load(strcat(outpath,'B.mat'));	
	
%--------------------------------------------------%
% Estimate the fluxes using Lagrange multipliers   %
%--------------------------------------------------%

        % Set the time counter
        disp('Time at the initiation of the L-BFGS algorithm:');
        disp(clock);


%-------------------------------------------------------------%
% Create function handles for the cost function and gradient  %
%-------------------------------------------------------------%

        % For the cost function
        f1 = @(shat) cost_gradient_fun(Z, R, X, B, Dinv, Einv, sigmaQ, Hpath, shat);

        % Create an empty flux estimate
        shat = [];


%-----------------------------------------%
% Set options for the L-BFGS-B algorithm  %
%-----------------------------------------%

        options = struct('HessUpdate','lbfgs','GradObj','on','Display','iter','MaxIter',maxit,'GradConstr',false);


%-------------------------------------------------%
% Check to see if the solution already exists.    %
%-------------------------------------------------%

	% If so, read in the fluxes from file
        outname = strcat(outpath,'fluxes_LBFGS',num2str(maxit),'.csv');

	fflag = exist(outname);

	if fflag == 2;
	disp('Fluxes found. Will read in existing file.');
	shat = csvread(outname);
	
	% Calculate the cost function
	disp('Calculate cost function');
	[costfun gradfun] = cost_gradient_fun(Z, R, X, B, Dinv, Einv, sigmaQ, Hpath, shat);	

	end;

	if fflag == 0;

%-----------------------%
% Run the algorithm     %
%-----------------------%

        disp('Run the L-BFGS-B algorithm');

        [shat,costfun,exitflag,gradient] = fminlbfgs(f1,shat0,options);


        % Print out time information
        disp('Time at the end of the L-BFGS algorithm:');
        disp(clock);
	
	
%------------------------------%
% Write the outputs to file    %
%------------------------------%

        disp('Writing outputs to file');
        dlmwrite(outname,full(shat),',');
	save(strcat(outpath,'fluxes_LBFGS',num2str(maxit),'.mat'),'shat');
	
        disp('Outputs written to file');
        disp(outname);	

	end; % End of fflag if statement


%-------------------------------------------------%
% Calculate RMSE compared to analytical solution  %
%-------------------------------------------------%

        % Calculate the root mean squared error
        % rmse = sqrt(mean((shat - shata).^2));

        % disp('Current root mean squared error');
        % disp('(relative to analytical solution)');
        % disp(num2str(rmse));

        % Save these metrics to a file
        % load('/home-4/smill191@jhu.edu/work/smiller/GIM_test/convergence_LBFGS.mat');
        convergence = [maxit costfun];
        % convergence = [convergence; temp];
        save(strcat(outpath,'convergence_LBFGS_',num2str(maxit),'.mat'),'convergence');


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
