function Q = Create_Q(CaseStudy,HybridOptions,deltamat,theta,n)
% This function is to create covariance matrix Q
% 
%   CaseStudy     - {'6wk', '1yr'}
%   HybridOptions - {'genhybr','meanest','direct'}
%   deltamat      - matrix that lists the distance (in kilometers) 
%                   between the center of each model grid box
%   theta         - covariance matrix parameters
%   n             - size of observation
%
% Latest Update by T.Cho, Nov. 28, 2021   

    switch CaseStudy
        case '6wk'
            ntimes =328;
        case '1yr'
            ntimes =2920;
    end

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

    % Create X matrix
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


    switch HybridOptions
        case 'genhybr'
            Q = kronMat(D,E);

        case 'meanest'
            nbeta = size(X,2);
            gamma = 10;
            Qbeta = (1/gamma^2) * eye(nbeta);
            Q = kronMat(D,E);
            Q = QtilMat(Qbeta,Q,X);
    end