% This script is for plotting relative errors of case study with various
% noise levels 
%
% Required package:
%   - panel (https://www.mathworks.com/matlabcentral/fileexchange/20003-panel)
%
% Latest Update by T.Cho, Nov. 28, 2021

%% Load data

    CaseStudy = '6wk';

    ntimes = 328;
    Hpath = strcat(CaseStudy,'/H/');
    load(strcat(Hpath,'distmat.mat'));
    
    theta = [ 2.000 10.046 555.420 9.854 ];
    disp('Covariance matrix parameters');
    disp(theta);

    % Create E
    % Spherical covariance model
    E_6wk                      = 1 - 1.5 .* (deltamat ./theta(3))  + 0.5 .* (deltamat.^3 ./ theta(3).^3);
    E_6wk(deltamat > theta(3)) = 0;

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
    D_6wk = 1 - 1.5 .* (days   ./theta(4))  + 0.5 .* (days.^3   ./ theta(4).^3);
    D_6wk(days   > theta(4)) = 0;


%% Plot 6 week case study
    n = 40;
    figure, 
    p = panel();
    p.pack(1,2);

    p(1,1).select(), cspy(D_6wk,'colormap','parula','levels',n,'caxis',linspace(0,1,n)), title('6wk Temporal')
    set(gca,'XAxisLocation','bottom'), xlabel('')
    p(1,2).select(), cspy(E_6wk,'colormap','parula','levels',n,'caxis',linspace(0,1,n)), title('6wk Spatial')
    set(gca,'XAxisLocation','bottom'), xlabel('')

    set(gcf,'Position',[300 300 800 375])
    p.fontsize = 15;
    p.margin = [10 15 15 10];