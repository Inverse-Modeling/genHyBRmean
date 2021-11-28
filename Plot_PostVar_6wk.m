% This script is for plotting the average posterior standard deviation over
% 6wk for both genHyBRs and genHyBRmean
%
% Required packages:
%   - panel (https://www.mathworks.com/matlabcentral/fileexchange/20003-panel)
%   - cmocean (https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)
%   - cspy_TC.m is modified version of cspy.m from (https://www.mathworks.com/matlabcentral/fileexchange/46551-cspy-m)
%     in order to add caxis option and add cmocean color.
%
% Latest Update by T.Cho, Nov. 28, 2021

%% Load data

    CaseStudy = '6wk';

    spath = strcat(CaseStudy,'/s/');
    Hpath = strcat(CaseStudy,'/H/');
    load(strcat(Hpath,'land_mask'));
    load(strcat(Hpath,'distmat.mat'));
    lon = -179.5:1:-10.5;
    lat = 10.5:1:79.5;

    % load observation vector and noise s.d (sigma)
    Zmat = load(strcat(spath,'Z_',num2str(50),'.mat'));
    Z = Zmat.Zn;

    % Set the covariance matrix parameters
    theta = [ 2.000 10.046 555.420 9.854 ];
    theta(1) = Zmat.sigma;
    sigma = theta(1);
    n = length(Z);

%% Compute approximation of posterior s.d. from genhybrs
    HybridOptions = 'genhybr';
    output_6wk_genhybr_opt = "Load Output of 6wk case study using genHyBRs with optimal regularization parameter for 50% noise"
    output_6wk_genhybr_dp = "Load Output of 6wk case study using genHyBRs with dp for 50% noise"

    % Create Prior Q
    Q_genhybr = Create_Q(CaseStudy,HybridOptions,deltamat,theta,n);

    PostVar  = estdiagpost(output_6wk_genhybr_opt.output, Q_genhybr, sigma^2, HybridOptions);
    PostSD = sqrt(PostVar);
    SD_genhybr_opt = GetAvgReconImage(PostSD, CaseStudy, HybridOptions);

    PostVar  = estdiagpost(output_6wk_genhybr_dp.output, Q_genhybr, sigma^2, HybridOptions);
    PostSD = sqrt(PostVar);
    SD_genhybr_DP = GetAvgReconImage(PostSD, CaseStudy, HybridOptions);

%% Compute approximation of posterior s.d. from genhybrmean
    HybridOptions = 'meanest';
    output_6wk_meanest_opt = "Load Output of 6wk case study using genHyBRmean with optimal regularization parameter for 50% noise"
    output_6wk_meanest_dp = "Load Output of 6wk case study using genHyBRmean with dp for 50% noise"

    % Create Prior Q
    Q_meanest = Create_Q(CaseStudy,HybridOptions,deltamat,theta,n);

    PostVar  = estdiagpost(output_6wk_meanest_opt.output, Q_meanest, sigma^2, HybridOptions);
    PostSD = sqrt(PostVar);
    SD_meanest_opt = GetAvgReconImage(PostSD, CaseStudy, HybridOptions);

    PostVar  = estdiagpost(output_6wk_meanest_dp.output, Q_meanest, sigma^2, HybridOptions);
    PostSD = sqrt(PostVar);
    SD_meanest_DP = GetAvgReconImage(PostSD, CaseStudy, HybridOptions);

%% Plot image of post s.d.

    % Only show for dp
    figure,
        p=panel();
        p.pack(1,2);

        p(1,1).select()
        Tit = 'genHyBRs-dp (average)';
        Plot_CO2_Recons(SD_genhybr_DP, Tit, [min(SD_genhybr_DP(find(SD_genhybr_DP(:)))), max(SD_genhybr_DP(:))], land_mask, lat, lon)
        Plot_noAfrica(SD_genhybr_DP,20)
        h = colorbar;
        h.Location = 'east';
        h.AxisLocation = 'out';
        xlabel(h, '\mu mol m^{-2} s^{-1}')

        p(1,2).select()
        Tit = 'genHyBRmean-dp (average)';
        Plot_CO2_Recons(SD_meanest_DP, Tit, [min(SD_meanest_DP(find(SD_meanest_DP(:)))), max(SD_meanest_DP(:))], land_mask, lat, lon)
        Plot_noAfrica(SD_meanest_DP,20)
        h = colorbar;
        h.Location = 'east';
        h.AxisLocation = 'out';
        xlabel(h, '\mu mol m^{-2} s^{-1}')

    set(gcf,'Position',[200 300 800 350])
    p.fontsize = 15;
    p.margin = [15 20 30 20];
