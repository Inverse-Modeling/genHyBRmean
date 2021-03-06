% This script is for plotting reconstructions of case study with various
% noise levels
%
% Required packages:
%   - panel (https://www.mathworks.com/matlabcentral/fileexchange/20003-panel)
%   - cmocean (https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)
%   - cspy_TC.m is modified version of cspy.m from (https://www.mathworks.com/matlabcentral/fileexchange/46551-cspy-m)
%     in order to add caxis option and add cmocean color.
%
% Latest Update by T.Cho, Nov. 28, 2021

%% Load grid and true average image
    % Add path
    addpath(genpath('functions'));
    
    CaseStudy = '6wk';

    Hpath = strcat(CaseStudy,'/H/');
    load(strcat(Hpath,'land_mask.mat'));
    lon = -179.5:1:-10.5;
    lat = 10.5:1:79.5;

    % Get True average image
    S_true_avg = GetAvgTrueImage(CaseStudy);

%% Plot
    
    % Colorbar Scale (y/n)
    ScalefromTrue = 'y';
    
    figure,
    p=panel();
    p.pack(2,4);
    caxis_scale = [-5, 2];

        % Plot True
            p(1,1).select()
            Tit = strcat(CaseStudy,' - True Average');
            % Plot image 
            Plot_CO2_Recons(S_true_avg,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(S_true_avg,20)       
            colorbar off;
        
        % Direct 
            M = "Load estimated fluxes of 6wk case study using direct method with 0.0653 regularization parameter for 50% noise"
            Shat = GetAvgReconImage(M,CaseStudy,'direct');
            S_relerr = norm(S_true_avg(:)-Shat(:))/norm(S_true_avg(:));
            switch ScalefromTrue
                case 'n'
                    caxis_scale = [min(Shat(:)), max(Shat(:))];
            end
            p(2,1).select()
            Tit = strcat('Direct fixed \lambda = 0.0653, (',num2str(S_relerr),')');
            % Plot image 
            Plot_CO2_Recons(Shat,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(Shat,20)
            colorbar off    

        % genhybr
            % Plot Optimal Reg.
            M = "Load estimated fluxes of 6wk case study using genHyBRs with optimal regularization parameter for 50% noise"
            Shat = GetAvgReconImage(M,CaseStudy,'genhybr');
            S_relerr = norm(S_true_avg(:)-Shat(:))/norm(S_true_avg(:));
            switch ScalefromTrue
                case 'n'
                    caxis_scale = [min(Shat(:)), max(Shat(:))];
            end
            p(1,2).select()
            Tit = strcat('genHyBRs','-opt (',num2str(S_relerr),')');
            % Plot image 
            Plot_CO2_Recons(Shat,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(Shat,20)
            colorbar off;

            % Plot DP Reg.
            M = "Load estimated fluxes of 6wk case study using genHyBRs with dp for 50% noise"
            Shat = GetAvgReconImage(M,CaseStudy,'genhybr');
            S_relerr = norm(S_true_avg(:)-Shat(:))/norm(S_true_avg(:));
            switch ScalefromTrue
                case 'n'
                    caxis_scale = [min(Shat(:)), max(Shat(:))];
            end
            p(1,3).select()
            Tit = strcat('genHyBRs','-dp (',num2str(S_relerr),')');
            % Plot image 
            Plot_CO2_Recons(Shat,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(Shat,20)
            colorbar off;
        
            % Plot no Reg.
            M = "Load estimated fluxes of 6wk case study using genHyBRs with no regularization parameter for 50% noise"
            Shat = GetAvgReconImage(M,CaseStudy,'genhybr');
            S_relerr = norm(S_true_avg(:)-Shat(:))/norm(S_true_avg(:));
            switch ScalefromTrue
                case 'n'
                    caxis_scale = [min(Shat(:)), max(Shat(:))];
            end
            p(1,4).select()
            Tit = strcat('genHyBRs','-none (',num2str(S_relerr),')');
            % Plot image 
            Plot_CO2_Recons(Shat,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(Shat,20)
            colorbar off;

        % meanest
            % Plot Optimal Reg.
            M = "Load estimated fluxes of 6wk case study using genHyBRmean with optimal regularization parameter for 50% noise"
            Shat = GetAvgReconImage(M,CaseStudy,'meanest');
            S_relerr = norm(S_true_avg(:)-Shat(:))/norm(S_true_avg(:));
            switch ScalefromTrue
                case 'n'
                    caxis_scale = [min(Shat(:)), max(Shat(:))];
            end
            p(2,2).select()
            Tit = strcat('genHyBRmean','-opt (',num2str(S_relerr),')');
            % Plot image 
            Plot_CO2_Recons(Shat,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(Shat,20)
            colorbar off;

            % Plot DP Reg.
            M = "Load estimated fluxes of 6wk case study using genHyBRmean with dp for 50% noise"
            Shat = GetAvgReconImage(M,CaseStudy,'meanest');
            S_relerr = norm(S_true_avg(:)-Shat(:))/norm(S_true_avg(:));
            switch ScalefromTrue
                case 'n'
                    caxis_scale = [min(Shat(:)), max(Shat(:))];
            end
            p(2,3).select()
            Tit = strcat('genHyBRmean','-dp (',num2str(S_relerr),')');
            % Plot image 
            Plot_CO2_Recons(Shat,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(Shat,20)
            colorbar off;

            % Plot no Reg.
            M = "Load estimated fluxes of 6wk case study using genHyBRmean with no regularization parameter for 50% noise"
            Shat = GetAvgReconImage(M,CaseStudy,'meanest');
            S_relerr = norm(S_true_avg(:)-Shat(:))/norm(S_true_avg(:));
            switch ScalefromTrue
                case 'n'
                    caxis_scale = [min(Shat(:)), max(Shat(:))];
            end
            p(2,4).select()
            Tit = strcat('genHyBRmean','-none (',num2str(S_relerr),')');
            % Plot image 
            Plot_CO2_Recons(Shat,Tit,caxis_scale,land_mask,lat,lon)
            % Remove Africa from image
            Plot_noAfrica(Shat,20)
            
    h = colorbar;
    xlabel(h, '\mu mol m^{-2} s^{-1}')        
    set(gcf,'Position',[200 200 1600 650])
    p.fontsize = 15;
    p.margin = [15 20 30 20];

