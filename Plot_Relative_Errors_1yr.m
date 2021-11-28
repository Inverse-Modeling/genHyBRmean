% This script is for plotting relative errors of case study with various
% noise levels 
%
% Required package:
%   - panel (https://www.mathworks.com/matlabcentral/fileexchange/20003-panel)
%
% Latest Update by T.Cho, Nov. 28, 2021
%% Load output
   % 50% added noise 
        % genHyBRs
        output_1yr_nlevel_50_genhybr_opt = "Load output of 1yr case study using genHyBRs with optimal reguralization parameter for 50% noise"
        output_1yr_nlevel_50_genhybr_dp = "Load output of 1yr case study using genHyBRs with dpfor 50% noise"
        output_1yr_nlevel_50_genhybr_lsqr = "Load output of 1yr case study using genHyBRs with no reguralization parameter for 50% noise"

        % genHyBRmean
        output_1yr_nlevel_50_meanest_opt = "Load output of 1yr case study using genHyBRmean with optimal reguralization parameter for 50% noise"
        output_1yr_nlevel_50_meanest_dp = "Load output of 1yr case study using genHyBRmean with dp for 50% noise"
        output_1yr_nlevel_50_meanest_lsqr = "Load output of 1yr case study using genHyBRmean with no regularization parameter for 50% noise" 

%% Plot figure    
    % Define color
    color1 = [0 0.4470 0.7410];
    color2 = [0.8500 0.3250 0.0980];
    color3 = [0.4660, 0.6740, 0.1880];

    figure, 
    p=panel();
    p.pack(1,1)

    % 50%
        nLevel = '50'; iter = 100;
        p(1,1).select(), hold on
    
        % genHyBRs
        p1 = plot(output_1yr_nlevel_50_genhybr_opt.output.Enrm,'Color',color1,'Linewidth',1.5);
        p2 = plot(output_1yr_nlevel_50_genhybr_dp.output.Enrm,'--','Color',color1,'Linewidth',1.5);
        p3 = plot(output_1yr_nlevel_50_genhybr_dp.output.iterations,...
            output_1yr_nlevel_50_genhybr_dp.output.Enrm(output_1yr_nlevel_50_genhybr_dp.output.iterations),...
            'o','Color',color1,'MarkerFaceColor',color1);
        p4 = plot(output_1yr_nlevel_50_genhybr_lsqr.output.Enrm,':','Color',color1,'LineWidth',1.5);

        % genHyBRmean
        p5 = plot(output_1yr_nlevel_50_meanest_opt.output.Enrm,'Color',color2,'Linewidth',1.5);
        p6 = plot(output_1yr_nlevel_50_meanest_dp.output.Enrm,'--','Color',color2,'Linewidth',1.5);
        p7 = plot(output_1yr_nlevel_50_meanest_dp.output.iterations,...
            output_1yr_nlevel_50_meanest_dp.output.Enrm(output_1yr_nlevel_50_meanest_dp.output.iterations),...
            'o','Color',color2,'MarkerFaceColor',color2);
        p8 = plot(output_1yr_nlevel_50_meanest_lsqr.output.Enrm,':','Color',color2,'LineWidth',1.5);

        hold off
        xlim([1 iter+1]); ylim([output_1yr_nlevel_50_meanest_opt.output.Enrm(iter)*.95 1]);
        legend([p1,p2,p4,p5,p6,p8],'genHyBRs - opt','genHyBRs - dp','genHyBRs - none', 'genHyBRmean - opt',...
                                  'genHyBRmean - dp','genHyBRmean - none',...
                                  'location','southoutside','Orientation','horizontal','NumColumns',3)
        title(['nLevel (', num2str(nLevel),'%)'])
        xlabel('iterations')
        ylabel('relative errors')
        set(gca,'FontSize',15)

p.fontsize = 15;
p.margin = [25 50 30 10];