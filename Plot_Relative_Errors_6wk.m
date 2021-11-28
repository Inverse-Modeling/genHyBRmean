% This script is for plotting relative errors of case study with various
% noise levels 
%
% Required package:
%   - panel (https://www.mathworks.com/matlabcentral/fileexchange/20003-panel)
%
% Latest Update by T.Cho, Nov. 28, 2021

%% Load output
% Load True to compute relative errors of direct method
    S_true = [];
    spath = strcat('6wk/s/');
    for i = 1:ntimes
        load(strcat(spath,'s_',num2str(i),'.mat'));
        S_true = [S_true; s(:)];
    end
    
    % 5% added noise
        % genHyBRs
        output_6wk_nlevel_5_genhybr_opt = "Load output of 6wk case study using genHyBRs with optimal reguralization parameter for 5% noise"
        output_6wk_nlevel_5_genhybr_dp = "Load output of 6wk case study using genHyBRs with dp for 5% noise"
        output_6wk_nlevel_5_genhybr_lsqr = "Load output of 6wk case study using genHyBRs with no regularization parameter for 5% noise" 

        % genHyBRmean
        output_6wk_nlevel_5_meanest_opt = "Load output of 6wk case study using genHyBRmean with optimal reguralization parameter for 5% noise"
        output_6wk_nlevel_5_meanest_dp = "Load output of 6wk case study using genHyBRmean with dp for 5% noise"
        output_6wk_nlevel_5_meanest_lsqr = "Load output of 6wk case study using genHyBRmean with no regularization parameter for 5% noise" 

        % direct method
        fluxes_6wk_nlevel_5_direct_reg_1_0000 = "Load output of 6wk case study using direct method with 1 as regularization parameter for 5% noise" 
        fluxes_6wk_nlevel_5_direct_reg_0_0084 = "Load output of 6wk case study using direct method with 0.0084 as regularization parameter for 5% noise" 
        dir_relerr_nlevel_5_reg_1_0000 = norm(S_true-fluxes_6wk_nlevel_5_direct_reg_1_0000)/norm(S_true);
        dir_relerr_nlevel_5_reg_0_0084 = norm(S_true-fluxes_6wk_nlevel_5_direct_reg_0_0084)/norm(S_true);
        
    % 10% added noise
        % genHyBRs
        output_6wk_nlevel_10_genhybr_opt = "Load output of 6wk case study using genHyBRs with optimal reguralization parameter for 10% noise"
        output_6wk_nlevel_10_genhybr_dp = "Load output of 6wk case study using genHyBRs with dp for 10% noise"
        output_6wk_nlevel_10_genhybr_lsqr = "Load output of 6wk case study using genHyBRs with no regularization parameter for 10% noise" 

        % genHyBRmean
        output_6wk_nlevel_10_meanest_opt = "Load output of 6wk case study using genHyBRmean with optimal reguralization parameter for 10% noise"
        output_6wk_nlevel_10_meanest_dp = "Load output of 6wk case study using genHyBRmean with dp for 10% noise"
        output_6wk_nlevel_10_meanest_lsqr = "Load output of 6wk case study using genHyBRmean with no regularization parameter for 10% noise" 

        % Direct method
        fluxes_6wk_nlevel_10_direct_reg_1_0000 = "Load output of 6wk case study using direct method with 1 as regularization parameter for 10% noise" 
        fluxes_6wk_nlevel_10_direct_reg_0_0164 = "Load output of 6wk case study using direct method with 0.0164 as regularization parameter for 10% noise" 
        dir_relerr_nlevel_10_reg_1_0000 = norm(S_true-fluxes_6wk_nlevel_10_direct_reg_1_0000)/norm(S_true);
        dir_relerr_nlevel_10_reg_0_0164 = norm(S_true-fluxes_6wk_nlevel_10_direct_reg_0_0164)/norm(S_true);    

    % 50% added noise
        % genHyBRs
        output_6wk_nlevel_50_genhybr_opt = "Load output of 6wk case study using genHyBRs with optimal reguralization parameter for 50% noise"
        output_6wk_nlevel_50_genhybr_dp = "Load output of 6wk case study using genHyBRs with dp for 50% noise"
        output_6wk_nlevel_50_genhybr_lsqr = "Load output of 6wk case study using genHyBRs with no regularization parameter for 50% noise" 

        % genHyBRmean
        output_6wk_nlevel_50_meanest_opt = "Load output of 6wk case study using genHyBRmean with optimal reguralization parameter for 50% noise"
        output_6wk_nlevel_50_meanest_dp = "Load output of 6wk case study using genHyBRmean with dp for 50% noise"
        output_6wk_nlevel_50_meanest_lsqr = "Load output of 6wk case study using genHyBRmean with no regularization parameter for 50% noise" 

        % Direct method
        fluxes_6wk_nlevel_50_direct_reg_1_0000 = "Load output of 6wk case study using direct method with 1 as regularization parameter for 50% noise" 
        fluxes_6wk_nlevel_50_direct_reg_0_0653 = "Load output of 6wk case study using direct method with 0.0653 as regularization parameter for 50% noise" 
        dir_relerr_nlevel_50_reg_1_0000 = norm(S_true-fluxes_6wk_nlevel_50_direct_reg_1_0000)/norm(S_true);
        dir_relerr_nlevel_50_reg_0_0653 = norm(S_true-fluxes_6wk_nlevel_50_direct_reg_0_0653)/norm(S_true);    

%% Plot figure        
    % Define color
    color1 = [0 0.4470 0.7410];
    color2 = [0.8500 0.3250 0.0980];
    color3 = [0.4660, 0.6740, 0.1880];

    figure, 
    p=panel();
    p.pack(1,3); 

    % 5%
        nLevel = '5'; iter = 300;
        p(1,1).select(), hold on

        % Direct method
        p9 = plot(linspace(1,iter+1,100), dir_relerr_nlevel_5_reg_1_0000*ones(1,100),':k','LineWidth',1.5);
        text(iter+1,dir_relerr_nlevel_5_reg_1_0000,'\lambda=1','fontsize',12)
        p10 = plot(linspace(1,iter+1,100), dir_relerr_nlevel_5_reg_0_0084*ones(1,100),':k','LineWidth',1.5);
        text(iter+1,dir_relerr_nlevel_5_reg_0_0084,'\lambda=0.0084','fontsize',12)

        % genHyBRs
        p1 = plot(output_6wk_nlevel_5_genhybr_opt.output.Enrm,'Color',color1,'Linewidth',1.5);
        p2 = plot(output_6wk_nlevel_5_genhybr_dp.output.Enrm,'--','Color',color1,'Linewidth',1.5);
        p3 = plot(output_6wk_nlevel_5_genhybr_dp.output.iterations,...
            output_6wk_nlevel_5_genhybr_dp.output.Enrm(output_6wk_nlevel_5_genhybr_dp.output.iterations),...
            'o','Color',color1,'MarkerFaceColor',color1);
        p4 = plot(output_6wk_nlevel_5_genhybr_lsqr.output.Enrm,':','Color',color1,'LineWidth',1.5);

        % genHyBRmean
        p5 = plot(output_6wk_nlevel_5_meanest_opt.output.Enrm,'Color',color2,'Linewidth',1.5);
        p6 = plot(output_6wk_nlevel_5_meanest_dp.output.Enrm,'--','Color',color2,'Linewidth',1.5);
        p7 = plot(output_6wk_nlevel_5_meanest_dp.output.iterations,...
            output_6wk_nlevel_5_meanest_dp.output.Enrm(output_6wk_nlevel_5_meanest_dp.output.iterations),...
            'o','Color',color2,'MarkerFaceColor',color2);
        p8 = plot(output_6wk_nlevel_5_meanest_lsqr.output.Enrm,':','Color',color2,'LineWidth',1.5);

        hold off
        xlim([1 iter+1]); ylim([.6,1.01]);
        title(['nLevel (', num2str(nLevel),'%)'])
        xlabel('iterations')
        ylabel('relative errors')
        set(gca,'FontSize',15)

    % 10%
        nLevel = '10'; iter = 150;
        p(1,2).select(), hold on

        % Direct method
        p9 = plot(linspace(1,iter+1,100), dir_relerr_nlevel_10_reg_1_0000*ones(1,100),':k','LineWidth',1.5);
        text(iter+1,dir_relerr_nlevel_10_reg_1_0000,'\lambda=1','fontsize',12)
        p10 = plot(linspace(1,iter+1,100), dir_relerr_nlevel_10_reg_0_0164*ones(1,100),':k','LineWidth',1.5);
        text(iter+1,dir_relerr_nlevel_10_reg_0_0164,'\lambda=0.0164','fontsize',12)

        % genHyBRs
        p1 = plot(output_6wk_nlevel_10_genhybr_opt.output.Enrm,'Color',color1,'Linewidth',1.5);
        p2 = plot(output_6wk_nlevel_10_genhybr_dp.output.Enrm,'--','Color',color1,'Linewidth',1.5);
        p3 = plot(output_6wk_nlevel_10_genhybr_dp.output.iterations,...
            output_6wk_nlevel_10_genhybr_dp.output.Enrm(output_6wk_nlevel_10_genhybr_dp.output.iterations),...
            'o','Color',color1,'MarkerFaceColor',color1);
        p4 = plot(output_6wk_nlevel_10_genhybr_lsqr.output.Enrm,':','Color',color1,'LineWidth',1.5);

        % genHyBRmean
        p5 = plot(output_6wk_nlevel_10_meanest_opt.output.Enrm,'Color',color2,'Linewidth',1.5);
        p6 = plot(output_6wk_nlevel_10_meanest_dp.output.Enrm,'--','Color',color2,'Linewidth',1.5);
        p7 = plot(output_6wk_nlevel_10_meanest_dp.output.iterations,...
            output_6wk_nlevel_10_meanest_dp.output.Enrm(output_6wk_nlevel_10_meanest_dp.output.iterations),...
            'o','Color',color2,'MarkerFaceColor',color2);
        p8 = plot(output_6wk_nlevel_10_meanest_lsqr.output.Enrm,':','Color',color2,'LineWidth',1.5);

        hold off
        xlim([1 iter+1]); ylim([.6,1.01]);
        title(['nLevel (', num2str(nLevel),'%)'])
        xlabel('iterations')
        legend([p1,p2,p4,p5,p6,p8,p9],'genHyBRs - opt','genHyBRs - dp','genHyBRs - none', 'genHyBRmean - opt',...
                                          'genHyBRmean - dp','genHyBRmean - none','direct method (fixed \lambda)',...
                                          'location','southoutside','Orientation','horizontal','NumColumns',3)
        set(gca,'FontSize',15)

    % 50%
        nLevel = '50'; iter = 50;
        p(1,3).select(), hold on

        % Direct method
        p9 = plot(linspace(1,iter+1,100), dir_relerr_nlevel_50_reg_1_0000*ones(1,100),':k','LineWidth',1.5);
        text(iter+1,dir_relerr_nlevel_50_reg_1_0000,'\lambda=1','fontsize',12)
        p10 = plot(linspace(1,iter+1,100), dir_relerr_nlevel_50_reg_0_0653*ones(1,100),':k','LineWidth',1.5);
        text(iter+1,dir_relerr_nlevel_50_reg_0_0653,'\lambda=0.0653','fontsize',12)

        % genHyBRs
        p1 = plot(output_6wk_nlevel_50_genhybr_opt.output.Enrm,'Color',color1,'Linewidth',1.5);
        p2 = plot(output_6wk_nlevel_50_genhybr_dp.output.Enrm,'--','Color',color1,'Linewidth',1.5);
        p3 = plot(output_6wk_nlevel_50_genhybr_dp.output.iterations,...
            output_6wk_nlevel_50_genhybr_dp.output.Enrm(output_6wk_nlevel_50_genhybr_dp.output.iterations),...
            'o','Color',color1,'MarkerFaceColor',color1);
        p4 = plot(output_6wk_nlevel_50_genhybr_lsqr.output.Enrm,':','Color',color1,'LineWidth',1.5);


        % genHyBRmean
        p5 = plot(output_6wk_nlevel_50_meanest_opt.output.Enrm,'Color',color2,'Linewidth',1.5);
        p6 = plot(output_6wk_nlevel_50_meanest_dp.output.Enrm,'--','Color',color2,'Linewidth',1.5);
        p7 = plot(output_6wk_nlevel_50_meanest_dp.output.iterations,...
            output_6wk_nlevel_50_meanest_dp.output.Enrm(output_6wk_nlevel_50_meanest_dp.output.iterations),...
            'o','Color',color2,'MarkerFaceColor',color2);
        p8 = plot(output_6wk_nlevel_50_meanest_lsqr.output.Enrm,':','Color',color2,'LineWidth',1.5);

        hold off
        xlim([1 iter+1]); ylim([.6,1.01]);

        title(['nLevel (', num2str(nLevel),'%)'])
        xlabel('iterations')
        set(gca,'FontSize',15)

    set(gcf,'Position',[100 100 800 400])
    p.fontsize = 20;
    p.margin = [25 50 30 10];