function S_fin = GetAvgReconImage(M,CaseStudy,HybridOptions,TakeOffTimes)
% This function is to get average reconstructions from estimated fluxes
% 
%   M             - estimated fluxes using HybridOptions
%   CaseStudy     - {'6wk', '1yr'}
%   HybridOptions - {'genhybr','meanest','direct'}
%   TakeOffTimes  - early time period not to be selected on average image
%
% Latest Update by T.Cho, Nov. 28, 2021   

if nargin < 4
        TakeOffTimes = 0;
    end

    switch CaseStudy
        case '6wk'
            ntimes = 328;
            switch HybridOptions
                case 'genhybr'
                    nbeta = 0;
                case 'meanest'
                    nbeta = 8;
                case 'direct'
                    nbeta = 0;
            end
        case '1yr'
            ntimes = 2920;
            switch HybridOptions
                case 'genhybr'
                    nbeta = 0;
                case 'meanest'
                    nbeta = 12;
                case 'direct'
                    nbeta = 0;    
            end
    end
    
    M = M(1:end-nbeta);

    S = reshape(M,3222,ntimes);
    S_takeofftimes = S(:,TakeOffTimes+1:end);
    S_sum = sum(S_takeofftimes,2);
    S_avg = S_sum / ntimes;

    lon = -179.5:1:-10.5;
    lat = 10.5:1:79.5;
    Hpath = strcat(CaseStudy,'/H/');
    load(strcat(Hpath,'land_mask.mat'));

    S_grid = zeros(length(lon),length(lat));
    S_grid(land_mask) = S_avg;
    S_fin = flipud(S_grid');

end