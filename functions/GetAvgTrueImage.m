function S_fin = GetAvgTrueImage(CaseStudy,TakeOffTimes)
% This function is to get average true image
% 
%   CaseStudy    - {'6wk', '1yr'}
%   TakeOffTimes - early time period not to be selected on average image
%
% Latest Update by T.Cho, Nov. 28, 2021

    if nargin < 2
        TakeOffTimes = 0;
    end
    
    spath = strcat(CaseStudy,'/s/');
    S_sum = zeros(3222,1);

    switch CaseStudy
        case '6wk'
            ntimes = 328;
        case '1yr'
            ntimes = 2920;
    end

    for i = (TakeOffTimes+1):ntimes
        load(strcat(spath,'s_',num2str(i),'.mat'));
        S_sum = S_sum + s(:);
    end
    S_avg = S_sum / (ntimes-TakeOffTimes);
    
    lon = -179.5:1:-10.5;
    lat = 10.5:1:79.5;
    Hpath = strcat(CaseStudy,'/H/');
    load(strcat(Hpath,'land_mask.mat'));
    
    S_grid = zeros(length(lon),length(lat));
    S_grid(land_mask) = S_avg;
    S_fin = flipud(S_grid');

end

