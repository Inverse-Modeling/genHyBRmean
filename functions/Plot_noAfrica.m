function Plot_noAfrica(SHAT,d)
% This function is to remove Africa from reconstructions
% 
%   SHAT - array before clipping
%      d - degree from east to be clipped
%
% Latest Update by T.Cho, Nov. 28, 2021

    LON = -179.5:1:-10.5-d;
    LAT = 10.5:1:79.5;
    % Xticks, Yticks 
    [m,n] = size(SHAT);
    num_ticks = 4;
    
    LAT_min = floor(min(LAT(:)));
    LAT_max = ceil(max(LAT(:)));
    LON_min = floor(min(LON(:)));
    LON_max = ceil(max(LON(:)));
    
    lat_ticks = linspace(LAT_max,LAT_min,num_ticks);
    lon_ticks = linspace(LON_min,LON_max,num_ticks);
    lat_ticks_location = linspace(0,m,num_ticks);
    lat_ticks_location(1) = 1;
    lon_ticks_location = linspace(0,n-d,num_ticks);
    lon_ticks_location(1) = 1;
    
    yticks(lat_ticks_location)
    yticklabels(num2cell(round(lat_ticks)))
    
    xticks(lon_ticks_location)
    xticklabels(num2cell(round(lon_ticks)))
    
    set(gca,'XAxisLocation','bottom','LineWidth',1)
    xlim([1 length(LON)])