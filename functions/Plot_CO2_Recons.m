function Plot_CO2_Recons(SHAT,TITLE_NAME,cscale,LAND,LAT,LON)
% This function is to generate map of estimated fluxes
% 
%   SHAT       - array representing average fluxes
%   TITLE_NAME - title of image
%   cscale     - color scale containing min and max values
%   LAND       - terrestrial regions of North America 
%   LAT        - range of lattitude
%   LON        - range of longitude
%
% Latest Update by T.Cho, Nov. 28, 2021

    % Need cspy_edited function
    cspy_edited(SHAT,'colormap','gray','levels',8,'caxis',cscale)
    xlabel('')
    
    % title
    title(TITLE_NAME,'FontSize',15), 

    % colorbar
    colorbar('FontSize',16)
    caxis([min(cscale(:)) max(cscale(:))])
    
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
    lon_ticks_location = linspace(0,n,num_ticks);
    lon_ticks_location(1) = 1;
    
    yticks(lat_ticks_location)
    yticklabels(num2cell(round(lat_ticks)))
    
    xticks(lon_ticks_location)
    xticklabels(num2cell(round(lon_ticks)))
    
    set(gca,'XAxisLocation','bottom','LineWidth',1)    
    
    % Boundary 
	S_BW = zeros(n,m);
    S_BW(LAND)= 1;
    S_BW = flipud(S_BW');
    BW = bwboundaries(S_BW,'holes');

    hold on
    for k = 1:length(BW)
        boundary = BW{k};
        plot(boundary(:,2),boundary(:,1), 'k', 'Linewidth',1.5)
    end
    hold off