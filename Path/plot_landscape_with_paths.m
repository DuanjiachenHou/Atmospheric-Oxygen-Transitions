function plot_landscape_with_paths_2D                                   
 

    clc; clear; format long;                                              

 
    landscapeFile   = 'data14.9_11.txt';  
  
    smoothingSigma  = 7;                                                   
    upsampleK       = 3;                                                      
    bgAlpha         = 0.70;                                                  
    useManualCLim   = true;                                                  
    climMin         = [];                                                     
    climMax         = 3;                                                     
    cmapName        = 'jet';                                                
    logClamp        = 0;                                                  

 
    xlimPR = [0, 4e17];                                                      
    ylimPSS = [1.36e16, 1.52e16];                                          

 
    resampleN        = 400;                                                
    pathSmoothSpan   = 1;                                                 
    lockEndpoints    = true;                                                
    pathLineWidth    = 2;                                                    
    markerSize       = 4;                                                    

 
    scalePath        = 1;                                                    
    shiftPR          = 0;                                                    
    shiftPSS         = 0;                                                   

 
    data = readmatrix(landscapeFile);                                       
    PSS  = data(:,1);                                                      
    PR   = data(:,2);                                                         
    uabs = data(:,3);                                                      
    yv = unique(PSS);                                                      
    xv = unique(PR);                                                         

 
    Z = reshape(-log(uabs), numel(yv), numel(xv));                     
    if smoothingSigma > 0                                                   
        Z = imgaussfilt(Z, smoothingSigma);                                    
    end
    if logClamp > 0                                                           
        Z = min(Z, logClamp);                                                 
    end

 
    if upsampleK > 1                                                         
        xv_hi = linspace(xv(1), xv(end), max(2, numel(xv)*upsampleK));          
        yv_hi = linspace(yv(1), yv(end), max(2, numel(yv)*upsampleK));       
        Z_hi  = imresize(Z, [numel(yv_hi), numel(xv_hi)], 'bicubic');          
    else                                                                       
        xv_hi = xv; yv_hi = yv; Z_hi = Z;                                     
    end

 
    dx = median(diff(xv_hi)); if isempty(dx), dx = 1; end                    
    dy = median(diff(yv_hi)); if isempty(dy), dy = 1; end                     
    xEdges = [xv_hi(1)-dx/2, xv_hi(end)+dx/2];                                
    yEdges = [yv_hi(1)-dy/2, yv_hi(end)+dy/2];                                

 
    figure('Name','Landscape 2D with Paths: Zoomed-In Version', ...          
           'NumberTitle','off');                                              
    ax = gca;                                                                  
    set(ax,'YDir','normal');                                                 
    set(ax,'SortMethod','childorder');                                        
    hold(ax,'on');                                                          
    colormap(ax, cmapName);                                                     

 
    hBG = imagesc(ax, 'XData', xEdges, 'YData', yEdges, 'CData', Z_hi);        
    set(hBG,'AlphaData', bgAlpha, ...                                          
            'HitTest','off', 'PickableParts','none');                        

 
 
ZminShow = min(Z_hi(:)); ZmaxShow = max(Z_hi(:));                       
if useManualCLim                                                         
    if isempty(climMin)                                                  
        climMinUse = ZminShow;                                         
    else                                                                   
        climMinUse = climMin;                                            
    end
    caxis(ax, [climMinUse, climMax]);                                   
    cbTicks = linspace(climMinUse, climMax, 4);                        
else                                                                     
    caxis(ax, [ZminShow, ZmaxShow]);                                   
    cbTicks = linspace(ZminShow, ZmaxShow, 4);                            
end

cb = colorbar(ax);                                                       
ylabel(cb,'-log(|u|)');                                                
set(cb,'Ticks',cbTicks);                                                

cbTickLabels = arrayfun(@(v) sprintf('%.0f', v), cbTicks, ...            
                        'UniformOutput', false);                        
set(cb,'TickLabels', cbTickLabels);                                       

set(cb,'FontSize',10);                                                   



 
    xRangeShow = [max(xlimPR(1), xv_hi(1)), min(xlimPR(2), xv_hi(end))];       
    yRangeShow = [max(ylimPSS(1), yv_hi(1)), min(ylimPSS(2), yv_hi(end))];    
    xlim(ax, xRangeShow);                                                     
    ylim(ax, yRangeShow);                                                       
 
 
 
xticks([0 1e17 2e17 3e17]);                     
yticks([1.36e16 1.40e16 1.44e16 1.48e16 1.52e16]);       
 

 
    set(ax,'FontSize',10,'LineWidth',1);     

 
    uistack(hBG,'bottom');                                                   

    xlabel(ax,'P_R (x)');                                                    
    ylabel(ax,'P_{SS} (y)');                                                  
 

 
    load('map_paths_0725.mat','data_table');                               
    colorFwd = [0.35, 0.00, 0.70];                                       
    colorRev = [0.85, 0.75, 1.00];                                             

    hTop = gobjects(0);                                                    
    for idx = 1:numel(data_table)                                         
        P = data_table(idx).path_fwd;                                         
        R = data_table(idx).path_rev;                                         

 
        PR_fwd_raw  = scalePath*(P(1,:) + shiftPR);                           
        PSS_fwd_raw = scalePath*(P(2,:) + shiftPSS);                      
        PR_rev_raw  = scalePath*(R(1,:) + shiftPR);                            
        PSS_rev_raw = scalePath*(R(2,:) + shiftPSS);                            

 
        [PR_fwd, PSS_fwd] = resample_curve_equal_arc(PR_fwd_raw, PSS_fwd_raw, ...
                                                     resampleN, lockEndpoints); 
        [PR_rev, PSS_rev] = resample_curve_equal_arc(PR_rev_raw, PSS_rev_raw, ...
                                                     resampleN, lockEndpoints); 

 
        if pathSmoothSpan > 1                                                
            win = max(3, 2*floor(pathSmoothSpan/2)+1);                      
            PR_fwd  = movmean(PR_fwd,  win);                                
            PSS_fwd = movmean(PSS_fwd, win);                               
            PR_rev  = movmean(PR_rev,  win);                                  
            PSS_rev = movmean(PSS_rev, win);                                
            if lockEndpoints                                                   
                PR_fwd(1) = PR_fwd_raw(1);  PSS_fwd(1) = PSS_fwd_raw(1);       
                PR_fwd(end) = PR_fwd_raw(end); PSS_fwd(end) = PSS_fwd_raw(end);
                PR_rev(1) = PR_rev_raw(1);  PSS_rev(1) = PSS_rev_raw(1);       
                PR_rev(end) = PR_rev_raw(end); PSS_rev(end) = PSS_rev_raw(end);
            end
        end

 
        hF = plot(ax, PSS_fwd, PR_fwd, '-', 'Color', colorFwd, ...          
                  'LineWidth', pathLineWidth);                                 
        hR = plot(ax, PSS_rev, PR_rev, '-', 'Color', colorRev, ...            
                  'LineWidth', pathLineWidth);                                  
        hFS = plot(ax, PSS_fwd(1),  PR_fwd(1),  'p', ...                      
                   'MarkerEdgeColor','k','MarkerFaceColor',colorFwd, ...
                   'MarkerSize',markerSize+2);                                
        hFE = plot(ax, PSS_fwd(end), PR_fwd(end),'h', ...                     
                   'MarkerEdgeColor','k','MarkerFaceColor',colorFwd, ...
                   'MarkerSize',markerSize+2);                                 
        hRS = plot(ax, PSS_rev(1),  PR_rev(1),  'p', ...                      
                   'MarkerEdgeColor','k','MarkerFaceColor',colorRev, ...
                   'MarkerSize',markerSize+2);                                  
        hRE = plot(ax, PSS_rev(end), PR_rev(end),'h', ...                      
                   'MarkerEdgeColor','k','MarkerFaceColor',colorRev, ...
                   'MarkerSize',markerSize+2);                                 

        hTop = [hTop, hF, hR, hFS, hFE, hRS, hRE];                              
    end

 
    if ~isempty(hTop)                                                          
        uistack(hTop,'top');                                                 
    end

 
    legend({'Forward Path','Reverse Path','Forward Start', ...                
            'Forward End','Reverse Start','Reverse End'}, ...
           'Location','best', 'FontSize',9);                                 
end

function [xs, ys] = resample_curve_equal_arc(x, y, N, lockEnds)           
 
    x = x(:)'; y = y(:)';                                                     
    if numel(x) < 2                                                            
        xs = x; ys = y; return;                                            
    end
    ds = hypot(diff(x), diff(y));                                           
    s  = [0, cumsum(ds)];                                                    
    if s(end) == 0                                                            
        xs = x; ys = y; return;                                               
    end
    su = linspace(0, s(end), max(N, numel(x)));                              
    xs = pchip(s, x, su);                                                    
    ys = pchip(s, y, su);                                                   
    if lockEnds                                                            
        xs(1) = x(1); ys(1) = y(1);                                          
        xs(end) = x(end); ys(end) = y(end);                                  
    end
end
