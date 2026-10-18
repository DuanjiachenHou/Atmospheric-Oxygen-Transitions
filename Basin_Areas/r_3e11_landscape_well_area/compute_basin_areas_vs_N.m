 
data     = readmatrix('data15.5_11_29.4.txt');         

 
x_flat   = data(:,2);                                                
y_flat   = data(:,1);                                                 
u_abs    = data(:,3);                                                   
F1_flat  = data(:,6);                                                 
F2_flat  = data(:,7);                                                 

 
xv = unique(x_flat);                                                  
yv = unique(y_flat);                                                 
[X, Y] = meshgrid(xv, yv);                                             
Z   = reshape(-log(u_abs), numel(yv), numel(xv));                      
F1  = reshape(F1_flat,   numel(yv), numel(xv));                         
F2  = reshape(F2_flat,   numel(yv), numel(xv));                          

 
threshold = 10;                                                   
Z_mask = Z;                                                          
Z_mask(Z_mask > threshold) = NaN;                                      
figure;                                                                
surf(X, Y, Z_mask, 'EdgeColor','none');                               
hold on;                                                               
colormap(jet); caxis([0, threshold]); colorbar;                     
title('挖空阈值展示（Z \le threshold 部分）');                            
xlabel('O'); ylabel('M'); zlabel('-log(u)');            

 

upper = 0;                                                          
 
Z_mask2 = min(Z, upper);                                                
figure;                                                              
surf(X, Y, Z_mask2, 'FaceColor','interp', 'EdgeColor','none');         
hold on;                                                               
lower = min(Z_mask2(:));                                               
zlim([lower, upper]);                                                  
colormap(jet); caxis([lower, upper]); colorbar;                       
shading interp;                                                       
title('钳制上限后的连续景观');                                           
xlabel('O'); ylabel('M'); zlabel('-log(u)');            

 
N = 1e-6;                                                             
followSign = +1;                                                      
 
odefun = @(t,xy) [ ...                                                   
    followSign * interp2(xv, yv, F1, xy(1), xy(2), 'linear', 0); ...     
    followSign * interp2(xv, yv, F2, xy(1), xy(2), 'linear', 0) ...     
];                                                                    

  
tspan   = [0, 100];                                                     
x_start = xv(1);                                                       
y_start = yv(1);                                                        
xMin = min(xv); xMax = max(xv); yMin = min(yv); yMax = max(yv);         
eventOut = @(t,xy) boxExitEvent(t,xy,xMin,xMax,yMin,yMax);              
opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'Events',eventOut);           

[~, XY] = ode45(odefun, tspan, [x_start; y_start], opts);            
x_path = XY(:,1);                                                      
y_path = XY(:,2);                                                        
z_path = interp2(xv, yv, Z,  x_path, y_path, 'linear', NaN);           

 
plot3(x_path, y_path, z_path, 'w-', 'LineWidth', 2);                    
scatter3(x_path(1), y_path(1), z_path(1), 60, 'g', 'filled');           
scatter3(x_path(end), y_path(end), z_path(end), 60, 'r', 'filled');      
view(30, 15);                                                            
title('Probability Landscape with Path (基于插值漂移场)');               

 
 
 
 
hcut = -0.3;                                                            
 
Zsub = Z <= hcut;                                                     
CC = bwconncomp(Zsub, 4);                                              
numBasins = CC.NumObjects;                                             
if numBasins < 2                                                       
    warning('hcut=%.4g 仅检测到 %d 个坑，请调整 hcut。', hcut, numBasins); 
end                                                                      

dx = mean(diff(xv));                                                   
dy = mean(diff(yv));                                                  
cellArea = dx * dy;                                                     

areas = zeros(numBasins,1);                                             
for k = 1:numBasins                                                    
    areas(k) = numel(CC.PixelIdxList{k}) * cellArea;                    
end                                                                      
[areasSorted, order] = sort(areas, 'descend');                          
pick = min(2, numBasins);                                               
chosenIdx = order(1:pick);                                              

fprintf('=== 切面 z = %.6g 的两态面积 ===\n', hcut);                     
for i = 1:pick                                                                             
    fprintf('坑 #%d（ID=%d）面积 = %.6e  (像素=%d, 单元=%.3e)\n', ...  
        i, chosenIdx(i), areasSorted(i), numel(CC.PixelIdxList{chosenIdx(i)}), cellArea); 
end                                                                    
if numBasins > 2                                                    
    fprintf('注意：总共检测到 %d 个坑，已按面积选取最大的两个。\n', numBasins); 
end                                                                     

 
figure;                                                                
contourf(X, Y, Z, 30, 'LineStyle','none');                            
hold on; colormap(jet); colorbar;                                       
contour(X, Y, Z, [hcut hcut], 'k--', 'LineWidth', 1.5);                  
L = labelmatrix(CC);                                                     
for i = 1:pick                                                          
    mask_i = (L == chosenIdx(i));                                      
    B = bwboundaries(mask_i);                                            
    maxLen = 0; maxb = 1;                                               
    for b = 1:numel(B)                                                  
        if size(B{b},1) > maxLen, maxLen = size(B{b},1); maxb = b; end   
    end                                                                  
    rrcc = B{maxb}; rr = rrcc(:,1); cc = rrcc(:,2);                      
    px = X(sub2ind(size(X), rr, cc));                                  
    py = Y(sub2ind(size(Y), rr, cc));                                    
    patch(px, py, i, 'FaceAlpha', 0.35, 'EdgeColor', 'none');            
    text(mean(px,'omitnan'), mean(py,'omitnan'), ...                  
        sprintf('坑#%d\\nA=%.2e', i, areasSorted(i)), ...               
        'Color','[0.85 0.85 0.85]','FontWeight','bold','HorizontalAlignment','center'); 
 
 
end                                                                  
xlabel('O'); ylabel('M');                             
title(sprintf('log(N)=15.5', hcut));               
hold off;                                                  

 
function [value, isterminal, direction] = boxExitEvent(~, xy, xmin, xmax, ymin, ymax) 
    x = xy(1); y = xy(2);                                              
 
    value = [x - xmin; xmax - x; y - ymin; ymax - y];                  
    isterminal = [1; 1; 1; 1];                                         
    direction = [0; 0; 0; 0];                                         
end                                                                  
