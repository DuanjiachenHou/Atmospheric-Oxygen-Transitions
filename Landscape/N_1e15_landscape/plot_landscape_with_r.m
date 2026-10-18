 
dataFile = 'data_r_11.0_29.5.txt'; 


data     = readmatrix(dataFile);                                          
x_flat   = data(:,1);                                                   
y_flat   = data(:,2);                                                  
u_abs    = data(:,3);                                                 
F1_flat  = data(:,6);                                                  
F2_flat  = data(:,7);                                                  

 
xv = unique(x_flat);                                                   
yv = unique(y_flat);                                                    
[X, Y] = meshgrid(xv, yv);                                      
Z   = reshape(-log(u_abs), numel(yv), numel(xv));                      
F1  = reshape(F1_flat,       numel(yv), numel(xv));                    
F2  = reshape(F2_flat,       numel(yv), numel(xv));                      

 
Z_top  = 0;                                                             
Z_clip = min(Z, Z_top);                                                

 
figure('Color','w');                                                    
ax = gca;                                                              
set(ax,'SortMethod','childorder');                                    

hSurf = surf(X, Y, Z_clip, ...                                          
    'FaceColor','interp', ...                                           
    'EdgeColor','none');                                                  
hold on;                                                              

 
 
z_min_data = min(Z_clip(:));                                           
zAxisMin   = min(z_min_data, -2);                                        
zAxisMax   = 2.000001;                                                
zlim([zAxisMin, zAxisMax]);                                            
ax.ZTick = [-2 0 2];                                                    

 
 
colormap(jet);                                                      
caxis([z_min_data, Z_top]);                                          
 
shading interp;                                               
box off;                                                            

 
set(ax,'FontSize',14,'LineWidth',1.2);                                  
xlabel('O','FontSize',18);                                           
ylabel('M','FontSize',18);                                         
zlabel('U','FontSize',18);                                               

 
pos = get(ax,'Position');                                               
pos(2) = pos(2) + 0.06;                                                  
set(ax,'Position',pos);                                                 

 
bottomTitle = 'log_{10}(r)=11.0';                                                  
text('Units','normalized', ...                                          
     'Position',[0.5, -0.08, 0], ...                                      
     'String', bottomTitle, ...                                         
     'HorizontalAlignment','center', ...                                 
     'VerticalAlignment','top', ...                                       
     'FontSize',18);                                                      

 
force = @(t,xy) [ ...                                                     
    interp2(X, Y, F1, xy(1), xy(2)); ...                                 
    interp2(X, Y, F2, xy(1), xy(2))  ...                               
];                                                                      

 
x_start = 1.47783e16;                                                    
y_start = 5.1745e13;                                                   
x_stop  = 1.47783e16;                                                  
y_stop  = 4.20673e18;                                                    

 
view(35, 20);                                                       
hold off;                                                           
