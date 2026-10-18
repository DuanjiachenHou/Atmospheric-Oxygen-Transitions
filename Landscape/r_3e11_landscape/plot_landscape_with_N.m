 
data     = readmatrix('data15.2_11_29.4.txt');


u_abs    = data(:,3);                                                 
F1_flat  = data(:,6);                                                  
F2_flat  = data(:,7);                                                 

x_flat   = data(:,2);                                                
y_flat   = data(:,1);                                                 

 
xv = unique(x_flat);                                                  
yv = unique(y_flat);                                                  
[X, Y] = meshgrid(xv, yv);                                            
Z  = reshape(-log(u_abs), numel(yv), numel(xv));                       

 
Z_top   = 0;                                                          
Z_clip  = min(Z, Z_top);                                               

 
figure;                                                               
ax = gca;                                                              
set(ax,'SortMethod','childorder');                                    

hSurf = surf(X, Y, Z_clip, ...                                       
             'FaceColor','interp', ...                                 
             'EdgeColor','none');                                      
hold on;                                                               

 
z_min_data = min(Z_clip(:));                                           
zAxisMin   = min(z_min_data, -3);                                     
zAxisMax   = 2;                                                 
zlim([zAxisMin, zAxisMax]);                                           

 ax.ZTick = [-2 0 2];                                           

 
colormap(jet);                                                        
caxis([z_min_data, Z_top]);                                            
 
shading interp;                                                       
box off;                                                              

 
set(ax,'FontSize',14,'LineWidth',1.2);                                 

xlabel('P_R','FontSize',18);                                          
ylabel('P_{SS}','FontSize',18);                                       
zlabel('U','FontSize',18);                                            

 
pos = get(ax,'Position');                                              
pos(2) = pos(2) + 0.06;                                               
set(ax,'Position',pos);                                               

 
bottomTitle = 'log_{10}(N)=15.3';                                     
text('Units','normalized', ...                                         
     'Position',[0.5, -0.08, 0], ...                                   
     'String', bottomTitle, ...                                       
     'HorizontalAlignment','center', ...                              
     'VerticalAlignment','top', ...                                   
     'FontSize',18);                                                   
view(35,20);                                                          
hold off;                                                              
