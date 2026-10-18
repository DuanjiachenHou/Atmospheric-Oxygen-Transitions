function plot_potential_well_area(r)                                                      
    clc; clear; close all;                                                                 

 
    log_r = [10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6];                   

 
    area_high = [2.84256E+35, 2.79828E+35, 2.76264E+35, 2.71296E+35, 2.64924E+35, ...        
                 2.60820E+35, 2.56932E+35, 2.50236E+35, 2.44512E+35, 2.35332E+35];          
    inc_high  = [4.42800E+33, 3.56400E+33, 4.96800E+33, 6.37200E+33, 4.10400E+33, ...       
                 3.88800E+33, 6.69600E+33, 5.72400E+33, 9.18000E+33, 2.35332E+35];           
    prev_high = area_high - inc_high;                                                          

 
    area_low = [4.64400E+33, 4.32000E+33, 3.88800E+33, 3.02400E+33, 2.64000E+33, ...          
                2.37600E+33, 2.26800E+33, 2.16000E+33, 2.05200E+33, 1.83600E+33];             
    inc_low  = [3.24000E+32, 4.32000E+32, 8.64000E+32, 3.84000E+32, 2.64000E+32, ...          
                1.08000E+32, 1.08000E+32, 1.08000E+32, 2.16000E+32, 1.83600E+33];              
    prev_low = area_low - inc_low;                                                          

 
    yLim_high = [2.30e35, 2.90e35];                                                     
    yLim_low  = [1.7e33, 5e33];                                                       

 
    figure;                                                                                  
    hold on; box on;                                                                        

    Yh = [prev_high(:), inc_high(:)];                                                  
    bh = bar(log_r, Yh, 'stacked', 'BarWidth', 0.70);                             

    bh(1).FaceColor = 'w';                                                                
    bh(1).EdgeColor = 'k';                                                                
    bh(1).LineStyle = '--';                                                            

    bh(2).FaceColor = [0 0.45 0.74];                                                    
    bh(2).EdgeColor = 'none';                                                              

    xlabel('log_{10}(r)','FontSize',12);                                                
    ylabel('Pot. well area','FontSize',12);                                                
    title('High-Oxygen Area (stack: baseline + incremental)','FontSize',14);             
    legend({'Baseline (=Area-Inc)','Incremental Area'},'Location','northwest');              

    xticks(log_r);                                                                         
    xticklabels(compose('%.1f', log_r));                                                   

    ylim(yLim_high);                                                                      

    for k = 1:numel(log_r)                                                                   
        text(log_r(k), area_high(k)*1.01, sprintf('%.2E', area_high(k)), ...               
            'HorizontalAlignment','center','FontSize',8);                                
    end                                                                                    

    set(gca,'FontSize',10,'LineWidth',1);                                               

 
    figure;                                                                                  
    hold on; box on;                                                                     

    Yl = [prev_low(:), inc_low(:)];                                                          
    bl = bar(log_r, Yl, 'stacked', 'BarWidth', 0.70);                                        

    bl(1).FaceColor = 'w';                                                                
    bl(1).EdgeColor = 'k';                                                              
    bl(1).LineStyle = '--';                                                       

    bl(2).FaceColor = [0.85 0.33 0.10];                                                  
    bl(2).EdgeColor = 'none';                                                                  
    xlabel('log_{10}(r)','FontSize',12);                                                        
    ylabel('Pot. well area','FontSize',12);                                                    
    title('Low-Oxygen Area (stack: baseline + incremental)','FontSize',14);                     
    legend({'Baseline (=Area-Inc)','Incremental Area'},'Location','northwest');                 

    xticks(log_r);                                                                             
    xticklabels(compose('%.1f', log_r));                                                       

    ylim(yLim_low);                                                                           

    for k = 1:numel(log_r)                                                                   
        text(log_r(k), area_low(k)*1.03, sprintf('%.2E', area_low(k)), ...                     
            'HorizontalAlignment','center','FontSize',8);                                      
    end                                                                                        

    set(gca,'FontSize',10,'LineWidth',1);                                                    
end                                                                                          
