function plot_potential_well_area(N) 
    clc; 
    r = 3.00E+11; 
    D = 10^29.4; 
    N_log10 = [15.0 15.1 15.2 15.3 15.4 15.5 15.6]; 
    area_high = [8.87E+33 3.57E+34 4.58E+34 8.36E+34 1.06E+35 1.64E+35 2.86E+35]; 
    area_low = [2.42E+32 4.90E+32 5.28E+32 9.36E+32 1.70E+33 2.69E+33 3.56E+33]; 
    prev_high = [0 area_high(1:end-1)]; 
    inc_high = area_high - prev_high; 
    prev_low = [0 area_low(1:end-1)]; 
    inc_low = area_low - prev_low; 
    crit1_log10 = 14.8999; 
    crit2_log10 = 15.6026; 
    fprintf('r = %.2E, D = %.6E\n', r, D); 

    figure; 
    hold on; box on; 
    plot(N_log10, area_high, 's-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'High-Oxygen Area'); 
    xline(crit1_log10, 'r--', 'LineWidth', 1.5, 'DisplayName', 'LCP'); 
    xline(crit2_log10, 'g--', 'LineWidth', 1.5, 'DisplayName', 'HCP'); 
    xlabel('log_{10}(N)', 'FontSize', 12); 
    ylabel('Potential Well Cross-Sectional Area', 'FontSize', 12); 
    legend('Location', 'northwest'); 
    xlim([14.85 15.65]); 
    set(gca, 'FontSize', 10, 'LineWidth', 1);

    figure; 
    hold on; box on; 
    bh = bar(N_log10, [prev_high(:), inc_high(:)], 'stacked'); 
    bh(1).FaceColor = 'w'; 
    bh(1).EdgeColor = 'k'; 
    bh(1).LineStyle = '--'; 
    bh(2).FaceColor = [0 0.45 0.74]; 
    bh(2).EdgeColor = 'none'; 
    xlabel('log_{10}(N)', 'FontSize', 12); 
    ylabel('Potential Well Cross-Sectional Area', 'FontSize', 12); 
    title('High-Oxygen Stable State Area', 'FontSize', 14); 
    legend({'Area at Previous N', 'Incremental Area'}, 'Location', 'northwest'); 
    for k = 1:numel(N_log10) 
        text(N_log10(k), area_high(k)*1.02, sprintf('%.2E', area_high(k)), 'HorizontalAlignment', 'center', 'FontSize', 8); 
    end 
    ylim([0 max(area_high)*1.15]); 
    set(gca, 'FontSize', 10, 'LineWidth', 1); 

    figure; 
    hold on; box on; 
    bl = bar(N_log10, [prev_low(:), inc_low(:)], 'stacked'); 
    bl(1).FaceColor = 'w'; 
    bl(1).EdgeColor = 'k'; 
    bl(1).LineStyle = '--'; 
    bl(2).FaceColor = [0.85 0.33 0.10]; 
    bl(2).EdgeColor = 'none'; 
    xlabel('log_{10}(N)', 'FontSize', 12); 
    ylabel('Potential Well Cross-Sectional Area', 'FontSize', 12); 
    title('Low-Oxygen Stable State Area', 'FontSize', 14); 
    legend({'Area at Previous N', 'Incremental Area'}, 'Location', 'northwest'); 
    for k = 1:numel(N_log10) 
        text(N_log10(k), area_low(k)*1.02, sprintf('%.2E', area_low(k)), 'HorizontalAlignment', 'center', 'FontSize', 8); 
    end 
    ylim([0 max(area_low)*1.15]); 
    set(gca, 'FontSize', 10, 'LineWidth', 1); 
end 

