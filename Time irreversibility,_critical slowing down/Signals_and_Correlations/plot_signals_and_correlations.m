clear
h = 0.01; 
t = 0:h:30000; 
N = length(t); 
Nc=floor(N/10);
num_runs = 100; 
tf = 500; 

 
max_lag_length = 2 * N - 1; 

 
autoCorr_times_lgx = zeros(1, num_runs);
autoCorr_times_lgy = zeros(1, num_runs);
crossCorr_differences = zeros(1, num_runs);

 
avg_autoCorrX = zeros(1, max_lag_length);
avg_autoCorrY = zeros(1, max_lag_length);
avg_crossCorrXY = zeros(1, max_lag_length);
avg_crossCorrYX = zeros(1, max_lag_length);


    x0 = 1.487e16;
    y0 = 2.0976e+18; 
[x1,y1,t1] = stablepoint(x0,y0);
N1=length(t1)-10;



for run = 1:num_runs
 
    x = zeros(1, N); 
    y = zeros(1, N); 
    x(1) = x1(N1);
    y(1) = y1(N1); 

    D = 1 * 10^20.4; 
    para = 10^(15.2); 

 
    for i = 1:N-1
        [Fx, Fy] = force_XY(x(i), y(i), para); 
        x(i+1) = x(i) + (Fx)*h + sqrt(2*D*h)*randn(); 
        y(i+1) = y(i) + (Fy)*h + sqrt(2*D*h)*randn(); 

        if x(i+1) < 0
            x(i+1) = -x(i+1); 
        end
        if y(i+1) < 0
            y(i+1) = -y(i+1); 
        end
    end

 
    lgx = log(x + 1e16);  
    lgy = log(y + 1e16);

 
    lgx_zero_mean = lgx - mean(lgx);
    lgy_zero_mean = lgy - mean(lgy);

    stdlgx=std(lgx_zero_mean);
    stdlgy=std(lgy_zero_mean);


    [autoCorrX, lagX] = xcorr(lgx_zero_mean, lgx_zero_mean(1:Nc)); 
    [autoCorrY, lagY] = xcorr(lgy_zero_mean, lgy_zero_mean(1:Nc)); 

 
    [crossCorrXY, lagXY] = xcorr(lgx_zero_mean, lgy_zero_mean(1:Nc)); 
    [crossCorrYX, lagYX] = xcorr(lgy_zero_mean, lgx_zero_mean(1:Nc)); 

 
    autoCorrX_padded = zeros(1, max_lag_length);
    autoCorrX_padded(1:length(autoCorrX)) = autoCorrX/Nc;
    
    autoCorrY_padded = zeros(1, max_lag_length);
    autoCorrY_padded(1:length(autoCorrY)) = autoCorrY/Nc;
    
    crossCorrXY_padded = zeros(1, max_lag_length);
    crossCorrXY_padded(1:length(crossCorrXY)) = crossCorrXY/Nc;
    
    crossCorrYX_padded = zeros(1, max_lag_length);
    crossCorrYX_padded(1:length(crossCorrYX)) = crossCorrYX/Nc;

 
    avg_autoCorrX = avg_autoCorrX + autoCorrX_padded;
    avg_autoCorrY = avg_autoCorrY + autoCorrY_padded;
    avg_crossCorrXY = avg_crossCorrXY + crossCorrXY_padded;
    avg_crossCorrYX = avg_crossCorrYX + crossCorrYX_padded;

 
    diff_near_zero = (crossCorrXY(lagXY >= 0 & lagXY <= tf) - crossCorrYX(lagYX >= 0 & lagYX <= tf)).^2;
    integral_value = trapz(lagXY(lagXY >= 0 & lagXY <= tf), diff_near_zero); 
    crossCorr_differences(run) = sqrt(integral_value / tf); 
end

 
avg_autoCorrX = avg_autoCorrX / num_runs;
avg_autoCorrY = avg_autoCorrY / num_runs;
avg_crossCorrXY = avg_crossCorrXY / num_runs;
avg_crossCorrYX = avg_crossCorrYX / num_runs;

 
positive_indices = (max_lag_length + 1)/2 : max_lag_length;
positive_lags = 0:(length(positive_indices) - 1); 

 
figure;
 
 


 
subplot(3,1,1);
plot(positive_lags, lgx_zero_mean/stdlgx, 'r', 'DisplayName', 'X Signal'); 
hold on;
plot(positive_lags, lgy_zero_mean/stdlgy, 'b', 'DisplayName', 'Y Signal'); 
hold off;

 
title('Standardized Signals of X and Y');

xlabel('Time Lag');
ylabel('Standardized Value');

 
legend('Location', 'best');


subplot(3,1,2);
plot(positive_lags, avg_autoCorrX(positive_indices)/stdlgx^2, 'r', 'DisplayName', 'Autocorrelation of X'); 
hold on;
plot(positive_lags, avg_autoCorrY(positive_indices)/stdlgy^2, 'b', 'DisplayName', 'Autocorrelation of Y'); 
hold off;
title('Average Autocorrelation Function');
xlabel('Time Lag');
ylabel('Autocorrelation');
legend;

subplot(3,1,3);
plot(positive_lags, avg_crossCorrXY(positive_indices), 'r', ...
     positive_lags, avg_crossCorrYX(positive_indices), 'b');
title('Average Cross-correlation Function');
legend('Forward C_{XY}', 'Reverse C_{YX}');
xlabel('Time Lag');
ylabel('Cross-correlation');

 
 