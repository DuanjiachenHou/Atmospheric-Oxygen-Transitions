function main
    
    stablePt1 = [1.47784e16, 4.90558E+12];
    stablePt2 = [1.47783e16, 3.72183E+17];

    
    xmin = [1.45e16, 1e12];
    xmax = [1.50e16, 5e17];  

    N           = 1*10^14.9;   
    numPoints   = 10000;      
    tolerance1  = 6e12;    
    tolerance2  = 3.5e17;     
    tspan       = [0, 1e9]; 

    
    options = odeset('RelTol',1e-2, 'AbsTol',[1e6, 1e6], ...
                     'Events', @eventsFunc);

    
    figure('Name','相图与随机撒点','NumberTitle','off');
    xVec = linspace(xmin(1), xmax(1), 400);
    yVec = linspace(xmin(2), xmax(2), 400);
    [X, Y] = meshgrid(xVec, yVec);

   
    [Fx, Fy] = force_XY(X, Y, N);

   
    contour(X, Y, Fx, [0 0], 'b');
    hold on;
    contour(X, Y, Fy, [0 0], 'k');

    xlabel('x1');
    ylabel('x2');
    title('Fx=0 & Fy=0 等高线 (相图)');

    
    plot(stablePt1(1), stablePt1(2), 'ro', 'MarkerSize', 10, 'LineWidth',2, ...
         'DisplayName','稳点1');
    plot(stablePt2(1), stablePt2(2), 'bo', 'MarkerSize', 10, 'LineWidth',2, ...
         'DisplayName','稳点2');
    legend('Fx=0','Fy=0','稳点1','稳点2','Location','best');
    ylim([0, xmax(2)]);

   
    rng('shuffle'); % 根据时间种子随机
    randomPoints = [ ...
        xmin(1) + (xmax(1)-xmin(1)) * rand(numPoints,1), ...
        10.^(log10(xmin(2)) + (log10(xmax(2)) - log10(xmin(2))) * rand(numPoints,1)) ...
    ];

    count1 = 0;  % 流向稳点1的计数
    count2 = 0;  % 流向稳点2的计数

   
    func = @(t,x) force_MO(t, x, N);

    
    for i = 1:numPoints
        x0 = randomPoints(i,:);

        
        [~, Xt, te, Xe, ie] = ode15s(func, tspan, x0, options);

        
        final = Xt(end,:);

        % 判断与两个稳点的距离
        dist1 = norm(final - stablePt1);
        dist2 = norm(final - stablePt2);

        if dist1 < dist2 && dist1 < tolerance1
            count1 = count1 + 1;
            plot(final(1), final(2), 'r.', 'MarkerSize', 6);
        elseif dist2 < dist1 && dist2 < tolerance2
            count2 = count2 + 1;
            plot(final(1), final(2), 'b.', 'MarkerSize', 6);
        else
         
            plot(final(1), final(2), 'k.', 'MarkerSize', 4);
        end
    end

    
    prob1 = count1 / numPoints;
    prob2 = count2 / numPoints;
    fprintf('---------------------------\n');
    fprintf('共撒点 %d 个\n', numPoints);
    fprintf('流向稳点1 (%.4g, %.4g) 的概率：%.2f%%\n', ...
            stablePt1(1), stablePt1(2), prob1 * 100);
    fprintf('流向稳点2 (%.4g, %.4g) 的概率：%.2f%%\n', ...
            stablePt2(1), stablePt2(2), prob2 * 100);
    fprintf('未收敛到任一稳点的点数：%d\n', numPoints - (count1 + count2));
    fprintf('---------------------------\n');

end



function f = force_MO(t,x,para)
N = para;
% N = 3.75e15;
% N = 1e15;
r = 3e11;
s = 2.03e-5;
% w = 6e-9;
dr = 1.36e19;
dd = 2.73e17;
a1 = 0.003; a2=-0.1655; a3=3.2305; a4=-25.8343; %a5=71.5398;
a5=71.6;
delta = x(2)/(dd+x(2));
gamma = x(2)/(dr+x(2));
omega = (1-gamma)*(1-delta);
p = log10(x(2));
fai = 10^(a1*p^4+a2*p^3+a3*p^2+a4*p+a5);
f = [0.5*omega*(N+r)-s*x(1)-0.5*fai*x(1)^0.7;
    omega*N-(1-omega)*r-s*x(1)-fai*x(1)^0.7];
end


function [Fx, Fy] = force_XY(X, Y, N)
    s   = 2.03e-5;
    dr  = 1.36e19;
    dd  = 2.73e17;
    a1  = 0.003; 
    a2  = -0.1655; 
    a3  = 3.2305; 
    a4  = -25.8343; 
    a5  = 71.6;

    Y(Y <= 0) = 1e-10;

    delta = Y ./ (dd + Y);
    gamma = Y ./ (dr + Y);
    omega = (1 - gamma) .* (1 - delta);

    p = log10(Y);
    p = max(min(p, 18), 13); % 限制 p 的范围

    fai = 10.^(a1.*p.^4 + a2.*p.^3 + a3.*p.^2 + a4.*p + a5);

    Fx = 0.5 .* omega .* (2*N) - s.*X - 0.5 .* fai .* (X.^0.7);
    Fy = omega .* N - (1 - omega).*N - s.*X - fai .* (X.^0.7);
end


function [value, isterminal, direction] = eventsFunc(t, x)
    stablePt1 = [1.47782e16, 2.24673e13];
    stablePt2 = [1.47783e16, 2.86041e18];
    tolerance1 = 6e12;   %1e14;  tolerance1  = 6e12;       
    tolerance2 = 3.5e17;   %1e16;  tolerance2  = 3.5e17;

    dist1 = norm(x - stablePt1);
    dist2 = norm(x - stablePt2);

    value = [dist1 - tolerance1; dist2 - tolerance2];
    isterminal = [1; 1]; 
    direction = [0; 0];   
end
