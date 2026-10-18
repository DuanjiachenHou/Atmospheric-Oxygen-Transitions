function pathway_table_map_bidirectional()

clc; clear; format long; 

stdata = [
   3e11, 1*10^14.9, 1.47784e16, 4.90558e12, 1.47783e16, 3.72183e17;

];


Np = 40;                  
parameters.lambda = 3;   


colorDeepPurple = [0.4, 0, 0.6];  
colorLightPurple = [0.8, 0.7, 1]; 


data_table = struct('r', {}, 'N', {}, 'path_fwd', {}, 'path_rev', {}, 'Xi', {}, 'Xf', {});


for idx = 1:size(stdata,1)  

    r  = stdata(idx,1);     
    N  = stdata(idx,2);     
    Xi = stdata(idx,3:4)';   
    Xf = stdata(idx,5:6)';   
    parameters.vol = N;      


    [~,f1] = fminsearch(@(x)V_eff(x,parameters), Xi); 
    [~,f2] = fminsearch(@(x)V_eff(x,parameters), Xf); 
    parameters.Eeff = -min([f1,f2]);               
    fprintf('Case %d: r=%.2e, N=%.2e, Eeff=%.3e\n', idx, r, N, parameters.Eeff); 


    parameters.Xi = Xi;     
    parameters.Xf = Xf;      

    x0_mat_fwd = zeros(2, Np); 
    for j = 1:Np
        x0_mat_fwd(:,j) = Xi + j*(Xf - Xi)/(Np+1);
    end
    x0_fwd = x0_mat_fwd(:);   
    opts = optimset('Display','off','MaxFunEvals',1e50,'MaxIter',500);
    [x_flat_fwd, ~, ~] = fminsearch(@(x) SHJ_WKB(reshape(x,2,Np), parameters), x0_fwd, opts);
    path_fwd = [ Xi, reshape(x_flat_fwd,2,[]), Xf ]; 


    parameters.Xi = Xf;     
    parameters.Xf = Xi;     
    x0_mat_rev = zeros(2, Np);
    for j = 1:Np
        x0_mat_rev(:,j) = Xf + j*(Xi - Xf)/(Np+1);
    end
    x0_rev = x0_mat_rev(:);  
    [x_flat_rev, ~, ~] = fminsearch(@(x) SHJ_WKB(reshape(x,2,Np), parameters), x0_rev, opts);
    path_rev = [ Xf, reshape(x_flat_rev,2,[]), Xi ]; 


    data_table(idx).r = r;              
    data_table(idx).N = N;              
    data_table(idx).Xi = Xi;             
    data_table(idx).Xf = Xf;             
    data_table(idx).path_fwd = path_fwd; 
    data_table(idx).path_rev = path_rev;


    figure('Name',sprintf('MAP 双向 r=%.2e, N=%.2e',r,N),'NumberTitle','off');
    plot(path_fwd(1,:), path_fwd(2,:), '-o', 'Color', colorDeepPurple, 'LineWidth',1); hold on;
    plot(path_rev(1,:), path_rev(2,:), '-s', 'Color', colorLightPurple, 'LineWidth',1);
    plot(Xi(1),Xi(2),'ks','MarkerFaceColor','g','DisplayName','Xi (低稳态)');
    plot(Xf(1),Xf(2),'ks','MarkerFaceColor','m','DisplayName','Xf (高稳态)');
    legend('正向: 低→高','反向: 高→低','Location','best'); xlabel('X'); ylabel('Y');
    title(sprintf('MAP 双向路径 (r=%.2e, N=%.2e)', r, N)); grid on;
    hold off;
end


save('map_paths_0725.mat', 'data_table'); 

end  


function Fout = F(x, parameters)
    [Fx,Fy] = force_XY(x(1), x(2), parameters.vol);
    Fout    = [Fx; Fy];
end

function V = V_eff(x, parameters)
    Fv = F(x, parameters);
    Dm = Dij(x, parameters);
    V  = 0.5 * Fv' * (Dm \ Fv);
end

function D = Dij(x, parameters)
    D = eye(2);
end

function [Fx, Fy] = force_XY(X, Y, N)
    s   = 2.03e-5; dr = 1.36e19; dd = 2.73e17;
    a1  = 0.003;    a2 = -0.1655; a3 = 3.2305;
    a4  = -25.8343; a5 = 71.6;
    delta = Y./(dd + Y);
    gamma = Y./(dr + Y);
    omega = (1 - gamma).*(1 - delta);
    p     = log10(Y);
    fai   = 10.^(a1.*p.^4 + a2.*p.^3 + a3.*p.^2 + a4.*p + a5);
    Fx    = 0.5 .* omega .* (2*N) - s.*X - 0.5 .* fai .* (X.^0.7);
    Fy    = omega .* N   - (1-omega).*N - s.*X - fai .* (X.^0.7);
end

function S = SHJ_WKB(x, parameters)
    Eeff   = parameters.Eeff;
    Xi     = parameters.Xi;
    Xf     = parameters.Xf;
    lambda = parameters.lambda;
    Npt    = size(x,2);

    X  = [Xi, x, Xf];
    Xc = (X(:,1:end-1) + X(:,2:end))/2;
    lsum = sum(vecnorm(X(:,2:end)-X(:,1:end-1)));
    lave = lsum/(Npt+1);

    s = zeros(Npt+1,1);
    for i = 1:Npt+1
        dX   = X(:,i+1) - X(:,i);
        dl0  = dX'*dX;
        dl   = dX'*(Dij(Xc(:,i),parameters)\dX);
        fdl  = dX'*(Dij(Xc(:,i),parameters)\F(Xc(:,i),parameters));
        Veff = V_eff(Xc(:,i), parameters);
        s(i) = sqrt(dl)*sqrt(Veff+Eeff) - fdl/2 + lambda*(sqrt(dl0)-lave)^2;
    end
    S = sum(s);
end
