
function buzsaki_sim

    y0 = [0 0 ];
    
    [T,Y] = ode45 (@sim_buzsaki,0:0.0001:0.1,y0);
    figure; plot(T,Y)

end




function Xpr = sim_buzsaki (t,X)

    V = X(1);
    mKDR = X(2);
    
    Iinj=1;

    % Parameters:
    gKDR = [20];
    KDR_V1 = [27];
    KDR_d1 = [11.5];
    KDR_V2 = [10];
    KDR_d2 = [10];
    E_KDR = [-100];
    IC = [0];
    IC_noise = [0.01];

    % Functions:
    minf = @(V) 1./(1+exp((-V-KDR_V1)/KDR_d1));
    mtau= @(V)  .25+4.35*exp(-abs(V+KDR_V2)/KDR_d2);
    aM = @(V)  minf(V) ./ mtau(V);
    bM = @(V) (1-minf(V))./mtau(V);
    IKDR = @(V,m) gKDR.*m.^4.*(V-E_KDR);
    
    % ODEs:
    mKDR_pr = aM(V).*(1-mKDR)-bM(V).*mKDR;
    V_pr = -IKDR(V,mKDR)+0;
    

    Xpr = [V_pr, mKDR_pr]';
    
end

