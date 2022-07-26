function [p, y, ci, pi, values] = lin_reg_2d(x_data,y_data,x,alpha)
%LIN_REG_2D Simple linear regression.
% Computes the simple linear regression of (x_data,y_data), the confidence
% and prediction interval with p-value alpha over the range x, as well as
% useful values to analyse the regression.
%
% Author: Antoine Hilhorst
%
% See also CHOW_TEST

    p = polyfit(x_data,y_data,1);
    y = polyval(p,x);    
    yf = polyval(p, x_data);
    N = length(y_data);
    resid = y_data-yf;
    SSresid = sum(resid.^2)/(N-2);
    ser = sqrt(SSresid);
    SStot = (N-1)*var(y_data);
    rsq_adj = 1 - (N-2)^2*SSresid/SStot/(N-1); 
    tresid = resid./(sqrt(SSresid*(1-(1/N+(x_data-mean(x_data)).^2/sum((x_data-mean(x_data)).^2)))));
    rpearson = (sum(x_data.*y_data)-N*mean(x_data)*mean(y_data))/(sqrt(sum((x_data-mean(x_data)).^2))*sqrt(sum((y_data-mean(y_data)).^2)));
    
    % confidence intervals
    ci = [p(2)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/N+((x-mean(x_data)).^2/sum((x_data-mean(x_data)).^2))))+...
        (p(1)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((x_data-mean(x_data)).^2))))*x;
        p(2)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/N+((x-mean(x_data)).^2/sum((x_data-mean(x_data)).^2))))+...
        (p(1)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((x_data-mean(x_data)).^2))))*x];
    
    % prediction intervals
    pi = [p(2)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1+1/N+((x-mean(x_data)).^2/sum((x_data-mean(x_data)).^2))))+...
        (p(1)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((x_data-mean(x_data)).^2))))*x;
        p(2)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1+1/N+((x-mean(x_data)).^2/sum((x_data-mean(x_data)).^2))))+...
        (p(1)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((x_data-mean(x_data)).^2))))*x];

    %% output values of interest
    values.rsq_adj = rsq_adj;
    values.rpearson = rpearson;
    values.resid = resid;
    values.tresid = tresid;
    values.ser = ser;
end
