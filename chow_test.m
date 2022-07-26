function [h, pval, outliers] = chow_test(x,y,start,palpha,h,pval,outliers)
%CHOW_TEST Chow statistical test.
% Performs a statistical test to identify the range over which the data
% provided follows a linear response.
%
% Author: Antoine Hilhorst
%
% See also LIN_REG_2D

    xtmp = x;
    ytmp = y;

    x(outliers)=[];
    y(outliers)=[];
    
    for i=length(ytmp)-start:-1:2
        p = polyfit(x(i:end),y(i:end),1);
        pstar = polyfit(x(i-1:end),y(i-1:end),1);

        yf = polyval(p, x(i:end));
        yfstar = polyval(pstar, x(i-1:end));

        SSR = sum((y(i:end)-yf).^2);
        SSRstar = sum((y(i-1:end)-yfstar).^2);
        
        n = length(x(i:end));
        k = 2;
        Fobs = (SSRstar-SSR)*(n-2*k)/SSR/k;
        
        h(i-1) = finv(1-palpha,1,n-2)>Fobs;
        pval(i-1) = 1-fcdf(Fobs,1,n-2);
        
        if ~h(i-1)            
            outliers = [outliers i-1];
            start = length(xtmp)-(i-1);
            [h, pval, outliers] = chow_test(xtmp,ytmp,start,palpha,h,pval,outliers);
            break
        end
    end
end