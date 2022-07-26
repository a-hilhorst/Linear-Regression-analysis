function [pyu, pyl, p, y, rsq_adj, rpearson, resid, tresid, ser] = lin_reg_2d(ewf1,yo,x,alpha)
    p = polyfit(ewf1,yo,1);
    y = polyval(p,x);    
    yf = polyval(p, ewf1);
    N = length(yo);
    resid = yo-yf;
    SSresid = sum(resid.^2)/(N-2);
    ser = sqrt(SSresid);
    SStot = (N-1)*var(yo);
%     rsq = 1 - (N-2)*SSresid/SStot;
    rsq_adj = 1 - (N-2)^2*SSresid/SStot/(N-1); % is it the correct formula?
%     rsq_adj = 1 - (N-1)*SSresid/SStot;
    tresid = resid./(sqrt(SSresid*(1-(1/N+(ewf1-mean(ewf1)).^2/sum((ewf1-mean(ewf1)).^2)))));
    rpearson = (sum(ewf1.*yo)-N*mean(ewf1)*mean(yo))/(sqrt(sum((ewf1-mean(ewf1)).^2))*sqrt(sum((yo-mean(yo)).^2)));
    
    % confidence intervals
    pyu = p(2)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/N+((x-mean(ewf1)).^2/sum((ewf1-mean(ewf1)).^2))))+...
        (p(1)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((ewf1-mean(ewf1)).^2))))*x;
    pyl = p(2)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/N+((x-mean(ewf1)).^2/sum((ewf1-mean(ewf1)).^2))))+...
        (p(1)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((ewf1-mean(ewf1)).^2))))*x;
    
    % prediction intervals
%     pyu = p(2)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1+1/N+((x-mean(ewf1)).^2/sum((ewf1-mean(ewf1)).^2))))+...
%         (p(1)+tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((ewf1-mean(ewf1)).^2))))*x;
%     pyl = p(2)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1+1/N+((x-mean(ewf1)).^2/sum((ewf1-mean(ewf1)).^2))))+...
%         (p(1)-tinv(1-alpha/200,N-2)*sqrt(SSresid*(1/sum((ewf1-mean(ewf1)).^2))))*x;
end
