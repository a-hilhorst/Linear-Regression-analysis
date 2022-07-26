%%%
% Example
%
% Author: Antoine Hilhorst
% Date: 2022
%%%

close all
%% generate random values
x0 = 5;
x1 = 25;
x = linspace(x0,x1,100);

% linear model y = alpha1*x + beta1 + epsilon
% epsilon is white noise
alpha1 = 20;
beta1 = 250;
noise_ampl1 = 15;
y = alpha1*x + beta1 + noise_ampl1*randn(size(x));

% let's assume that y follows a different model below x0
xx = linspace(0,x0,10);
beta2 = 10;
alpha2 = alpha1+(beta1-beta2)/x0; % defined such that yy(x0) = y(x0)
noise_ampl2 = noise_ampl1;
yy = alpha2*xx + beta2 + noise_ampl2*randn(size(xx));

x = [xx(1:end-1) x];
y = [yy(1:end-1) y];

% let's add outliers
n_outliers = 3; % adds 1 to n_outliers outliers
outliers = round(length(x)*rand(n_outliers,1));
y(outliers) = mean(y) + std(y)*randn(size(outliers));

%% 
figure('Color','w','units','normalized','outerposition',[0 0 .5 1])
hold on

colors = lines(5);
plot(x,alpha1*x+beta1,'k-','Linewidth',1)
plot(x,y,'o','Color',colors(1,:),'Linewidth',2)

%% identifying x0 and outliers using Chow Test
start = 5;
palpha = 0.05;

h = zeros(size(x));
pval = zeros(size(x));
outliers_chow = [];

[h, pval, outliers_chow] = chow_test(x,y,start,palpha,h,pval,outliers_chow);

plot(x(outliers_chow),y(outliers_chow),'x','Color',colors(2,:),'Linewidth',2,'MarkerSize',10)

%% linear regression of the data without the outliers
% with confidence and prediction intervals

% removing the outliers
x(outliers_chow) = [];
y(outliers_chow) = [];

[p1, ~, ci, ~, val] = lin_reg_2d(x,y,x,0.05);
% plot(x,y_lr,'-','Color',colors(3,:),'Linewidth',2)
plot(x,ci(1,:),'--','Color',colors(3,:),'Linewidth',3)
plot(x,ci(2,:),'--','Color',colors(3,:),'Linewidth',3)
% plot(x,pi(1,:),'-.','Color',colors(3,:),'Linewidth',3)
% plot(x,pi(2,:),'-.','Color',colors(3,:),'Linewidth',3)

%% identifying outliers from Studentized residuals
outliers_tresid = (val.tresid>1.96) | (val.tresid<-1.96);

plot(x(outliers_tresid),y(outliers_tresid),'s','Color',colors(2,:),'Linewidth',2,'MarkerSize',10)

% removing the outliers
x(outliers_tresid) = [];
y(outliers_tresid) = [];

[p2, ~, ci, ~, ~] = lin_reg_2d(x,y,x,0.05);
% plot(x,y_lr,'-','Color',colors(4,:),'Linewidth',2)
plot(x,ci(1,:),'--','Color',colors(4,:),'Linewidth',3)
plot(x,ci(2,:),'--','Color',colors(4,:),'Linewidth',3)
% plot(x,pi(1,:),'-.','Color',colors(4,:),'Linewidth',3)
% plot(x,pi(2,:),'-.','Color',colors(4,:),'Linewidth',3)

%% identifying outliers from the CI
outliers_ci = (y>ci(1,:)) | (y<ci(2,:));

plot(x(outliers_ci),y(outliers_ci),'v','Color',colors(2,:),'Linewidth',2,'MarkerSize',10)

% removing the outliers
x(outliers_ci) = [];
y(outliers_ci) = [];

[p3, y_lr, ci, pi, ~] = lin_reg_2d(x,y,x,0.05);
plot(x,y_lr,'-','Color',colors(5,:),'Linewidth',2)
plot(x,ci(1,:),'--','Color',colors(5,:),'Linewidth',3)
plot(x,ci(2,:),'--','Color',colors(5,:),'Linewidth',3)
plot(x,pi(1,:),'-.','Color',colors(5,:),'Linewidth',3)
plot(x,pi(2,:),'-.','Color',colors(5,:),'Linewidth',3)

%% display linear regression estimators
disp([p1; ...
    p2; ...
    p3; ...
    alpha1 beta1])



