clear; 
close all; 
clc; 
%% 加载数据以及数据预处理
load ('.\data\d00.mat')
load ('.\data\d05_te.mat')
d0=d00; 
d01=d05te; 
X_train=d0(:,[1:22,42:52]); 
Y_train=d0(:,35); 
X_test=d01(:,[1:22,42:52]); 
Y_test=d01(:,35); 
%% 数据标准化
[X_train, x_mean,x_v] = autos(X_train);
[Y_train, y_mean,y_v] = autos(Y_train);
X_test=(X_test-repmat(x_mean,size(X_test,1),1))./(repmat(x_v,size(X_test,1),1));
Y_test=(Y_test-repmat(y_mean,size(X_test,1),1))./(repmat(y_v,size(X_test,1),1));
%% pls建模
A = 3;
[T,P,Q,B,W]=pls_NIPALS(X_train,Y_train,A);
R=W*inv(P'*W); 
n1=size(X_train,1); % 训练数据的样本数
Ab=inv(T'*T/(n1-1));
n2=size(X_test,1); 
for i=1:n2 
    t_new=R'*X_test(i,:)';%故障数据 
    tsquare2(i)=t_new'*Ab*t_new; 
end
%% 确定控制限 Control Limit
alpha=0.95;%显著性水平
level1=A*(n1*n1-1)/(n1*(n1-A))*finv(alpha,A,(n1-A)); %T2控制限
 for i=1:n1 
     xr_old = X_train(i,:)'-P*R'*X_train(i,:)';%正常数据 
     SPE_old(i) = xr_old'*xr_old; 
 end
 S1=var(SPE_old); 
 mio1=mean(SPE_old);
 V1=2*mio1^2/S1; 
 level2=S1/(2*mio1)*chi2inv(alpha,V1); %SPE控制限
%% 监控结果可视化
t = 1:n2;
figure
set(gcf,'color','white')
subplot(2,1,1)
plot(t,tsquare2,t,ones(n2,1)*level1,'r--'); 
xlabel('Sample Number'); 
ylabel('T^2'); % title('PLS-fault data');
 for i=1:n2 
     xr_new=X_test(i,:)'-P*R'*X_test(i,:)';%故障数据 
     SPE(i)=xr_new'*xr_new; 
 end
legend('SPE','95% upper limit','Location','best')
subplot(2,1,2)
plot(t,SPE,t,ones(n2,1)*level2,'r--');
xlabel('Sample Number'); 
ylabel('SPE'); 
legend('Hotelling''s T^2','95% upper limit','Location','best')
%% 故障检测率 
% % T2 fault detection rate 
% t2_detect=0; 
% for i=161:n2 
%    if tsquare2(i)>=level1
        % t2_detect=t2_detect+1;
     % end
% end 
% t2_detect_rate=100*t2_detect/800
% % SPE fault detection rate
% spe_detect=0;
% for i=161:n2
    % if SPE2(i)>=level2 
    % spe_detect=spe_detect+1; 
    % end 
 % end 
 % spe_detect_rate=100*spe_detect/800 
 %% 误报率和漏报率 
 % 误报率计算 
 % T2 false alarm rate
%  t2_false=0; 
%  for i=1:160 
% 	if tsquare2(i)>level1 
% 		t2_false=t2_false+1; 
% 		end 
% end 
% t2_false_rate=100*t2_false/160 
% % SPE false alarm rate 
% spe_false=0; 
% for i=1:160 
% 	if SPE2(i)>=level2 
% 		spe_false=spe_false+1; 
% 	end 
% end 
% spe_false_rate=100*spe_false/160 
% % 漏报率计算 
% % T2 missed detection rate 
% t2_miss=0; 
% for i=161:n2 
% if tsquare2(i)<level1 
% 	t2_miss=t2_miss+1; 
% 	end 
% end 
% t2_miss_rate=100*t2_miss/800 
% % SPE missed detection rate 
% spe_miss=0; 
% for i=161:n2 
% 	if SPE2(i)<level2 
% 		spe_miss=spe_miss+1; 
% 	end 
% end 
% spe_miss_rate=100*spe_miss/800