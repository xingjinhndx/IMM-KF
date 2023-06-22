clear
close all
clc
r_m_std_ALL = logspace(-1, 1, 40);
r_dm_std_ALL = logspace(-1, 1, 40);
draw_imm_result = 1;
% determine tha parameter 
nnnn = 14;
mmmm = 33;
r_m_std = r_m_std_ALL(nnnn);
r_dm_std = r_dm_std_ALL(mmmm);
r_m_std=0.46416;
r_dm_std = 4.3755;
plot_versus_std(nnnn,mmmm,r_m_std,r_dm_std,draw_imm_result);


function []=plot_versus_std(nnnn,mmmm,r_m_std,r_dm_std,draw_imm_result)
%
close all;
T=0.1;
load("data_forfilter_row187_col51.mat")
data_line_CY = data_forfilter_row187_col51.data_line_CY;
data_line_FY = data_forfilter_row187_col51.data_line_FY;
data_line_real = data_forfilter_row187_col51.data_line_real;

Z_m = data_line_CY;

Z_dm = data_line_FY;
X = data_line_real;
X00 = [X(1); 0; 0];
N = size(X, 2);

%% initial
RMSE_ca_cv=zeros(3,N);
RMSE_measurement_m=zeros(1,N);
RMSE_measurement_dm=zeros(1,N);

% calculate the noise
RMSE_measurement_m  = abs(X(1,:) - Z_m).^2;
RMSE_measurement_dm  = abs(X(1,:) - Z_dm).^2;
H=[1 0 0];
u_input=[1/4 1/4 1/4 1/4]';
temp_p=diag([100^2 10^2 10]);
model_number=4;% model num
x1=[];
p1=[];
for i=1:model_number
    x1=[x1 X00];
    p1=[p1 temp_p];
end
u_output_ca_cv=zeros(model_number,1);

for m=2:N
 
    [X_estimate_ca_cv(:,m),P_ca_cv(:,:,m),x_model_filter_ca_cv,p_model_filter_ca_cv,u_output_ca_cv(:,m)]=imm_ca_cv2(Z_m(:,m),Z_dm(:,m),x1,p1,u_input, T,r_m_std,r_dm_std);
    x1=x_model_filter_ca_cv;
    p1=p_model_filter_ca_cv;
    u_input=u_output_ca_cv(:,m);
    %
    RMSE_ca_cv(:,m)=RMSE_ca_cv(:,m)+(X_estimate_ca_cv(:,m)-X(:,m)).^2;
end



%% 

RMSE_ca_cv=sqrt(RMSE_ca_cv);
RMSE_measurement_m=sqrt(RMSE_measurement_m);
RMSE_measurement_dm=sqrt(RMSE_measurement_dm);

line_w = 2.2;
mark_s = 7.5;
%color used in this paper
color_0 = [242/255 204/255 142/255];
color_1 = [130/255 178/255 154/255];
color_2 = [255/255 183/255 3/255];
color_3 = [191/255 30/255 46/255];
color_4 = [14/255 96/255 107/255]; 
if draw_imm_result
figure;
p =plot(Z_m,'--o','color',[130/255 178/255 154/255],'linewidth',line_w,'MarkerSize',mark_s-2,'MarkerFaceColor',[130/255 178/255 154/255]);
hold on
plot(Z_dm,'--o','color',color_4,'linewidth',line_w,'MarkerSize',mark_s-2,'MarkerFaceColor',color_4);
hold on
X_estimate_ca_cv(:,1)=Z_m(:,1); 
plot(X_estimate_ca_cv(1,:),'-s','color',color_3,'linewidth',line_w,'MarkerSize',mark_s,'MarkerFaceColor',color_3);
hold on
plot(X(1,:),'-s','color',color_2,'linewidth',line_w,'MarkerSize',mark_s, 'MarkerFaceColor',color_2);
hold on;
legend('CYGNSS SSH','FY-3E SSH','IMM-KF results','Insitu SSH');
xlabel('Date (day)');
ylabel('MSS in the specific grid (m)');
grid on;
hold off;
xticks([1 13 25 37 49])
xticklabels({'August 1st','August 5th','August 9th','August 13th','August 17th'})


xlim([0,52]);
ylim([63.5,73.5]);
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.5);
set(gcf,'Units','centimeter','Position',[1 1 45 12]);

end



u_output_ca_cv(:, 1) = [1/4 1/4 1/4 1/4]';

figure;
color_p_1 = [223 74 104]/255;
color_p_2 = [230 111 81]/255;
color_p_3 = [42 157 140]/255;
color_p_4 = [102 155 187]/255;


subplot(4,1,1)
plot(u_output_ca_cv(1,:),'-s','color',color_p_1,'LineWidth',line_w);
xticks([1 13 25 37 49])
xticklabels({'August 1st','August 5th','August 9th','August 13th','August 17th'})
xlim([0,52]);
ylim([0 1]);
legend('CY-CA');
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',line_w);
grid minor

hold on;
subplot(4,1,2)
plot(u_output_ca_cv(2,:),'-s','color',color_p_2,'LineWidth',line_w);
xticks([1 13 25 37 49])
xticklabels({'August 1st','August 5th','August 9th','August 13th','August 17th'})
xlim([0,52]);
ylim([0 1]);

legend('CY-CV');
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.5);
grid minor

subplot(4,1,3)
plot(u_output_ca_cv(3,:),'-s','color',color_p_3,'LineWidth',line_w);
xticks([1 13 25 37 49])
xticklabels({'August 1st','August 5th','August 9th','August 13th','August 17th'})
xlim([0,52]);
ylim([0 1]);
legend('FY-CA');

set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',line_w);
grid minor

subplot(4,1,4)
plot(u_output_ca_cv(4,:),'-s','color',color_p_4,'LineWidth',line_w);
xlim([0,52]);
ylim([0 1]);
legend('FY-CV');
xlabel('Date (day)');
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',line_w);
grid minor
xticks([1 13 25 37 49])
xticklabels({'August 1st','August 5th','August 9th','August 13th','August 17th'})
xlim([0,52]);
ylim([0 1]);
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.5);



set(gcf,'Units','centimeter','Position',[1 5 45 20]);
ylabel('Probability of the 4 models','position',[-2.1,2.8]);
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.5);
grid;

%%%RMSE
RMSE_m = sqrt(mean((Z_m - X).^2))
RMSE_dm = sqrt(mean((Z_dm - X).^2))
RMSE_filtered = sqrt(mean((X_estimate_ca_cv(1,:) - X).^2))
end
%%











