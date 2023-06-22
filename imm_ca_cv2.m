%imm
function [X_estimate,P,x_model_filter,p_model_filter,u_output]=imm_ca_cv2(z_m,z_dm,x1,p1,u_input,T,r_m_std,r_dm_std)
%
%

u_output=zeros(4,1); %4
fai=[];
Q=[];


%
F_ca_m=[1 T T^2/2;0 1 T;0 0 1];
G_ca_m=[T^2/2,T,1]';
Q_ca_m=10;
Q_ca_m=G_ca_m*Q_ca_m*G_ca_m';
%
F_cv_m=[1 T 0;0 1 0;0 0 0];
    
G_cv_m=[T^2/2,T,0]';
Q_cv_m=10;
Q_cv_m=G_cv_m*Q_cv_m*G_cv_m';
%
F_ca_dm=[1 T T^2/2;0 1 T;0 0 1];
G_ca_dm=[T^2/2,T,1]';
Q_ca_dm=10;
Q_ca_dm=G_ca_dm*Q_ca_dm*G_ca_dm';
%
F_cv_dm=[1 T 0;0 1 0;0 0 0];
    
G_cv_dm=[T^2/2,T,0]';
Q_cv_dm=10;
Q_cv_dm=G_cv_dm*Q_cv_dm*G_cv_dm';


fai=[F_ca_m F_cv_m F_ca_dm F_cv_dm]; %变成4维度 [3*3 3*3 3*3 3*3]
Q=[Q_ca_m Q_cv_m Q_ca_dm Q_cv_dm];%变成4维度 [3*3 3*3 3*3 3*3]
h=[1 0 0];
% r=0.5^2; % zzb：可改！！！！！�?% 独立�?4份！！！！！！�?�米级用�?个， dm级的用一个�??
r_m = r_m_std^2;
r_dm = r_dm_std^2;
r = [r_m r_m r_dm r_dm];
%模型参数设置完毕

dim=3;%状�?�维�? x x' x''
model_number=4;%模型�? CV CA % �?4 = m dm × CV CA�?
% MarkoProb=[0.85 0.05 0.05 0.05;0.05 0.85 0.05 0.05;0.05 0.05 0.85 0.05;0.05 0.05 0.05 0.85]; 
MarkoProb=[0.95 0.05;0.05 0.95];%转移矩阵概率 4*4
MarkoProb = eye(4); %%% 应该要转移，暂定
predictProb=MarkoProb'*u_input;
mixedProb=[];
mixedInitX=zeros(dim,model_number);
mixedInitP=zeros(dim,dim*model_number); %3 *2 �? 3*4
mixedInitC=zeros(dim,dim*model_number); %3 *2 �? 3*4
for i=1:model_number
    for j=1:model_number
        mixedProb(i,j)=MarkoProb(i,j)*u_input(i)/predictProb(j); % 4 * 4
    end
end


mixedInitX=x1*mixedProb;


for j=1:model_number
    for i=1:model_number
        mixedInitP(:,(j-1)*dim+1:j*dim)=mixedInitP(:,(j-1)*dim+1:j*dim)+(p1(:,(i-1)*dim+1:i*dim)+(x1(:,i)-mixedInitX(:,j))*(x1(:,i)-mixedInitX(:,j))')*mixedProb(i,j);
    end
end
for i=1:model_number
    temp1=(i-1)*dim+1;
    temp2=i*dim;
    x_model_predict(:,i)=fai(:,temp1:temp2)*mixedInitX(:,i);
%   p_model_predict(:,temp1:temp2)=fai(:,temp1:temp2)*mixedInitP(:,temp1:temp2)*fai(:,temp1:temp2)'+g(:,(i-1)*2+1:i*2)*q(:,(i-1)*2+1:i*2)*g(:,(i-1)*2+1:i*2)';
   
    p_model_predict(:,temp1:temp2)=fai(:,temp1:temp2)*mixedInitP(:,temp1:temp2)*fai(:,temp1:temp2)'+Q(:,temp1:temp2);
    % 卡尔曼滤�?
    s(:,i)=h*p_model_predict(:,temp1:temp2)*h'+r(i);
    k(:,i)=p_model_predict(:,temp1:temp2)*h'*inv(s(:,i));
    if i == 1 || i == 2
        z = z_m;
    elseif i == 3 || i == 4
        z = z_dm;
    else
        clear z
    end
    v_model_filter(:,i)=z-h*x_model_predict(:,i);
    x_model_filter(:,i)=x_model_predict(:,i)+k(:,i)*v_model_filter(:,i);
    p_model_filter(:,temp1:temp2)=p_model_predict(:,temp1:temp2)-k(:,i)*s(:,i)*k(:,i)';
    % like function 极大似然函数�?
    likelyhood(i)=(det(2*pi*s(:,i)))^(-0.5)*exp(-0.5*v_model_filter(:,i)'*inv(s(:,i))*v_model_filter(:,i));
%     A=det(2*pi*s(:,(i-1)*2+1:i*2));
%     B=exp(-0.5*v_model_filter(:,i)'*inv(s(:,(i-1)*2+1:i*2))*v_model_filter(:,i));
    u_output(i)=likelyhood(i)*predictProb(i);

end

u_output=(u_output/sum(u_output));
X_estimate=x_model_filter*u_output;
% X_estimate=real(X_estimate);%%%%仿真表明估计可能出现虚部，�?�过三角分解可消除此现象，为�?略起见直接去掉虚�?
P=zeros(dim);
for i=1:model_number
    temp1=(i-1)*dim+1;
    temp2=i*dim;
    P=P+(p_model_filter(:,temp1:temp2)+(x_model_filter(:,i)-X_estimate)*(x_model_filter(:,i)-X_estimate)')*u_output(i);
end
   

    
    
    


