clc;
clear;

%% parameters
w = 1;
a = 0.2; %0.165
b = 0.2; %0.2
c = 5.7;		
x(1)=1;
y(1)=1;
z(1)=1;
h=0.01;
n = 3e4;
%% RK-4
for i=1:n
    K1=-(w*y(i)+z(i));
    L1=w*x(i)+a*y(i);
    M1=b+x(i)*z(i)-c*z(i);

    K2=-((y(i)+h/2*L1)+(z(i)+h/2*M1));
    L2=(x(i)+h/2*K1)+a*(y(i)+h/2*L1);
    M2=b+(x(i)+h/2*K1)*(z(i)+h/2*M1)-c*(z(i)+h/2*M1);

    K3=-((y(i)+h/2*L2)+(z(i)+h/2*M2));
    L3=(x(i)+h/2*K2)+a*(y(i)+h/2*L2);
    M3=b+(x(i)+h/2*K2)*(z(i)+h/2*M2)-c*(z(i)+h/2*M2); 

    K4=-((y(i)+h*L3)+(z(i)+h*M3));
    L4=(x(i)+h*K3)+a*(y(i)+h*L3);
    M4=b+(x(i)+h*K3)*(z(i)+h*M3)-c*(z(i)+h*M3);  

    x(i+1)=x(i)+h/6*(K1+2*K2+2*K3+K4);
    y(i+1)=y(i)+h/6*(L1+2*L2+2*L3+L4);
    z(i+1)=z(i)+h/6*(M1+2*M2+2*M3+M4);
end
clear K* L* M*
%% phase space
figure(1);
plot3(x(2001:3000),y(2001:3000),z(2001:3000));
xlabel('X(t)');
ylabel('Y(t)');
zlabel('Z(t)');
grid;

%% phase space with changing color

% 提取时间序列的维度
dimension1 = x(2001:3000);
dimension2 = y(2001:3000);
dimension3 = z(2001:3000);
% 定义颜色映射
cmap = hot(length(dimension1)+200); 

% 绘制相空间图
figure;
hold on;
for t = 1:length(dimension1)
    colorIndex = t ;
    scatter3(dimension1(t), dimension2(t),dimension3(t), [], cmap(colorIndex, :), 'filled');
end
hold off;

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('Z');
% 添加颜色刻度条
colormap(gca, cmap);
c = colorbar;
c.Label.String = 'Time';
view(-50,30)

%% classical RP(Fig 1)
dt = 1;
p = 10e5;
q = p + dt * 1e4;
X = [x(p:dt:q);y(p:dt:q);z(p:dt:q)];
M = X';

d1 = pdist(M);
d1 = squareform(d1);
figure('NumberTitle','off','Name','imagesc法','Color','w');
imagesc(d1);
colormap('jet');
colorbar;
set(gca,'YDir','normal'); 
