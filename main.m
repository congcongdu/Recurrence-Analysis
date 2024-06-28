%% step1 get matrix D and reverse
% Fig1 right
%%run GetRosslerData.m to get time series: x,y,z

st = 100000;    %start point
len = 5000;
l =5000;        %time span
ed = st+2*l;
m = zeros(len,l);

%cut data
xx = x(st+1:ed);
yy = y(st+1:ed);
zz = z(st+1:ed);

%calculate the distance
for i = 1:l
    for j = 1:l
        m(j,i) = sqrt((xx(i+j)-xx(j))^2 + (yy(i+j)-yy(j))^2 + (zz(i+j)-zz(j))^2);
    end
end


figure(1)
imagesc(m)
colormap("jet")
colorbar;
%caxis([0,45]);
set(gcf,'Color','White')
set(gca,'YDir','normal','FontSize',18); 
%% Epilepsy data
% Fig 4 ABC
xx = ictal;
for i = 1:512
    for j = 1:512
        m(j,i) = abs(xx(i+j)-xx(j));
    end
end


figure(2)
imagesc(m)
colormap("jet")
colorbar;
%caxis([0,45]);
set(gcf,'Color','White')
set(gca,'YDir','normal','FontSize',18); 
%% %%%%%% calculate all the mentioned index in our paper
rD = m;
%% 1. sum
lsum = sum(rD);
hsum = sum(rD');
figure(2)
set(gcf,'Color','White')
subplot 121
plot(lsum,'LineWidth',2)
xlabel('Time interval')
set(gca,'FontSize',18); 
subplot 122
plot(hsum,'LineWidth',2)
xlabel('Time')
set(gca,'FontSize',18); 

%% 2. RQA——RR(recurrence rate)
%[recur_pt,recur_dist,N_y1,N_y2] = CRP(y1,y2,'RR',5); 
% crqa_stat = CRQA(recur_pt,N_y1,N_y2)
e = 15;
y1 = [xx;yy;zz]; %classic system
%y1 = xx;       %Epilepsy data
[recur_pt,recur_dist,N_y1,N_y2] = CRP(y1,y1,'RR',e);
for i = 1:length(recur_dist)
    for j = 1:length(recur_dist)
        if recur_dist(i,j) <= e
            recur_dist(i,j) = 1;
        else
            recur_dist(i,j) = 0;
        end
    end
end
RR = sum(recur_dist)./(length(recur_dist) - 1);

len = length(hsum);

figure(3) % fig 2 B
plot(normalize(RR(1:len)),'LineWidth',2)
hold on
plot(normalize(hsum),'LineWidth',2)
xlim([0,len+1])
legend('RR','sum')
set(gcf,'Color','White')
set(gca,'FontSize',18); 

corrcoef(normalize(RR(1:len)),normalize(hsum)) 
%% Fig 2 C-D
figure(4)
set(gcf,'Color','White')
subplot 121
scatter3(xx(1:len),yy(1:len),zz(1:len),10,normalize(RR(1:len)),"filled")
colormap("jet")
colorbar
set(gca,'FontSize',18); 
subplot 122
scatter3(xx(1:len),yy(1:len),zz(1:len),10,normalize(hsum),"filled")
colormap("jet")
colorbar
set(gca,'FontSize',18); 
%% 3. FFT
%% get distance series
xx = x(st+1:end);
yy = y(st+1:end);
zz = z(st+1:end);

for i = 1:100000
    x1(i) = sqrt((xx(i)-xx(1))^2+(yy(i)-yy(1))^2+(yy(i)-yy(1))^2);
end
%% FFT: get PF and Amplitude(Fig 5 A-B)
xn = lsum; 
dt = h;% h
Fs = 1/dt;
t = 1:len;
frei=40;
NN=len;  
xn=xn-mean(xn);
XN=fft(xn,NN*frei)/NN/frei;  
f0=1/(dt*NN*frei);  
f=[0:ceil((NN*frei-1)/2)]*f0;  
A=abs(XN);  
B=2*A(1:ceil((NN*frei-1)/2)+1);
aa=plot(f(2:end),B(2:end),'LineWidth',2);
%xlim([0,0.005])

% aa.Color(4) = 0.4;
set(gcf,'Color','White')
set(gca,'FontSize',18); 
% set(gca,'yscale','log')
[val,idx]=max(B);
PF = f(idx);