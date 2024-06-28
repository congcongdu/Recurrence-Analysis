%% for calculating the chaos degree(CD)
%% get a long time serie
% ll = 1000000;
% x = normrnd(0,1,[1,ll]);
%AA=A;
mm = 300;
x = AA(mm,:);

%% use FFT to get the period
l = 20000; %length
p = 1; %start point

dis = zeros(1,l);
O = zeros(1,l);

    dt = 1;
    X = x(p:dt:p+(l-1));
%     Y = y(p:dt:p+(l-1));
%     Z = z(p:dt:p+(l-1));
    for j = 1:l/dt   
        dis(j) = abs(X(1)-X(j));%sqrt((X(1)-X(j))^2+(Y(1)-Y(j))^2+(Z(1)-Z(j))^2);
    end
    O = dis;
    dt = 1;%h;
    Fs = 1/dt;
    t = 1:l;
    xn = O; 
    frei=40;
    NN=l;  
    xn=xn-mean(xn);
    XN=fft(xn,NN*frei)/NN/frei; 
    f0=1/(dt*NN*frei);  
    f=[0:ceil((NN*frei-1)/2)]*f0;  
    A=abs(XN); 
    B=2*A(1:ceil((NN*frei-1)/2)+1);
    aa=plot(f(2:end),B(2:end),'r','LineWidth',2)
    xlim([0,0.005])
aa.Color(4) = 0.4;
%set(gca,'yscale','log')

[max_num,max_idx] = max(B);
cd = f(max_idx);

%% calculate SDR
cd = 1/cd;
d1 = disp(AA);


%% calculate CD

m=mapminmax(m);
%chen 0.6036T rossler 5.8565T lorenz 0.714T
m2=m(:,1:121);
mm = repmat(m2,1,100);
Distances=abs(mm-m);
