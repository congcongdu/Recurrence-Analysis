%% Revised RP
%%get_rossler_data.m need to be run to get the time series: x,y,z
X = x;
x = X(:,1);
y = X(:,2);
z = X(:,3);
%%
st = 10200;    %start point
ll=8000;
l =8000; %time span
%ed = st+2*l;
m = zeros(ll,l);

%cut data
xx = x(st+1:end);
yy = y(st+1:end);
zz = z(st+1:end);

%calculate the distance
for i = 1:l
    for j = 1:ll
        m(j,i) = sqrt((xx(i+j)-xx(j))^2 + (yy(i+j)-yy(j))^2 + (zz(i+j)-zz(j))^2);
    end
end

%plot
figure(2)
imagesc(m)
colormap("jet")
colorbar;
%caxis([0,45]);
set(gca,'YDir','normal','FontSize',18); 

