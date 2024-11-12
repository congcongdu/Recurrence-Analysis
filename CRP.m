function [recur_pt, recur_dist, N_x, N_y] = CRP(x,y,varargin)
% Compute a cross recurrence plot.
% Authors: Hui Yang, Adam Meyers, and Mohammed Buqammaz
% Affiliation:
% The Pennsylvania State University
% 310 Leohard Building, University Park, PA
% Corresponding Email: yanghui@gmail.com
% If you find this demo useful, please cite the following paper:
% [1] H. Yang, "Multiscale Recurrence Quantification Analysis of Spatial Vectorcardiogram (VCG)
% Signals," IEEE Transactions on Biomedical Engineering, Vol. 58, No. 2, p339-347, 2011
% DOI: 10.1109/TBME.2010.2063704
% [2] A. Meyers, M. Buqammaz, and H. Yang, "Cross Recurrence Analysis for
% Pattern Matching of Multidimensional Physiological Signals," Chaos: An
% Interdisciplinary Journal of Nonlinear Science (Feature Article),
% Vol. 30, No. 12, p.123125, DOI: 10.1063/5.0030838
% Cross recurrence plot (CRP) is a method of bivariate timeseries analysis 
% capable of assessing the nonlinear interrela-tionships between time 
% series, and cross recurrence quan-tification analysis (CRQA) of the 
% CRP quantifies such in-terrelationships. 
% INPUTS
% (1) x
%     - 1st signal represented as an N_x X dim_x matrix, where N_x is the
%       length of signal x and dim_x is the dimension of x
% (2) y
%     - 2nd signal represented as an N_y X dim_y matrix, where N_y is the
%       length of signal y and dim_y is the dimension of y
% (3) threshold option
%     - optional name-value pairs to specify the threshold for calculating
%       CRP
%     - Name-value pairs:
%     (a) ('epsilon',epsilon)
%         - epsilon is a positive real number
%         - specifies exact threshold that should be used to calculate CRP
%     (b) ('RR',RR)
%         - RR is a real number between 0 and 100
%         - specifies that the threshold epsilon should be computed so as
%           to give a specified recurrence rate RR
%
% OUTPUTS
% - If nargout = 0 and nargin = 2, then unthresholded CRP is plotted.
% - If nargout = 0 and nargin > 2, then both unthresholded and threshold
%   CRPs are plotted.
% (1) recur_pt
%     - Nx2 matrix where N is # of recurrence points
%     - entry (i,1) gives row index and entry (i,2) gives column index of
%       each ith recurrence point in CRP
%     - signal y corresponds to rows & signal x corresponds to columns to
%       reflect that x is typically plotted on the x-axis and y on the
%       y-axis in the CRP
%     - columns are indexed starting from bottom to top to be consistent
%       with how they are presented visually in CRP; e.g., a recurrence
%       point occurring in the most bottom-left part of the CRP is indexed
%       as (1,1), not as (N_y,1)
[mx,nx] = size(x);
[my,ny] = size(y);
if mx>nx
    N_x = mx;
else
    x=x';
    N_x = nx;
end
if my>ny
    N_y = my;
else
    y=y';
    N_y = ny;
end
recur_dist = zeros(N_y,N_x);
for i=1:N_y
    for j=1:N_x
        recur_dist(i,j) = norm(y(i,:)-x(j,:));
    end
end
% Unthresholded recurrence plot
if nargout < 1
    figure('Position',[100 100 550 400]);
    imagesc(recur_dist);
    colormap Jet;
    colorbar;
    axis image;
    xlabel('Time Index','FontSize',10,'FontWeight','bold');
    ylabel('Time Index','FontSize',10,'FontWeight','bold');
    title('CRP','FontSize',10,'FontWeight','bold');
    get(gcf,'CurrentAxes');
    set(gca,'LineWidth',2,'FontSize',10,'FontWeight','bold','YDir','normal');
end
if nargin > 2
    
    if strcmp(varargin{1},'epsilon')
        epsilon = varargin{2};
    elseif strcmp(varargin{1},'RR')
        RR = varargin{2}/100;
        epsilon = quantile(recur_dist(:),RR);
    end
    
    [row,col] = find(recur_dist <= epsilon);
    recur_pt = [row col];
    
    % Thresholded recurrence plot
    if nargout < 1
        figure('Position',[100 100 550 400]);
        plot(row,col,'k.','MarkerSize',2);
        xlim([0,N_x]);
        ylim([0,N_y]);
        xlabel('Time Index','FontSize',10,'FontWeight','bold');
        ylabel('Time Index','FontSize',10,'FontWeight','bold');
        title('CRP','FontSize',10,'FontWeight','bold');
        get(gcf,'CurrentAxes');
        set(gca,'LineWidth',2,'FontSize',10,'FontWeight','bold');
    end
    
end
end
