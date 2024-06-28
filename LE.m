%% calculate the LLE (Fig 5 D)
% load an 1-dimensional data
signal = ictal; 

% set parameters
tau = 5; % time delay
dim = 2; % embed dimension
T = length(signal); % data length


LE = lyapunov_exponent2(signal, tau, dim, T);

% output
fprintf('LLE = %f\n', LE);


function [LE] = lyapunov_exponent2(signal, tau, m, T)
    N = length(signal); 
    X = zeros(N - (m-1)*tau, m); 
    
    for i = 1:m
        X(:,i) = signal((1+(i-1)*tau):(N-(m-i)*tau));
    end
    
    J = zeros(N - (m-1)*tau, 1);
    D = zeros(N - (m-1)*tau, 1);
    LE = 0;
    
    for i = 1:(N - (m-1)*tau)
        distances = sqrt(sum((X - X(i,:)).^2, 2));
        [~, nearest] = sort(distances);
        J(i) = nearest(2); 
        
        for t = 1:T
            if i+t <= N-(m-1)*tau && J(i)+t <= N-(m-1)*tau
                D(i) = D(i) + (signal(i+t) - signal(J(i)+t))^2;
            end
        end
        
        D(i) = sqrt(D(i)/T);
        
        if D(i) == 0 || isnan(D(i))
            fprintf('Warn! \n');
            continue;
        end
    end
    
    D(D==0 | isnan(D)) = []; % 移除零和NaN
    
    if isempty(D)
        fprintf('Cannot get Lyapunov exponent. \n');
        LE = NaN;
        return;
    end
    
    LE = mean(log(D));
    
    if isnan(LE)
        fprintf('Lyapunov exponent = NaN. \n');
    end
end