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
            fprintf('警告：出现零或NaN，跳过该点。\n');
            continue;
        end
    end
    
    D(D==0 | isnan(D)) = []; % 移除零和NaN
    
    if isempty(D)
        fprintf('无法计算Lyapunov exponent，所有D值都是零或NaN。\n');
        LE = NaN;
        return;
    end
    
    LE = mean(log(D));
    
    if isnan(LE)
        fprintf('Lyapunov exponent计算结果为NaN。\n');
    end
end
