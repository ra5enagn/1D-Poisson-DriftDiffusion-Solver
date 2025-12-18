function [output] = LU_Decomposition(a,b,c,f,n_max)

    d = zeros(n_max,1);
    v = zeros(n_max,1);
    
    d(1) = b(1);
    for i = 2: n_max
        d(i) = b(i) - a(i)*c(i-1)/d(i-1);
    end

    % Solution of Lv = f %

    v(1) = f(1);
    for i = 2: n_max
        v(i) = f(i) - a(i)*v(i-1)/d(i-1)  ;
    end

    % Solution of U*fi = v %

    output(n_max) = v(n_max)/d(n_max);
    for i = n_max-1:-1:1
        output(i) = (v(i)-c(i)*output(i+1))/d(i);
    end
    
end