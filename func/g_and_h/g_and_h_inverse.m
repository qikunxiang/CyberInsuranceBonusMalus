function [z, iter] = g_and_h_inverse(x, g, h, tol)
%G_AND_H_INVERSE Inverse transform to recover the z values
% Inputs:
%       x: input to the cdf
%       g: the g parameter
%       h: the h parameter
%       tol: numerical tolerance, default is 1e-8
% Outputs:
%       z: corresponding z value for the standard normal distribution
%       iter: number of iterations

n = length(x);
z = ones(n, 1);

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-8;
end

if g ~= 0
    t_func = @(v)((exp(g * v) - 1) / g .* exp(0.5 * h * v .^ 2));

    d_func = @(v)((exp(g * v) .* (1 + h / g * v) - h / g * v) ...
        .* exp(0.5 * h * v .^ 2));
else
    t_func = @(v)(v .* exp(0.5 * h * v .^ 2));

    d_func = @(v)((1 + h * v .^ 2) .* exp(0.5 * h * v .^ 2));
end

diff0 = t_func(z) - x;
diff = diff0;
iter = 0;
exclude_list = false(n, 1);

while max(abs(diff(~exclude_list))) > tol
    z = z - diff ./ d_func(z);
    diff = t_func(z) - x;
    
    exclude_list = exclude_list | abs(diff) > max(1, abs(diff0)) ...
        | isinf(diff) | isnan(diff);
    iter = iter + 1;
end

if any(exclude_list)
    % the remaining divergent cases are solved first by bisection
    nr = sum(exclude_list);
    xr = x(exclude_list);
    
    z_range = 10 * repmat([-1, 1], nr, 1);
    f_range = t_func(z_range);
    
    while any(xr > f_range(:, 2) | xr < f_range(:, 1))
        z_range = z_range * 2;
        f_range = t_func(z_range);
        iter = iter + 1;
    end
    
    while max(f_range(:, 2) - f_range(:, 1)) > 1e-3
        z_mid = mean(z_range, 2);
        f_mid = t_func(z_mid);
        left_list = f_mid <= xr;
        right_list = f_mid > xr;
        z_range(left_list, 1) = z_mid(left_list);
        z_range(right_list, 2) = z_mid(right_list);
        f_range = t_func(z_range);
        iter = iter + 1;
    end
    
    zr = mean(z_range, 2);
    diff0 = t_func(zr) - xr;
    diff = diff0;
    while max(abs(diff)) > tol
        zr = zr - diff ./ d_func(zr);
        diff = t_func(zr) - xr;
        iter = iter + 1;
        if any(isnan(diff) | isinf(diff) | abs(diff) > abs(diff0) + 1)
            error('failed to converge');
        end
    end
    
    z(exclude_list) = zr;
end

assert(max(abs(t_func(z) - x)) < tol, 'unexpected error');

end

