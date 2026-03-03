clc; clear; close all

% Fixed parameters
Nc = 1000*5;
nu = 0.001;
log_tol = -2;

% Range of sig_x scaling factors
%factor_values = logspace(log10(0.1), log10(0.001), 50);
factor_values = logspace(log10(0.1), log10(0.0001), 50);
logJ_list = [];
logFactor_list = [];

for i = 1:length(factor_values)
    factor = factor_values(i);
    J = PIELM_solver(Nc, nu, factor);
    logJ = log10(J);

    logJ_list(end+1) = logJ;
    logFactor_list(end+1) = log10(factor);

    if logJ < log_tol
        fprintf('Breaking at factor = %.4g because log10(J) = %.4f\n', factor, logJ);
        break;
    end
end

% Save data to CSV (both columns together)
data_to_save = [logFactor_list(:), logJ_list(:)];
filename = sprintf('logFactor_logJ_Nc%d.csv', Nc);
writematrix(data_to_save, filename);

% Plot
f3 = figure;
plot(logFactor_list, logJ_list, '-ko', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('$$\log_{10}(\eta)$$','Interpreter','latex','FontSize',24);
ylabel('$$\log_{10}(J)$$','Interpreter','latex','FontSize',24);
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');

%----------------------------------------
function J = PIELM_solver(Nc, nu, factor)
% Solves: u_x - nu u_xx = 0,  x in [0,1], u(0)=0, u(1)=1   %%% CHANGED

xL = 0; xR = 1;

% Collocation Points
X_pde = linspace(xL, xR, Nc)';
X_bc  = [xL; xR];

% Gaussian RBFs
NN = round(0.75 * Nc);
alpha_star = linspace(xL, xR, NN)';
sig_x = factor * ones(size(alpha_star));

% Parameters for RBFs
m = 1 ./ (sqrt(2) * sig_x);
alpha = -m .* alpha_star;

% PDE Residuals
LHS_PDE = zeros(Nc, NN);
RHS_PDE = zeros(Nc, 1);
for k = 1:Nc
    x = X_pde(k);

    z_sqr = (m' * x + alpha').^2;
    Phi   = exp(-z_sqr);

    % Enforce: u_x - nu u_xx = 0  %%% CHANGED
    % phi_x  = -2*Phi .* (m' .* (m'*x + alpha'))
    % phi_xx =  2*Phi .* (m'.^2) .* (1 - 2*z_sqr)
    LHS_PDE(k,:) = -2 * Phi .* (m' .* (m' * x + alpha')) ...
                   - 2 * nu * Phi .* (m'.^2) .* (1 - 2 * z_sqr);  %%% CHANGED
end

% Boundary Conditions
LHS_BC = zeros(2, NN);
RHS_BC = [0; 1];
for k = 1:2
    x = X_bc(k);
    z_sqr = (m' * x + alpha').^2;
    LHS_BC(k,:) = exp(-z_sqr);
end

% System Assembly
H = [LHS_PDE; LHS_BC];
b = [RHS_PDE; RHS_BC];

% Solve (LS)
c = lsqminnorm(H, b);   % more stable than pinv (recommended)

% Residual
J = norm(H * c - b, Inf);

% Compare (now consistent with u_x - nu u_xx = 0)
u_exact = EXACT_SOLUTION(X_pde, nu);
u_pielm = gauss_rbf_eval(X_pde, m, alpha, c);

% Plot PIELM vs Exact
f1 = figure(1); clf
plot(X_pde, u_pielm, '-r', X_pde, u_exact, '--b', 'LineWidth', 3);
axis([0 1 0 1]);
xlabel('$$x$$','Interpreter','latex','FontSize',30);
ylabel('$$u$$','Interpreter','latex','FontSize',30);
title('PIELM Solution (fixed kernel)','Interpreter','latex','FontSize',20);
legend('PIELM','Exact','Interpreter','latex','FontSize',15,'Location','best');
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');

end

%----------------------------------------
function u = gauss_rbf_eval(X, m, alpha, c)
% u(x) = sum_i c_i exp(-(m_i x + alpha_i)^2)
X = X(:);
u = zeros(size(X));
for k = 1:numel(X)
    z_sqr = (m' * X(k) + alpha').^2;
    u(k) = exp(-z_sqr) * c;
end
end

%----------------------------------------
function u_exact = EXACT_SOLUTION(X, nu)
% Exact for: u_x - nu u_xx = 0, u(0)=0, u(1)=1
% u(x) = (e^{x/nu} - 1) / (e^{1/nu} - 1)

X = X(:);
overflow_threshold = 1 / log(realmax);

if nu > overflow_threshold
    u_exact = expm1(X ./ nu) ./ expm1(1 / nu);
else
    exponent = (X - 1) ./ nu;
    threshold = -log(eps(class(X)));
    u_exact = exp(exponent);
    u_exact(exponent < -threshold) = 0;
    u_exact(X == 1) = 1;
end

if any(abs(X(:)-1) < eps(class(X)) & X(:) ~= 1)
    u_exact(abs(X-1) < eps(class(X))) = 1;
end
end