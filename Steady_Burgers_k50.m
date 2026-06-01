%{
Two-level KAPI-style Bayesian optimization for RBF allocation
in steady Burgers/X-TFC.

Sharp-gradient version with k = 50.

Problem:
    -nu*u''(x) + u(x)*u'(x) = f(x),   x in [0, pi]

Boundary conditions:
    u(0)  = 0
    u(pi) = 0

Manufactured exact solution:
    u_exact(x) = sin(x) * [1 + A*tanh(k*(x - xc))]

Two-level KAPI idea:
    Instead of one adaptive RBF cloud, use two adaptive components:

        1. Local component:
           resolves the sharp transition.

        2. Wide component:
           covers the surrounding smooth region.

    The optimized KAPI vector is

        w = [mu_loc, tau_loc, sigma_loc, ...
             mu_wide, tau_wide, sigma_wide]

Figures produced:
    1. F_Burgers_KAPI_BayesOpt_History.png
    2. F_Burgers_KAPI_RBF_Distribution.png
    3. F_Burgers_KAPI_Solution_Error.png

Requires:
    Statistics and Machine Learning Toolbox for bayesopt.
%}

clear; close all; clc;
format long
rng(42,'twister');

%% ------------------------------------------------------------------------
%  Problem parameters
% -------------------------------------------------------------------------

x0 = 0;
xf = pi;

z0 = 0;
zf = 1;

c = (zf - z0)/(xf - x0);
c2 = c^2;

% Burgers viscosity
nu_burg = 0.02;

% Sharp-gradient exact solution parameters
A_sharp = 0.35;
k_sharp = 50;
xc = pi/2;

% Boundary values
u0 = 0;
uf = 0;

% True sharp-layer location in normalized z-coordinate
zc_true = (xc - x0)/(xf - x0);
layer_width_z = 1/(k_sharp*pi);

fprintf('\n============================================================\n');
fprintf('Two-level KAPI BayesOpt Burgers/X-TFC solver: k = 50\n');
fprintf('============================================================\n');
fprintf('True layer location z_c        = %.6f\n', zc_true);
fprintf('Estimated layer width in z     = %.6e\n', layer_width_z);
fprintf('Burgers viscosity nu           = %.6e\n', nu_burg);
fprintf('Sharpness parameter k          = %.6f\n', k_sharp);

%% ------------------------------------------------------------------------
%  Numerical resolution
% -------------------------------------------------------------------------

% Important:
% Do not only increase RBF count. For k = 50, the main change is the
% two-level distribution. The following sizes are moderate but sufficient
% for testing the idea.
N_train_bg = 380;
N_train_cluster = 500;

N_val_bg = 900;
N_val_cluster = 450;

N_test = 2500;

L_bg = 160;
L_loc = 240;
L_wide = 280;

% Inner nonlinear solver options
inner_opts.IterMax = 80;
inner_opts.IterTol = 1e-10;

% Conservative nonlinear update for sharper transition
inner_opts.alpha0 = 0.5;
inner_opts.alpha_min = 1e-5;

% Stabilization for narrow RBFs
inner_opts.ridge = 1e-8;

inner_opts.verbose = false;

fprintf('Training points background      = %d\n', N_train_bg);
fprintf('Training points clustered       = %d\n', N_train_cluster);
fprintf('Validation points background    = %d\n', N_val_bg);
fprintf('Validation points clustered     = %d\n', N_val_cluster);
fprintf('Test points                     = %d\n', N_test);
fprintf('Background RBF centers          = %d\n', L_bg);
fprintf('Local adaptive RBF centers      = %d\n', L_loc);
fprintf('Wide adaptive RBF centers       = %d\n', L_wide);

%% ------------------------------------------------------------------------
%  Training/collocation points
% -------------------------------------------------------------------------

% Background Chebyshev-like points
theta = linspace(0, pi, N_train_bg)';
z_bg = 0.5*(1 - cos(theta));

% Local collocation enrichment near the sharp layer.
% For k = 50, the layer width is about 6.37e-3 in z.
% We use a wider enrichment window but keep many points inside it.
scl = linspace(-1, 1, N_train_cluster)';
cluster_radius = 10*layer_width_z;
z_cluster = zc_true + cluster_radius*tanh(2.5*scl)/tanh(2.5);
z_cluster = min(max(z_cluster, z0), zf);

z_train = unique([z_bg; z_cluster; z0; zf]);
z_train = sort(z_train);
N_train = length(z_train);

x_train = x0 + (1/c)*(z_train - z0);

[u_exact_train, u_exact_x_train, u_exact_xx_train] = manufactured_solution( ...
    x_train, A_sharp, k_sharp, xc);

f_train = -nu_burg*u_exact_xx_train + u_exact_train.*u_exact_x_train;

fprintf('Actual number of training points = %d\n', N_train);

%% ------------------------------------------------------------------------
%  Validation points: mixed global + local residual validation
% -------------------------------------------------------------------------

% Global validation points
z_val_bg = linspace(z0, zf, N_val_bg)';

% Local validation points near transition
scl_val = linspace(-1, 1, N_val_cluster)';
val_cluster_radius = 10*layer_width_z;
z_val_cluster = zc_true + val_cluster_radius*tanh(2.5*scl_val)/tanh(2.5);
z_val_cluster = min(max(z_val_cluster, z0), zf);

z_val = unique([z_val_bg; z_val_cluster; z0; zf]);
z_val = sort(z_val);
N_val = length(z_val);

x_val = x0 + (1/c)*(z_val - z0);

[u_exact_val, u_exact_x_val, u_exact_xx_val] = manufactured_solution( ...
    x_val, A_sharp, k_sharp, xc);

f_val = -nu_burg*u_exact_xx_val + u_exact_val.*u_exact_x_val;

fprintf('Actual number of validation points = %d\n', N_val);

%% ------------------------------------------------------------------------
%  Test points
% -------------------------------------------------------------------------

z_test = linspace(z0, zf, N_test)';
x_test = x0 + (1/c)*(z_test - z0);

[u_exact_test, u_exact_x_test, u_exact_xx_test] = manufactured_solution( ...
    x_test, A_sharp, k_sharp, xc);

f_test = -nu_burg*u_exact_xx_test + u_exact_test.*u_exact_x_test;

%% ------------------------------------------------------------------------
%  Fixed background RBF distribution
% -------------------------------------------------------------------------

% Background RBF centers and widths provide global coverage.
centers_bg = linspace(z0, zf, L_bg)';

% Background width is not meant to resolve the layer.
sigma_bg = 0.035*ones(L_bg,1);

% Fixed random numbers for two adaptive distributions.
% Keeping these fixed makes J(w) deterministic.
r_center_loc = randn(L_loc,1);
r_width_loc = 0.75 + 0.50*rand(L_loc,1);

r_center_wide = randn(L_wide,1);
r_width_wide = 0.75 + 0.50*rand(L_wide,1);

%% ------------------------------------------------------------------------
%  Bayesian optimization variables
% -------------------------------------------------------------------------

% Local component:
%   near the sharp transition, with small spread and narrow widths.
%
% Wide component:
%   broader enrichment around the transition region and smooth branches.

vars = [
    optimizableVariable('mu_loc', [0.45, 0.55])
    optimizableVariable('tau_loc', [0.002, 0.030], 'Transform', 'log')
    optimizableVariable('sigma_loc', [0.001, 0.012], 'Transform', 'log')

    optimizableVariable('mu_wide', [0.35, 0.65])
    optimizableVariable('tau_wide', [0.040, 0.180], 'Transform', 'log')
    optimizableVariable('sigma_wide', [0.008, 0.040], 'Transform', 'log')
];

fprintf('\nBayesOpt search bounds:\n');
fprintf('mu_loc       in [0.4500, 0.5500]\n');
fprintf('tau_loc      in [0.0020, 0.0300] log-transform\n');
fprintf('sigma_loc    in [0.0010, 0.0120] log-transform\n');
fprintf('mu_wide      in [0.3500, 0.6500]\n');
fprintf('tau_wide     in [0.0400, 0.1800] log-transform\n');
fprintf('sigma_wide   in [0.0080, 0.0400] log-transform\n');

%% ------------------------------------------------------------------------
%  Objective function handle
% -------------------------------------------------------------------------

objFcn = @(T) bayes_objective_wrapper(T, ...
    z_train, f_train, ...
    z_val, f_val, ...
    z_test, u_exact_test, f_test, ...
    centers_bg, sigma_bg, ...
    r_center_loc, r_width_loc, ...
    r_center_wide, r_width_wide, ...
    nu_burg, u0, uf, c, c2, z0, zf, inner_opts);

%% ------------------------------------------------------------------------
%  Run Bayesian optimization
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('Starting Bayesian optimization\n');
fprintf('============================================================\n');

time_bayes = tic;

results = bayesopt(objFcn, vars, ...
    'MaxObjectiveEvaluations', 75, ...
    'NumSeedPoints', 20, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'IsObjectiveDeterministic', true, ...
    'UseParallel', false, ...
    'Verbose', 1, ...
    'PlotFcn', []);

elapsed_bayes = toc(time_bayes);

fprintf('\nBayesian optimization finished in %.4f s\n', elapsed_bayes);

%% ------------------------------------------------------------------------
%  Extract best parameters
% -------------------------------------------------------------------------

bestX = results.XAtMinObjective;

best_w = [
    bestX.mu_loc, bestX.tau_loc, bestX.sigma_loc, ...
    bestX.mu_wide, bestX.tau_wide, bestX.sigma_wide
];

fprintf('\n============================================================\n');
fprintf('Best two-level KAPI distribution from BayesOpt\n');
fprintf('============================================================\n');
fprintf('Best objective log10(Jval) = %.8f\n', results.MinObjective);
fprintf('mu_loc                    = %.8f\n', best_w(1));
fprintf('tau_loc                   = %.8f\n', best_w(2));
fprintf('sigma_loc                 = %.8f\n', best_w(3));
fprintf('mu_wide                   = %.8f\n', best_w(4));
fprintf('tau_wide                  = %.8f\n', best_w(5));
fprintf('sigma_wide                = %.8f\n', best_w(6));
fprintf('True layer location zc    = %.8f\n', zc_true);
fprintf('Estimated layer width z   = %.8e\n', layer_width_z);
fprintf('============================================================\n');

%% ------------------------------------------------------------------------
%  Recompute final solution using best parameters
% -------------------------------------------------------------------------

[Jlog_best, info_best] = KAPI_objective_bayes( ...
    best_w, ...
    z_train, f_train, ...
    z_val, f_val, ...
    z_test, u_exact_test, f_test, ...
    centers_bg, sigma_bg, ...
    r_center_loc, r_width_loc, ...
    r_center_wide, r_width_wide, ...
    nu_burg, u0, uf, c, c2, z0, zf, inner_opts, true);

xi_best = info_best.xi;
centers_best = info_best.centers;
sigma_best = info_best.sigma;

[Ft, Fdt, Fddt, Ct, Cdt, Cddt] = build_xtfc_rbf_matrices( ...
    z_test, centers_best, sigma_best, u0, uf, c, c2, z0, zf);

u_pred_test = Ft*xi_best + Ct;
ux_pred_test = Fdt*xi_best + Cdt;
uxx_pred_test = Fddt*xi_best + Cddt;

R_test = -nu_burg*uxx_pred_test + u_pred_test.*ux_pred_test - f_test;
err_test = abs(u_pred_test - u_exact_test);

test_rmse = sqrt(mean((u_pred_test-u_exact_test).^2));
test_maxerr = max(err_test);
test_res_rmse = sqrt(mean(R_test.^2));

fprintf('\nFinal test-grid diagnostics:\n');
fprintf('Validation residual RMSE = %.8e\n', info_best.val_res_rmse);
fprintf('Solution RMSE            = %.8e\n', test_rmse);
fprintf('Max abs error            = %.8e\n', test_maxerr);
fprintf('Residual RMSE            = %.8e\n', test_res_rmse);
fprintf('Inner iterations         = %d\n', info_best.inner_iter);

%% ------------------------------------------------------------------------
%  Plots
% -------------------------------------------------------------------------

plot_BayesOpt_History_Burgers_TwoLevel( ...
    results, zc_true, best_w, ...
    'F_Burgers_KAPI_BayesOpt_History.png');

scatter_centers_vs_widths_burgers_twolevel( ...
    info_best.centers_bg, info_best.sigma_bg, ...
    info_best.centers_loc, info_best.sigma_loc, ...
    info_best.centers_wide, info_best.sigma_wide, ...
    best_w, zc_true, ...
    'F_Burgers_KAPI_RBF_Distribution.png');

plot_Burgers_vs_Exact_Overlay( ...
    x_test, u_exact_test, u_pred_test, nu_burg, ...
    'F_Burgers_KAPI_Solution_Error.png');

%% ========================================================================
%  Local functions
% ========================================================================

function objective_value = bayes_objective_wrapper(T, ...
    z_train, f_train, ...
    z_val, f_val, ...
    z_test, u_exact_test, f_test, ...
    centers_bg, sigma_bg, ...
    r_center_loc, r_width_loc, ...
    r_center_wide, r_width_wide, ...
    nu_burg, u0, uf, c, c2, z0, zf, inner_opts)

    persistent eval_counter_local

    if isempty(eval_counter_local)
        eval_counter_local = 0;
    end

    eval_counter_local = eval_counter_local + 1;

    w = [
        T.mu_loc, T.tau_loc, T.sigma_loc, ...
        T.mu_wide, T.tau_wide, T.sigma_wide
    ];

    [Jlog, info] = KAPI_objective_bayes( ...
        w, ...
        z_train, f_train, ...
        z_val, f_val, ...
        z_test, u_exact_test, f_test, ...
        centers_bg, sigma_bg, ...
        r_center_loc, r_width_loc, ...
        r_center_wide, r_width_wide, ...
        nu_burg, u0, uf, c, c2, z0, zf, inner_opts, false);

    fprintf(['Eval %3d | log10(Jval)= %+9.4f | Jval= %.3e | ', ...
             'test RMSE= %.3e | ', ...
             'muL= %.5f | tauL= %.5f | sigL= %.5f | ', ...
             'muW= %.5f | tauW= %.5f | sigW= %.5f | inner=%d\n'], ...
        eval_counter_local, Jlog, info.val_res_rmse, info.test_sol_rmse, ...
        info.wclip(1), info.wclip(2), info.wclip(3), ...
        info.wclip(4), info.wclip(5), info.wclip(6), ...
        info.inner_iter);

    objective_value = Jlog;
end

function [Jlog, info] = KAPI_objective_bayes( ...
    w, ...
    z_train, f_train, ...
    z_val, f_val, ...
    z_test, u_exact_test, f_test, ...
    centers_bg, sigma_bg, ...
    r_center_loc, r_width_loc, ...
    r_center_wide, r_width_wide, ...
    nu_burg, u0, uf, c, c2, z0, zf, inner_opts, verbose_final)

    [centers, sigma, wclip, comp] = KAPI_rbf_distribution_bayes_twolevel( ...
        w, centers_bg, sigma_bg, ...
        r_center_loc, r_width_loc, ...
        r_center_wide, r_width_wide);

    [F, Fd, Fdd, C, Cd, Cdd] = build_xtfc_rbf_matrices( ...
        z_train, centers, sigma, u0, uf, c, c2, z0, zf);

    [xi, inner_iter] = solve_inner_burgers_nonlinear_ls_bayes( ...
        F, Fd, Fdd, C, Cd, Cdd, f_train, nu_burg, inner_opts);

    [Fv, Fdv, Fddv, Cv, Cdv, Cddv] = build_xtfc_rbf_matrices( ...
        z_val, centers, sigma, u0, uf, c, c2, z0, zf);

    uv = Fv*xi + Cv;
    udv = Fdv*xi + Cdv;
    uddv = Fddv*xi + Cddv;

    Rv = -nu_burg*uddv + uv.*udv - f_val;
    val_res_rmse = sqrt(mean(Rv.^2));

    [Ft, ~, ~, Ct, ~, ~] = build_xtfc_rbf_matrices( ...
        z_test, centers, sigma, u0, uf, c, c2, z0, zf);

    utest = Ft*xi + Ct;

    test_sol_rmse = sqrt(mean((utest - u_exact_test).^2));
    test_max_err = max(abs(utest - u_exact_test));

    % Small regularization to discourage pathological distributions.
    reg = 1e-8*sum(wclip.^2);

    J = val_res_rmse + reg;
    Jlog = log10(J + 1e-16);

    if verbose_final
        fprintf('\nRecomputed final candidate:\n');
        fprintf('    log10(Jval) = %.8f\n', Jlog);
        fprintf('    Jval RMSE   = %.8e\n', val_res_rmse);
        fprintf('    test RMSE   = %.8e\n', test_sol_rmse);
        fprintf('    test maxerr = %.8e\n', test_max_err);
    end

    if nargout > 1
        info.wclip = wclip;
        info.val_res_rmse = val_res_rmse;
        info.test_sol_rmse = test_sol_rmse;
        info.test_max_err = test_max_err;
        info.inner_iter = inner_iter;
        info.xi = xi;
        info.centers = centers;
        info.sigma = sigma;

        info.centers_bg = comp.centers_bg;
        info.sigma_bg = comp.sigma_bg;
        info.centers_loc = comp.centers_loc;
        info.sigma_loc = comp.sigma_loc;
        info.centers_wide = comp.centers_wide;
        info.sigma_wide = comp.sigma_wide;
    end
end

function [centers, sigma, wclip, comp] = KAPI_rbf_distribution_bayes_twolevel( ...
    w, centers_bg, sigma_bg, ...
    r_center_loc, r_width_loc, ...
    r_center_wide, r_width_wide)

    % Clip bounds
    mu_loc = min(max(w(1), 0.45), 0.55);
    tau_loc = min(max(w(2), 0.002), 0.030);
    sigma_loc = min(max(w(3), 0.001), 0.012);

    mu_wide = min(max(w(4), 0.35), 0.65);
    tau_wide = min(max(w(5), 0.040), 0.180);
    sigma_wide = min(max(w(6), 0.008), 0.040);

    wclip = [
        mu_loc, tau_loc, sigma_loc, ...
        mu_wide, tau_wide, sigma_wide
    ];

    % Local adaptive component
    centers_loc = mu_loc + tau_loc*r_center_loc;
    centers_loc = reflect_to_unit_interval(centers_loc);

    sigma_loc_vec = sigma_loc*r_width_loc;
    sigma_loc_vec = max(0.0007, min(0.018, sigma_loc_vec));

    % Wide adaptive component
    centers_wide = mu_wide + tau_wide*r_center_wide;
    centers_wide = reflect_to_unit_interval(centers_wide);

    sigma_wide_vec = sigma_wide*r_width_wide;
    sigma_wide_vec = max(0.004, min(0.070, sigma_wide_vec));

    centers = [
        centers_bg(:);
        centers_loc(:);
        centers_wide(:)
    ];

    sigma = [
        sigma_bg(:);
        sigma_loc_vec(:);
        sigma_wide_vec(:)
    ];

    comp.centers_bg = centers_bg(:);
    comp.sigma_bg = sigma_bg(:);
    comp.centers_loc = centers_loc(:);
    comp.sigma_loc = sigma_loc_vec(:);
    comp.centers_wide = centers_wide(:);
    comp.sigma_wide = sigma_wide_vec(:);
end

function x = reflect_to_unit_interval(x)

    x = mod(x,2);
    idx = x > 1;
    x(idx) = 2 - x(idx);
    x = max(0,min(1,x));
end

function [F, Fd, Fdd, C, Cd, Cdd] = build_xtfc_rbf_matrices( ...
    z, centers, sigma, u0, uf, c, c2, z0, zf)

    z = z(:);
    n = length(z);

    dz = zf - z0;

    om1 = (zf - z)./dz;
    om2 = (z - z0)./dz;

    om1d = -ones(n,1)/dz;
    om2d =  ones(n,1)/dz;

    om1dd = zeros(n,1);
    om2dd = zeros(n,1);

    [H, Hz, Hzz] = eval_gaussian_rbf_basis(z, centers, sigma);

    h0 = H(1,:);
    hf = H(end,:);

    F_z = H - om1.*h0 - om2.*hf;
    Fd_z = Hz - om1d.*h0 - om2d.*hf;
    Fdd_z = Hzz - om1dd.*h0 - om2dd.*hf;

    F = F_z;
    Fd = c*Fd_z;
    Fdd = c2*Fdd_z;

    C = om1*u0 + om2*uf;
    Cd = c*(om1d*u0 + om2d*uf);
    Cdd = c2*(om1dd*u0 + om2dd*uf);
end

function [xi, iter] = solve_inner_burgers_nonlinear_ls_bayes( ...
    F, Fd, Fdd, C, Cd, Cdd, f, nu_burg, opts)

    L = size(F,2);
    xi = zeros(L,1);

    R = residual_burgers(xi, F, Fd, Fdd, C, Cd, Cdd, f, nu_burg);
    res_old = norm(R);

    for iter = 1:opts.IterMax

        u   = F*xi   + C;
        ud  = Fd*xi  + Cd;
        udd = Fdd*xi + Cdd;

        R = -nu_burg*udd + u.*ud - f;

        J = -nu_burg*Fdd + u.*Fd + ud.*F;

        ridge = opts.ridge;

        dxi = (J.'*J + ridge*eye(L))\(J.'*R);

        alpha = opts.alpha0;
        accepted = false;

        while alpha >= opts.alpha_min
            xi_trial = xi - alpha*dxi;

            R_trial = residual_burgers( ...
                xi_trial, F, Fd, Fdd, C, Cd, Cdd, f, nu_burg);

            if norm(R_trial) < norm(R)
                xi = xi_trial;
                accepted = true;
                break;
            end

            alpha = 0.5*alpha;
        end

        if ~accepted
            xi = xi - opts.alpha_min*dxi;
        end

        R_new = residual_burgers(xi, F, Fd, Fdd, C, Cd, Cdd, f, nu_burg);
        res_new = norm(R_new);

        if opts.verbose
            fprintf('    inner %3d | residual %.6e | step %.3e\n', ...
                iter, res_new, norm(dxi));
        end

        if res_new < opts.IterTol
            break;
        end

        if abs(res_old - res_new) < opts.IterTol*max(1,res_old)
            break;
        end

        res_old = res_new;
    end
end

function R = residual_burgers(xi, F, Fd, Fdd, C, Cd, Cdd, f, nu_burg)

    u   = F*xi   + C;
    ud  = Fd*xi  + Cd;
    udd = Fdd*xi + Cdd;

    R = -nu_burg*udd + u.*ud - f;
end

function [H, Hz, Hzz] = eval_gaussian_rbf_basis(z, centers, sigma)

    z = z(:);
    centers = centers(:).';
    sigma = sigma(:).';

    R = (z - centers)./sigma;

    H = exp(-R.^2);

    Hz = -(2./sigma).*R.*H;

    Hzz = ((4*R.^2 - 2)./(sigma.^2)).*H;
end

function [u, ux, uxx] = manufactured_solution(x, A, k, xc)

    q = tanh(k*(x - xc));
    r = 1 - q.^2;

    u = sin(x).*(1 + A*q);

    ux = cos(x).*(1 + A*q) ...
         + A*k*sin(x).*r;

    uxx = -sin(x).*(1 + A*q) ...
          + 2*A*k*cos(x).*r ...
          - 2*A*k^2*sin(x).*q.*r;
end

%% ========================================================================
%  Paper-style plotting functions
% ========================================================================

function plot_BayesOpt_History_Burgers_TwoLevel(results, zc_true, best_w, filename)

    loss_history = results.ObjectiveTrace;
    vars = results.XTrace;
    w_vals = table2array(vars);

    best_loss_history = zeros(size(loss_history));
    best_params = zeros(size(w_vals));

    min_val = Inf;
    min_idx = 1;

    for i = 1:length(loss_history)
        if loss_history(i) < min_val
            min_val = loss_history(i);
            min_idx = i;
        end

        best_loss_history(i) = min_val;
        best_params(i,:) = w_vals(min_idx,:);
    end

    fig = figure('Units', 'inches', ...
                 'Position', [1, 1, 8, 6], ...
                 'Color', 'w', ...
                 'Name', 'Two-Level KAPI BayesOpt History');

    subplot(2,1,1)
    plot(best_loss_history, 'b-', 'LineWidth', 2);
    xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$\log_{10}(J_{\mathrm{val}})$', ...
        'Interpreter', 'latex', 'FontSize', 20);
    title('Two-Level KAPI--ELM Bayesian Optimization History', ...
        'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'FontSize', 18, 'LineWidth', 1.2);
    grid on;
    box on;

    subplot(2,1,2)
    semilogy(best_params(:,1), 'r-', 'LineWidth', 2); hold on;
    semilogy(best_params(:,2), 'r--', 'LineWidth', 2);
    semilogy(best_params(:,3), 'r:', 'LineWidth', 2);

    semilogy(best_params(:,4), 'b-', 'LineWidth', 2);
    semilogy(best_params(:,5), 'b--', 'LineWidth', 2);
    semilogy(best_params(:,6), 'b:', 'LineWidth', 2);

    semilogy(zc_true*ones(size(best_params(:,1))), ...
        'k--', 'LineWidth', 1.5);

    xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Parameter Values', 'Interpreter', 'latex', 'FontSize', 20);

    legend({ ...
        '$\mu_{\mathrm{loc}}$', '$\tau_{\mathrm{loc}}$', '$\sigma_{\mathrm{loc}}$', ...
        '$\mu_{\mathrm{wide}}$', '$\tau_{\mathrm{wide}}$', '$\sigma_{\mathrm{wide}}$', ...
        '$z_c$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 12, ...
        'Location', 'eastoutside');

    set(gca, 'FontSize', 16, 'LineWidth', 1.2);
    grid on;
    box on;

    sgtitle(sprintf(['Best local: $\\mu=%.4f$, $\\tau=%.4f$, $\\sigma=%.4f$; ', ...
                     'wide: $\\mu=%.4f$, $\\tau=%.4f$, $\\sigma=%.4f$'], ...
        best_w(1), best_w(2), best_w(3), best_w(4), best_w(5), best_w(6)), ...
        'Interpreter', 'latex', 'FontSize', 13);

    if nargin < 4 || isempty(filename)
        filename = 'F_Burgers_KAPI_BayesOpt_History.png';
    end

    exportgraphics(fig, filename, 'Resolution', 300);
end

function scatter_centers_vs_widths_burgers_twolevel( ...
    centers_bg, sigma_bg, centers_loc, sigma_loc, centers_wide, sigma_wide, ...
    best_w, zc_true, save_path)

    fig = figure('Units', 'inches', ...
                 'Position', [1, 1, 8, 4], ...
                 'Color', 'w', ...
                 'Name', 'Two-Level RBF Centers vs Widths');

    scatter(centers_bg, sigma_bg, 8, [0.2 0.2 0.2], 'filled'); hold on;
    scatter(centers_loc, sigma_loc, 10, 'r', 'filled');
    scatter(centers_wide, sigma_wide, 10, 'b', 'filled');

    xline(zc_true, 'k--', 'LineWidth', 2);

    xlim([0, 1.05]);
    grid on;

    xlabel('$\alpha^\ast$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$\sigma$', 'Interpreter', 'latex', 'FontSize', 20);

    title('Two-Level Optimized RBF Distribution', ...
        'Interpreter', 'latex', 'FontSize', 18);

    legend({'Background', 'Local adaptive', 'Wide adaptive', '$z_c$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 14, ...
        'Location', 'northeast');

    set(gca, ...
        'FontSize', 20, ...
        'LineWidth', 1.2, ...
        'Box', 'off');

    if nargin > 8 && ~isempty(save_path)
        print(fig, save_path, '-dpng', '-r300');
        fprintf('Saved RBF center-width plot to: %s\n', save_path);
    end
end

function plot_Burgers_vs_Exact_Overlay(X_test, u_exact, u_kapi, nu, filename)

    X_test = X_test(:);
    u_exact = u_exact(:);
    u_kapi = u_kapi(:);

    abs_error = abs(u_exact - u_kapi);

    fig = figure('Units', 'inches', ...
                 'Position', [1, 1, 8, 4], ...
                 'Color', 'w', ...
                 'Name', 'Burgers Solution Comparison');

    yyaxis left
    hTrue = plot(X_test, u_exact, 'b-', 'LineWidth', 2); hold on;
    hPred = plot(X_test, u_kapi, 'r--', 'LineWidth', 2);

    ylabel_str = ['$u(x;\, \nu=' sprintf('%.2f', nu) ')$'];
    ylabel(ylabel_str, 'Interpreter', 'latex', 'FontSize', 20);

    ylim([min([u_exact; u_kapi])-0.1, ...
          max([u_exact; u_kapi])+0.1]);

    ax = gca;
    ax.YColor = 'k';

    yyaxis right
    hErr = area(X_test, abs_error, ...
        'FaceColor', [0.1 0.6 0.1], ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none');

    ylabel('$|u_{\mathrm{exact}}-u_{\mathrm{pred}}|$', ...
        'Interpreter', 'latex', 'FontSize', 20);

    ax = gca;
    ax.YColor = [0.1 0.6 0.1];

    if max(abs_error) > 0
        ylim([0, 1.05*max(abs_error)]);
    end

    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);

    legend([hErr, hPred, hTrue], ...
        {'Abs. Error', 'Predicted', 'True'}, ...
        'Location', 'northwest', ...
        'Interpreter', 'latex', ...
        'FontSize', 18);

    grid on;
    box on;
    set(gca, 'FontSize', 20, 'LineWidth', 1.2);

    if nargin < 5 || isempty(filename)
        filename = 'F_Burgers_KAPI_Solution_Error.png';
    end

    print(fig, filename, '-dpng', '-r300');
end