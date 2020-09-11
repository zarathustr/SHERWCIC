clear all
close all
clc

format long g
load data1

len = size(RA, 3);
dt = 1 / 20;
time = dt * (1 : len);
iter = 100;
A = zeros(4, 4, len);
B = zeros(4, 4, len);
for i = 1 : len
    A(1 : 3, 1 : 3, i) = RA(:, :, i);
    A(1 : 3, 4, i) = tA(:, i);
    A(4, 4, i) = 1;
    
    B(1 : 3, 1 : 3, i) = RC(:, :, i);
    B(1 : 3, 4, i) = tC(:, i);
    B(4, 4, i) = 1;
end

[U, Y, f] = unified_fmincon_AX_YB(A, B, iter);
RY_ = Y(1 : 3, 1 : 3);
tY_ = Y(1 : 3, 4);

iter = 128;
[X_, W_, rho_, fs] = opt_fminunc_par(iter, RY_, tY_, RA, tA, RA_dot, tA_dot, RC, tC, RC_dot, tC_dot, omega);

X_

W_

rho_

sfs = sort(fs);
index = sfs(1);
fss = index;
for i = 1 : iter
    if(abs(sfs(i) - index) > 1e-4)
        index = sfs(i);
        fss = [fss, index];
    end
end
fss


figure(1);
plot(time, omega(1, :), '-.', 'LineWidth', 2); hold on
plot(time, omega(2, :), '--', 'LineWidth', 2); hold on
plot(time, omega(3, :), '-', 'LineWidth', 2); hold off
xlabel('Time (s)');
ylabel('Angular Rate (rad/s)');
title('$\bf{\omega}$', 'Interpreter', 'LaTeX');
xlim([0 max(time)]);
legend('X', 'Y', 'Z');

figure(2);
plot(time, tA(1, :), '-.', 'LineWidth', 2); hold on
plot(time, tA(2, :), '--', 'LineWidth', 2); hold on
plot(time, tA(3, :), '-', 'LineWidth', 2); hold off
xlabel('Time (s)');
ylabel('m/s');
title('$\bf{t_A}$', 'Interpreter', 'LaTeX');
xlim([0 max(time)]);
legend('X', 'Y', 'Z');

figure(3);
plot(time, tC(1, :), '-.', 'LineWidth', 2); hold on
plot(time, tC(2, :), '--', 'LineWidth', 2); hold on
plot(time, tC(3, :), '-', 'LineWidth', 2); hold off
xlabel('Time (s)');
ylabel('m/s');
title('$\bf{t_C}$', 'Interpreter', 'LaTeX');
xlim([0 max(time)]);
legend('X', 'Y', 'Z');


[e1, e2, e3] = dcm2angle(RC, 'XYZ');
euler_RC = [e1, e2, e3].' * 180 / pi;
figure(4);
plot(time, euler_RC(1, :), '-.', 'LineWidth', 2); hold on
plot(time, euler_RC(2, :), '--', 'LineWidth', 2); hold on
plot(time, euler_RC(3, :), '-', 'LineWidth', 2); hold off
xlabel('Time (s)');
ylabel('deg');
title('$\bf{R_C}$', 'Interpreter', 'LaTeX');
xlim([0 max(time)]);
legend('Roll', 'Pitch', 'Yaw');


[e1, e2, e3] = dcm2angle(RA, 'XYZ');
euler_RA = [e1, e2, e3].' * 180 / pi;
figure(5);
plot(time, euler_RA(1, :), '-.', 'LineWidth', 2); hold on
plot(time, euler_RA(2, :), '--', 'LineWidth', 2); hold on
plot(time, euler_RA(3, :), '-', 'LineWidth', 2); hold off
xlabel('Time (s)');
ylabel('deg');
title('$\bf{R_A}$', 'Interpreter', 'LaTeX');
xlim([0 max(time)]);
legend('Roll', 'Pitch', 'Yaw');


figure(6);
plot(time, reshape(RC_dot(1, 1, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(1, 2, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(1, 3, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(2, 1, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(2, 2, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(2, 3, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(3, 1, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(3, 2, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RC_dot(3, 3, :), [len, 1]), '-.', 'LineWidth', 2); hold on
xlabel('Time (s)');
title('$\bf{\dot{R}_C}$', 'Interpreter', 'LaTeX');
xlim([0 max(time)]);


figure(7);
plot(time, reshape(RA_dot(1, 1, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(1, 2, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(1, 3, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(2, 1, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(2, 2, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(2, 3, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(3, 1, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(3, 2, :), [len, 1]), '-.', 'LineWidth', 2); hold on
plot(time, reshape(RA_dot(3, 3, :), [len, 1]), '-.', 'LineWidth', 2); hold on
xlabel('Time (s)');
title('$\bf{\dot{R}_A}$', 'Interpreter', 'LaTeX');
xlim([0 max(time)]);


function [X, W, rho, fs] = opt_fminunc_par(iter, RY, tY, RA, tA, RA_dot, tA_dot, RC, tC, RC_dot, tC_dot, omega)
    syms RX11 RX12 RX13 real
    syms RX21 RX22 RX23 real
    syms RX31 RX32 RX33 real
    syms tX0 tX1 tX2 real
    syms RW11 RW12 RW13 real
    syms RW21 RW22 RW23 real
    syms RW31 RW32 RW33 real
    syms tW0 tW1 tW2 real
    syms rho real

    RX = [
        RX11, RX12, RX13;
        RX21, RX22, RX23;
        RX31, RX32, RX33;
        ];
    RW = [
        RW11, RW12, RW13;
        RW21, RW22, RW23;
        RW31, RW32, RW33;
        ];
    tX = [tX0; tX1; tX2];
    tW = [tW0; tW1; tW2];
    len = size(RA, 3);
    J = vpa(expand(J_func_fast(RX, RY, RW, tX, tY, tW, rho, RA, tA, RA_dot, tA_dot, RC, tC, RC_dot, tC_dot, omega)), 32);
    matlabFunction(J, 'File', 'J_func_fminunc_par_.m');
    
    Xs = zeros(4, 4, iter);
    Ws = zeros(4, 4, iter);
    rhos = zeros(iter, 1);
    fs = zeros(iter, 1);
    parfor i = 1 : iter
        x0 = [randn(12, 1); 1];
        opt = optimoptions('fminunc', 'Display', 'off', ...
            'Algorithm', 'quasi-newton', ...
            'SpecifyObjectiveGradient', false, ...
            'OptimalityTolerance', 1e-200, ...
            'UseParallel', false);
            
        [x1, f] = fminunc(@J_func_fminunc_par, x0, opt);

        R = expm(times_(x1(1 : 3), 3));
        t = x1(4 : 6);
        X = [R, t;
             zeros(1, 3), 1];
        R = expm(times_(x1(7 : 9), 3));
        t = x1(10 : 12);
        W = [R, t;
             zeros(1, 3), 1];
        rhos(i) = x1(13); 
        Xs(:, :, i) = X;
        Ws(:, :, i) = W;
        fs(i) = f;
    end
    
    [minimum, idx] = sort(fs);
    X = Xs(:, :, idx(1));
    W = Ws(:, :, idx(1));
    rho = rhos(idx(1), 1);
end


function J = J_func_fminunc_par(x1)
RX = expm(times_(x1(1 : 3), 3));
tX = x1(4 : 6);
RW = expm(times_(x1(7 : 9), 3));
tW = x1(10 : 12); 
rho = x1(13);

J = J_func_fminunc_par_(RW(1, 1), RW(1, 2), RW(1, 3), ...
                        RW(2, 1), RW(2, 2), RW(2, 3), ...
                        RW(3, 1), RW(3, 2), RW(3, 3), ...
                        RX(1, 1), RX(1, 2), RX(1, 3), ...
                        RX(2, 1), RX(2, 2), RX(2, 3), ...
                        RX(3, 1), RX(3, 2), RX(3, 3), ...
                        rho, ...
                        tW(1), tW(2), tW(3), ...
                        tX(1), tX(2), tX(3));
end


function J = J_func(RX, RY, RCB, tX, tY, tCB, rho, RA, tA, RA_dot, tA_dot, RC, tC, RC_dot, tC_dot, omega)
len = size(RA, 3);
J = 0;
for i = 1 : len
    JJ = 0;
    res = RA(:, :, i) * RX - RY * RC(:, :, i) * RCB;
    JJ = JJ + trace(res.' * res);
    res = RA_dot(:, :, i) * RX - RY * RC(:, :, i) * RCB * skew(omega(:, i));
    JJ = JJ + trace(res.' * res);
    res = RA(:, :, i) * tX + tA(:, i) - RY * (RC(:, :, i) * tCB + rho * tC(:, i)) - tY;
    JJ = JJ + res.' * res;
    res = RA_dot(:, :, i) * tX + tA_dot(:, i) - RY * (RC_dot(:, :, i) * tCB + rho * tC_dot(:, i));
    JJ = JJ + res.' * res;
    J = J + 1 / len * JJ;
end
end


function J = J_func_fast(RX, RY, RCB, tX, tY, tCB, rho, RA, tA, RA_dot, tA_dot, RC, tC, RC_dot, tC_dot, omega)
% for i = 1 : 3
%     for j = 1 : 3
%         str = sprintf('RX%d%d = RX(%d, %d);', i, j, i, j);
%         eval(str);
%         
%         str = sprintf('RY%d%d = RY(%d, %d);', i, j, i, j);
%         eval(str);
%         
%         str = sprintf('RW%d%d = RCB(%d, %d);', i, j, i, j);
%         eval(str);
%     end
%     str = sprintf('tX%d = tX(%d);', i - 1, i);
%     eval(str);
%     
%     str = sprintf('tY%d = tY(%d);', i - 1, i);
%     eval(str);
%     
%     str = sprintf('tW%d = tCB(%d);', i - 1, i);
%     eval(str);
% end

% for i = 1 : 3
%     for j = 1 : 3
%         fprintf('RX%d%d = RX(%d, %d);\n', i, j, i, j);
%         fprintf('RY%d%d = RY(%d, %d);\n', i, j, i, j);
%         fprintf('RW%d%d = RCB(%d, %d);\n', i, j, i, j);
%     end
%     fprintf('tX%d = tX(%d);\n', i - 1, i);
%     fprintf('tY%d = tY(%d);\n', i - 1, i);
%     fprintf('tW%d = tCB(%d);\n', i - 1, i);
% end

RX11 = RX(1, 1);
RY11 = RY(1, 1);
RW11 = RCB(1, 1);
RX12 = RX(1, 2);
RY12 = RY(1, 2);
RW12 = RCB(1, 2);
RX13 = RX(1, 3);
RY13 = RY(1, 3);
RW13 = RCB(1, 3);
tX0 = tX(1);
tY0 = tY(1);
tW0 = tCB(1);
RX21 = RX(2, 1);
RY21 = RY(2, 1);
RW21 = RCB(2, 1);
RX22 = RX(2, 2);
RY22 = RY(2, 2);
RW22 = RCB(2, 2);
RX23 = RX(2, 3);
RY23 = RY(2, 3);
RW23 = RCB(2, 3);
tX1 = tX(2);
tY1 = tY(2);
tW1 = tCB(2);
RX31 = RX(3, 1);
RY31 = RY(3, 1);
RW31 = RCB(3, 1);
RX32 = RX(3, 2);
RY32 = RY(3, 2);
RW32 = RCB(3, 2);
RX33 = RX(3, 3);
RY33 = RY(3, 3);
RW33 = RCB(3, 3);
tX2 = tX(3);
tY2 = tY(3);
tW2 = tCB(3);

mon = [ RX11^2, RX11*RX21, RX11*RX31, RW11*RX11, RW21*RX11, RW31*RX11, RW12*RX11, RW22*RX11, RW32*RX11, RW13*RX11, RW23*RX11, RW33*RX11, RX21^2, RX21*RX31, RW11*RX21, RW21*RX21, RW31*RX21, RW12*RX21, RW22*RX21, RW32*RX21, RW13*RX21, RW23*RX21, RW33*RX21, RX31^2, RW11*RX31, RW21*RX31, RW31*RX31, RW12*RX31, RW22*RX31, RW32*RX31, RW13*RX31, RW23*RX31, RW33*RX31, RX12^2, RX12*RX22, RX12*RX32, RW11*RX12, RW21*RX12, RW31*RX12, RW12*RX12, RW22*RX12, RW32*RX12, RW13*RX12, RW23*RX12, RW33*RX12, RX22^2, RX22*RX32, RW11*RX22, RW21*RX22, RW31*RX22, RW12*RX22, RW22*RX22, RW32*RX22, RW13*RX22, RW23*RX22, RW33*RX22, RX32^2, RW11*RX32, RW21*RX32, RW31*RX32, RW12*RX32, RW22*RX32, RW32*RX32, RW13*RX32, RW23*RX32, RW33*RX32, RX13^2, RX13*RX23, RX13*RX33, RW11*RX13, RW21*RX13, RW31*RX13, RW12*RX13, RW22*RX13, RW32*RX13, RW13*RX13, RW23*RX13, RW33*RX13, RX23^2, RX23*RX33, RW11*RX23, RW21*RX23, RW31*RX23, RW12*RX23, RW22*RX23, RW32*RX23, RW13*RX23, RW23*RX23, RW33*RX23, RX33^2, RW11*RX33, RW21*RX33, RW31*RX33, RW12*RX33, RW22*RX33, RW32*RX33, RW13*RX33, RW23*RX33, RW33*RX33, RW11^2, RW11*RW21, RW11*RW31, RW11*RW12, RW11*RW22, RW11*RW32, RW11*RW13, RW11*RW23, RW11*RW33, RW21^2, RW21*RW31, RW12*RW21, RW21*RW22, RW21*RW32, RW13*RW21, RW21*RW23, RW21*RW33, RW31^2, RW12*RW31, RW22*RW31, RW31*RW32, RW13*RW31, RW23*RW31, RW31*RW33, RW12^2, RW12*RW22, RW12*RW32, RW12*RW13, RW12*RW23, RW12*RW33, RW22^2, RW22*RW32, RW13*RW22, RW22*RW23, RW22*RW33, RW32^2, RW13*RW32, RW23*RW32, RW32*RW33, RW13^2, RW13*RW23, RW13*RW33, RW23^2, RW23*RW33, RW33^2, tX0^2, tX0*tX1, tX0*tX2, tW0*tX0, tW1*tX0, tW2*tX0, rho*tX0, tX0, tX1^2, tX1*tX2, tW0*tX1, tW1*tX1, tW2*tX1, rho*tX1, tX1, tX2^2, tW0*tX2, tW1*tX2, tW2*tX2, rho*tX2, tX2, tW0^2, tW0*tW1, tW0*tW2, rho*tW0, tW0, tW1^2, tW1*tW2, rho*tW1, tW1, tW2^2, rho*tW2, tW2, rho^2, rho, 1];
len = size(RA, 3);
all = zeros(1, 180);

for n = 1 : len
    w1 = omega(1, n);
    w2 = omega(2, n);
    w3 = omega(3, n);
    RA11 = RA(1, 1, n); 
    RAd11 = RA_dot(1, 1, n);
    RC11 = RC(1, 1, n);
    RCd11 = RC_dot(1, 1, n);
    RA12 = RA(1, 2, n);
    RAd12 = RA_dot(1, 2, n);
    RC12 = RC(1, 2, n);
    RCd12 = RC_dot(1, 2, n);
    RA13 = RA(1, 3, n);
    RAd13 = RA_dot(1, 3, n);
    RC13 = RC(1, 3, n);
    RCd13 = RC_dot(1, 3, n);
    tA0 = tA(1, n);
    tAd0 = tA_dot(1, n);
    tC0 = tC(1, n);
    tCd0 = tC_dot(1, n);
    RA21 = RA(2, 1, n);
    RAd21 = RA_dot(2, 1, n);
    RC21 = RC(2, 1, n);
    RCd21 = RC_dot(2, 1, n);
    RA22 = RA(2, 2, n);
    RAd22 = RA_dot(2, 2, n);
    RC22 = RC(2, 2, n);
    RCd22 = RC_dot(2, 2, n);
    RA23 = RA(2, 3, n);
    RAd23 = RA_dot(2, 3, n);
    RC23 = RC(2, 3, n);
    RCd23 = RC_dot(2, 3, n);
    tA1 = tA(2, n);
    tAd1 = tA_dot(2, n);
    tC1 = tC(2, n);
    tCd1 = tC_dot(2, n);
    RA31 = RA(3, 1, n);
    RAd31 = RA_dot(3, 1, n);
    RC31 = RC(3, 1, n);
    RCd31 = RC_dot(3, 1, n);
    RA32 = RA(3, 2, n);
    RAd32 = RA_dot(3, 2, n);
    RC32 = RC(3, 2, n);
    RCd32 = RC_dot(3, 2, n);
    RA33 = RA(3, 3, n);
    RAd33 = RA_dot(3, 3, n);
    RC33 = RC(3, 3, n);
    RCd33 = RC_dot(3, 3, n);
    tA2 = tA(3, n);
    tAd2 = tA_dot(3, n);
    tC2 = tC(3, n);
    tCd2 = tC_dot(3, n);
    
%     for i = 1 : 3
%         for j = 1 : 3
%             fprintf('RA%d%d = RA(%d, %d, %d); \n', i, j, i, j, n);
%             fprintf('RAd%d%d = RA_dot(%d, %d, %d); \n', i, j, i, j, n);
%             fprintf('RC%d%d = RC(%d, %d, %d); \n', i, j, i, j, n);
%             fprintf('RCd%d%d = RC_dot(%d, %d, %d); \n', i, j, i, j, n);
%         end
%         fprintf('tA%d = tA(%d, %d);\n', i - 1, i, n);
%         fprintf('tAd%d = tA_dot(%d, %d);\n', i - 1, i, n);
%         fprintf('tC%d = tC(%d, %d);\n', i - 1, i, n);
%         fprintf('tCd%d = tC_dot(%d, %d);\n', i - 1, i, n);
%     end
%     for i = 1 : 3
%         for j = 1 : 3
%         
%             str = sprintf('RA%d%d = RA(%d, %d, %d);', i, j, i, j, n);
%             eval(str);
%         
%             str = sprintf('RAd%d%d = RA_dot(%d, %d, %d);', i, j, i, j, n);
%             eval(str);
%         
%             str = sprintf('RC%d%d = RC(%d, %d, %d);', i, j, i, j, n);
%             eval(str);
%         
%             str = sprintf('RCd%d%d = RC_dot(%d, %d, %d);', i, j, i, j, n);
%             eval(str);
%         end
%         str = sprintf('tA%d = tA(%d, %d);', i - 1, i, n);
%         eval(str);
%         
%         str = sprintf('tAd%d = tA_dot(%d, %d);', i - 1, i, n);
%         eval(str);
%         
%         str = sprintf('tC%d = tC(%d, %d);', i - 1, i, n);
%         eval(str);
%         
%         str = sprintf('tCd%d = tC_dot(%d, %d);', i - 1, i, n);
%         eval(str);
%     end
    coeff = [ RAd11^2 + RAd21^2 + RAd31^2 + 1, 2*RAd11*RAd12 + 2*RAd21*RAd22 + 2*RAd31*RAd32, 2*RAd11*RAd13 + 2*RAd21*RAd23 + 2*RAd31*RAd33, - 2*RA11*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA21*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA31*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA11*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA21*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA31*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA11*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA21*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA31*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RAd11*w3*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd21*w3*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd31*w3*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd11*w3*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd21*w3*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd31*w3*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd11*w3*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd21*w3*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd31*w3*(RC13*RY31 + RC23*RY32 + RC33*RY33), 2*RAd11*w2*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd21*w2*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd31*w2*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd11*w2*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd21*w2*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd31*w2*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd11*w2*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd21*w2*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd31*w2*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd12^2 + RAd22^2 + RAd32^2 + 1, 2*RAd12*RAd13 + 2*RAd22*RAd23 + 2*RAd32*RAd33, - 2*RA12*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA22*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA32*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA12*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA22*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA32*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA12*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA22*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA32*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RAd12*w3*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd22*w3*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd32*w3*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd12*w3*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd22*w3*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd32*w3*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd12*w3*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd22*w3*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd32*w3*(RC13*RY31 + RC23*RY32 + RC33*RY33), 2*RAd12*w2*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd22*w2*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd32*w2*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd12*w2*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd22*w2*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd32*w2*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd12*w2*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd22*w2*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd32*w2*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd13^2 + RAd23^2 + RAd33^2 + 1, - 2*RA13*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA23*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA33*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA13*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA23*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA33*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA13*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA23*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA33*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RAd13*w3*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd23*w3*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd33*w3*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd13*w3*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd23*w3*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd33*w3*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd13*w3*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd23*w3*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd33*w3*(RC13*RY31 + RC23*RY32 + RC33*RY33), 2*RAd13*w2*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd23*w2*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd33*w2*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd13*w2*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd23*w2*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd33*w2*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd13*w2*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd23*w2*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd33*w2*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd11^2 + RAd21^2 + RAd31^2 + 1, 2*RAd11*RAd12 + 2*RAd21*RAd22 + 2*RAd31*RAd32, 2*RAd11*RAd13 + 2*RAd21*RAd23 + 2*RAd31*RAd33, 2*RAd11*w3*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd21*w3*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd31*w3*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd11*w3*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd21*w3*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd31*w3*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd11*w3*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd21*w3*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd31*w3*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RA11*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA21*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA31*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA11*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA21*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA31*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA11*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA21*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA31*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RAd11*w1*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd21*w1*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd31*w1*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd11*w1*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd21*w1*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd31*w1*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd11*w1*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd21*w1*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd31*w1*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd12^2 + RAd22^2 + RAd32^2 + 1, 2*RAd12*RAd13 + 2*RAd22*RAd23 + 2*RAd32*RAd33, 2*RAd12*w3*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd22*w3*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd32*w3*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd12*w3*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd22*w3*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd32*w3*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd12*w3*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd22*w3*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd32*w3*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RA12*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA22*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA32*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA12*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA22*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA32*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA12*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA22*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA32*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RAd12*w1*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd22*w1*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd32*w1*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd12*w1*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd22*w1*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd32*w1*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd12*w1*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd22*w1*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd32*w1*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd13^2 + RAd23^2 + RAd33^2 + 1, 2*RAd13*w3*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd23*w3*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd33*w3*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd13*w3*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd23*w3*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd33*w3*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd13*w3*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd23*w3*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd33*w3*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RA13*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA23*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA33*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA13*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA23*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA33*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA13*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA23*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA33*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RAd13*w1*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd23*w1*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd33*w1*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd13*w1*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd23*w1*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd33*w1*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd13*w1*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd23*w1*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd33*w1*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd11^2 + RAd21^2 + RAd31^2 + 1, 2*RAd11*RAd12 + 2*RAd21*RAd22 + 2*RAd31*RAd32, 2*RAd11*RAd13 + 2*RAd21*RAd23 + 2*RAd31*RAd33, - 2*RAd11*w2*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd21*w2*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd31*w2*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd11*w2*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd21*w2*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd31*w2*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd11*w2*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd21*w2*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd31*w2*(RC13*RY31 + RC23*RY32 + RC33*RY33), 2*RAd11*w1*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd21*w1*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd31*w1*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd11*w1*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd21*w1*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd31*w1*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd11*w1*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd21*w1*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd31*w1*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RA11*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA21*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA31*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA11*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA21*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA31*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA11*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA21*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA31*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd12^2 + RAd22^2 + RAd32^2 + 1, 2*RAd12*RAd13 + 2*RAd22*RAd23 + 2*RAd32*RAd33, - 2*RAd12*w2*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd22*w2*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd32*w2*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd12*w2*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd22*w2*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd32*w2*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd12*w2*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd22*w2*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd32*w2*(RC13*RY31 + RC23*RY32 + RC33*RY33), 2*RAd12*w1*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd22*w1*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd32*w1*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd12*w1*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd22*w1*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd32*w1*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd12*w1*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd22*w1*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd32*w1*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RA12*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA22*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA32*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA12*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA22*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA32*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA12*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA22*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA32*(RC13*RY31 + RC23*RY32 + RC33*RY33), RAd13^2 + RAd23^2 + RAd33^2 + 1, - 2*RAd13*w2*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RAd23*w2*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RAd33*w2*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RAd13*w2*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RAd23*w2*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RAd33*w2*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RAd13*w2*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RAd23*w2*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RAd33*w2*(RC13*RY31 + RC23*RY32 + RC33*RY33), 2*RAd13*w1*(RC11*RY11 + RC21*RY12 + RC31*RY13) + 2*RAd23*w1*(RC11*RY21 + RC21*RY22 + RC31*RY23) + 2*RAd33*w1*(RC11*RY31 + RC21*RY32 + RC31*RY33), 2*RAd13*w1*(RC12*RY11 + RC22*RY12 + RC32*RY13) + 2*RAd23*w1*(RC12*RY21 + RC22*RY22 + RC32*RY23) + 2*RAd33*w1*(RC12*RY31 + RC22*RY32 + RC32*RY33), 2*RAd13*w1*(RC13*RY11 + RC23*RY12 + RC33*RY13) + 2*RAd23*w1*(RC13*RY21 + RC23*RY22 + RC33*RY23) + 2*RAd33*w1*(RC13*RY31 + RC23*RY32 + RC33*RY33), - 2*RA13*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA23*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA33*(RC11*RY31 + RC21*RY32 + RC31*RY33), - 2*RA13*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA23*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA33*(RC12*RY31 + RC22*RY32 + RC32*RY33), - 2*RA13*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA23*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA33*(RC13*RY31 + RC23*RY32 + RC33*RY33), w2^2 + w3^2 + 1, 0, 0, -2*w1*w2, 0, 0, -2*w1*w3, 0, 0, w2^2 + w3^2 + 1, 0, 0, -2*w1*w2, 0, 0, -2*w1*w3, 0, w2^2 + w3^2 + 1, 0, 0, -2*w1*w2, 0, 0, -2*w1*w3, w1^2 + w3^2 + 1, 0, 0, -2*w2*w3, 0, 0, w1^2 + w3^2 + 1, 0, 0, -2*w2*w3, 0, w1^2 + w3^2 + 1, 0, 0, -2*w2*w3, w1^2 + w2^2 + 1, 0, 0, w1^2 + w2^2 + 1, 0, w1^2 + w2^2 + 1, RAd11^2 + RAd21^2 + RAd31^2 + 1, 2*RAd11*RAd12 + 2*RAd21*RAd22 + 2*RAd31*RAd32, 2*RAd11*RAd13 + 2*RAd21*RAd23 + 2*RAd31*RAd33, - 2*RA11*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA21*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA31*(RC11*RY31 + RC21*RY32 + RC31*RY33) - 2*RAd11*(RCd11*RY11 + RCd21*RY12 + RCd31*RY13) - 2*RAd21*(RCd11*RY21 + RCd21*RY22 + RCd31*RY23) - 2*RAd31*(RCd11*RY31 + RCd21*RY32 + RCd31*RY33), - 2*RA11*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA21*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA31*(RC12*RY31 + RC22*RY32 + RC32*RY33) - 2*RAd11*(RCd12*RY11 + RCd22*RY12 + RCd32*RY13) - 2*RAd21*(RCd12*RY21 + RCd22*RY22 + RCd32*RY23) - 2*RAd31*(RCd12*RY31 + RCd22*RY32 + RCd32*RY33), - 2*RA11*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA21*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA31*(RC13*RY31 + RC23*RY32 + RC33*RY33) - 2*RAd11*(RCd13*RY11 + RCd23*RY12 + RCd33*RY13) - 2*RAd21*(RCd13*RY21 + RCd23*RY22 + RCd33*RY23) - 2*RAd31*(RCd13*RY31 + RCd23*RY32 + RCd33*RY33), - 2*RA11*(RY11*tC0 + RY12*tC1 + RY13*tC2) - 2*RA21*(RY21*tC0 + RY22*tC1 + RY23*tC2) - 2*RA31*(RY31*tC0 + RY32*tC1 + RY33*tC2) - 2*RAd11*(RY11*tCd0 + RY12*tCd1 + RY13*tCd2) - 2*RAd21*(RY21*tCd0 + RY22*tCd1 + RY23*tCd2) - 2*RAd31*(RY31*tCd0 + RY32*tCd1 + RY33*tCd2), 2*RAd11*tAd0 + 2*RAd21*tAd1 + 2*RAd31*tAd2 + 2*RA11*(tA0 - tY0) + 2*RA21*(tA1 - tY1) + 2*RA31*(tA2 - tY2), RAd12^2 + RAd22^2 + RAd32^2 + 1, 2*RAd12*RAd13 + 2*RAd22*RAd23 + 2*RAd32*RAd33, - 2*RA12*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA22*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA32*(RC11*RY31 + RC21*RY32 + RC31*RY33) - 2*RAd12*(RCd11*RY11 + RCd21*RY12 + RCd31*RY13) - 2*RAd22*(RCd11*RY21 + RCd21*RY22 + RCd31*RY23) - 2*RAd32*(RCd11*RY31 + RCd21*RY32 + RCd31*RY33), - 2*RA12*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA22*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA32*(RC12*RY31 + RC22*RY32 + RC32*RY33) - 2*RAd12*(RCd12*RY11 + RCd22*RY12 + RCd32*RY13) - 2*RAd22*(RCd12*RY21 + RCd22*RY22 + RCd32*RY23) - 2*RAd32*(RCd12*RY31 + RCd22*RY32 + RCd32*RY33), - 2*RA12*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA22*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA32*(RC13*RY31 + RC23*RY32 + RC33*RY33) - 2*RAd12*(RCd13*RY11 + RCd23*RY12 + RCd33*RY13) - 2*RAd22*(RCd13*RY21 + RCd23*RY22 + RCd33*RY23) - 2*RAd32*(RCd13*RY31 + RCd23*RY32 + RCd33*RY33), - 2*RA12*(RY11*tC0 + RY12*tC1 + RY13*tC2) - 2*RA22*(RY21*tC0 + RY22*tC1 + RY23*tC2) - 2*RA32*(RY31*tC0 + RY32*tC1 + RY33*tC2) - 2*RAd12*(RY11*tCd0 + RY12*tCd1 + RY13*tCd2) - 2*RAd22*(RY21*tCd0 + RY22*tCd1 + RY23*tCd2) - 2*RAd32*(RY31*tCd0 + RY32*tCd1 + RY33*tCd2), 2*RAd12*tAd0 + 2*RAd22*tAd1 + 2*RAd32*tAd2 + 2*RA12*(tA0 - tY0) + 2*RA22*(tA1 - tY1) + 2*RA32*(tA2 - tY2), RAd13^2 + RAd23^2 + RAd33^2 + 1, - 2*RA13*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*RA23*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*RA33*(RC11*RY31 + RC21*RY32 + RC31*RY33) - 2*RAd13*(RCd11*RY11 + RCd21*RY12 + RCd31*RY13) - 2*RAd23*(RCd11*RY21 + RCd21*RY22 + RCd31*RY23) - 2*RAd33*(RCd11*RY31 + RCd21*RY32 + RCd31*RY33), - 2*RA13*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*RA23*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*RA33*(RC12*RY31 + RC22*RY32 + RC32*RY33) - 2*RAd13*(RCd12*RY11 + RCd22*RY12 + RCd32*RY13) - 2*RAd23*(RCd12*RY21 + RCd22*RY22 + RCd32*RY23) - 2*RAd33*(RCd12*RY31 + RCd22*RY32 + RCd32*RY33), - 2*RA13*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*RA23*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*RA33*(RC13*RY31 + RC23*RY32 + RC33*RY33) - 2*RAd13*(RCd13*RY11 + RCd23*RY12 + RCd33*RY13) - 2*RAd23*(RCd13*RY21 + RCd23*RY22 + RCd33*RY23) - 2*RAd33*(RCd13*RY31 + RCd23*RY32 + RCd33*RY33), - 2*RA13*(RY11*tC0 + RY12*tC1 + RY13*tC2) - 2*RA23*(RY21*tC0 + RY22*tC1 + RY23*tC2) - 2*RA33*(RY31*tC0 + RY32*tC1 + RY33*tC2) - 2*RAd13*(RY11*tCd0 + RY12*tCd1 + RY13*tCd2) - 2*RAd23*(RY21*tCd0 + RY22*tCd1 + RY23*tCd2) - 2*RAd33*(RY31*tCd0 + RY32*tCd1 + RY33*tCd2), 2*RAd13*tAd0 + 2*RAd23*tAd1 + 2*RAd33*tAd2 + 2*RA13*(tA0 - tY0) + 2*RA23*(tA1 - tY1) + 2*RA33*(tA2 - tY2), RCd11^2 + RCd21^2 + RCd31^2 + 1, 2*RCd11*RCd12 + 2*RCd21*RCd22 + 2*RCd31*RCd32, 2*RCd11*RCd13 + 2*RCd21*RCd23 + 2*RCd31*RCd33, 2*RC11*tC0 + 2*RC21*tC1 + 2*RC31*tC2 + 2*RCd11*tCd0 + 2*RCd21*tCd1 + 2*RCd31*tCd2, - 2*tAd0*(RCd11*RY11 + RCd21*RY12 + RCd31*RY13) - 2*tAd1*(RCd11*RY21 + RCd21*RY22 + RCd31*RY23) - 2*tAd2*(RCd11*RY31 + RCd21*RY32 + RCd31*RY33) - 2*(tA0 - tY0)*(RC11*RY11 + RC21*RY12 + RC31*RY13) - 2*(tA1 - tY1)*(RC11*RY21 + RC21*RY22 + RC31*RY23) - 2*(tA2 - tY2)*(RC11*RY31 + RC21*RY32 + RC31*RY33), RCd12^2 + RCd22^2 + RCd32^2 + 1, 2*RCd12*RCd13 + 2*RCd22*RCd23 + 2*RCd32*RCd33, 2*RC12*tC0 + 2*RC22*tC1 + 2*RC32*tC2 + 2*RCd12*tCd0 + 2*RCd22*tCd1 + 2*RCd32*tCd2, - 2*tAd0*(RCd12*RY11 + RCd22*RY12 + RCd32*RY13) - 2*tAd1*(RCd12*RY21 + RCd22*RY22 + RCd32*RY23) - 2*tAd2*(RCd12*RY31 + RCd22*RY32 + RCd32*RY33) - 2*(tA0 - tY0)*(RC12*RY11 + RC22*RY12 + RC32*RY13) - 2*(tA1 - tY1)*(RC12*RY21 + RC22*RY22 + RC32*RY23) - 2*(tA2 - tY2)*(RC12*RY31 + RC22*RY32 + RC32*RY33), RCd13^2 + RCd23^2 + RCd33^2 + 1, 2*RC13*tC0 + 2*RC23*tC1 + 2*RC33*tC2 + 2*RCd13*tCd0 + 2*RCd23*tCd1 + 2*RCd33*tCd2, - 2*tAd0*(RCd13*RY11 + RCd23*RY12 + RCd33*RY13) - 2*tAd1*(RCd13*RY21 + RCd23*RY22 + RCd33*RY23) - 2*tAd2*(RCd13*RY31 + RCd23*RY32 + RCd33*RY33) - 2*(tA0 - tY0)*(RC13*RY11 + RC23*RY12 + RC33*RY13) - 2*(tA1 - tY1)*(RC13*RY21 + RC23*RY22 + RC33*RY23) - 2*(tA2 - tY2)*(RC13*RY31 + RC23*RY32 + RC33*RY33), tC0^2 + tC1^2 + tC2^2 + tCd0^2 + tCd1^2 + tCd2^2, - 2*tAd0*(RY11*tCd0 + RY12*tCd1 + RY13*tCd2) - 2*tAd1*(RY21*tCd0 + RY22*tCd1 + RY23*tCd2) - 2*tAd2*(RY31*tCd0 + RY32*tCd1 + RY33*tCd2) - 2*(tA0 - tY0)*(RY11*tC0 + RY12*tC1 + RY13*tC2) - 2*(tA1 - tY1)*(RY21*tC0 + RY22*tC1 + RY23*tC2) - 2*(tA2 - tY2)*(RY31*tC0 + RY32*tC1 + RY33*tC2), (tA0 - tY0)^2 + (tA1 - tY1)^2 + (tA2 - tY2)^2 + tAd0^2 + tAd1^2 + tAd2^2];
    all = all + 1 / len * coeff;
end
J = all * mon.';
end