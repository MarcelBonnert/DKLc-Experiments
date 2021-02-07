% clear all

%% Load necessary data
load('normalizationData.mat')
load('pumpCharacteristic.mat')

%% Pump charactersistic
voltage2vdot_P1 = @(x) m1*x+b1;
voltage2vdot_P2 = @(x) m2*x+b2;

vdot2voltage_P1 = @(x) x/m1-b1/m1;
vdot2voltage_P2 = @(x) x/m2-b2/m2;

%%
Ts = 5;

order_x = 3;
order_u = 2;

path_full = '../../experiments/threeTank/models/';

content = dir(path_full);

runs = {...
        'S1_n_6';... 
        'S2_n_6';... 
        'S3_n_4';... 
        'S4_n_8'}; 
    
Q = [ones(1, order_x); ones(1, order_x); ones(1, order_x); ones(1, order_x)];
R = eye(order_u);

Q(:,2) = Q(:,2)*10;
Q(:,1) = Q(:,1)*1;
Q(:,3) = Q(:,3)*1;

% first row - mean; second row: variance
normData_x = [  mean_data_ds(3:end);
                var_data_ds(3:end)/1000]/1000;
normData_u = [  mean_data_ds(1:2);
                var_data_ds(1:2)];

xap = mean_data_ds(3:end);
uap = mean_data_ds(1:2);

t = 0:5:1200;
x0 = [120 60 130]/1000';

f1 = figure;
hold on
t5 = zeros(1,2);

for i = 1:length(runs)
    Kz = readmatrix([path_full runs{i} ...
            '/K.csv'], 'Delimiter', ' ')';
    Kw = readmatrix([path_full runs{i} ...
            '/L.csv'], 'Delimiter', ' ')';
    
    stateEncoderWeightFiles = dir([path_full runs{i} ...
       '/state_encoderWeights*' ]);
    stateEncoderWeight_cell = {};
    stateEncoderBiasFiles = dir([path_full runs{i} ...
       '/state_encoderBiases*' ]);
    stateEncoderBias_cell = {};
    for j = 1:length(stateEncoderWeightFiles)
       stateEncoderWeight_cell{end+1} = readmatrix([path_full runs{i} ...
           '/' stateEncoderWeightFiles(j).name ], 'Delimiter', ' ');
       stateEncoderBias_cell{end+1} = readmatrix([path_full runs{i} ...
           '/' stateEncoderBiasFiles(j).name ], 'Delimiter', ' ')';
    end
    
    z0 = g_z(zeros(1, order_x), stateEncoderWeight_cell, stateEncoderBias_cell)';
    
    if contains(runs{i}, 'S3') || contains(runs{i}, 'S4')
        hasLinearState = false;
        stateDecoderWeightFiles = dir([path_full runs{i} ...
           '/state_decoderWeights*' ]);
        stateDecoderWeight_cell = {};
        stateDecoderBiasFiles = dir([path_full runs{i} ...
           '/state_decoderBiases*' ]);
        stateDecoderBias_cell = {};
        for j = 1:length(stateDecoderWeightFiles)
           stateDecoderWeight_cell{end+1} = readmatrix(...
               [path_full runs{i} '/' ...
               stateDecoderWeightFiles(j).name ], 'Delimiter', ' ');
           stateDecoderBias_cell{end+1} = readmatrix(...
               [path_full runs{i} '/' ...
               stateDecoderBiasFiles(j).name ], 'Delimiter', ' ')';
        end
        xapn = (xap/1000 - normData_x(1,:))./(3*sqrt(normData_x(2,:)));

        [Qlqr, xi] = ...
            getQz(xapn, ...
            stateEncoderWeight_cell, stateEncoderBias_cell, ...
            diag(Q(i,:)), Kz, Kw, Ts);

    else
       hasLinearState = true;
       
       Qlqr = [ diag(Q(i,:)) zeros(order_x, size(Kz,2)-order_x);
                zeros(size(Kz,1)-order_x, size(Kz,2))];
        z0 = [zeros(order_x,1); z0];
    end

    if contains(runs{i}, 'S2') || contains(runs{i}, 'S4')
        hasNonlinearInput = true;

        inputEncoderWeightFiles = dir([path_full runs{i}, ...
           '/input_encoderWeights*' ]);
        inputEncoderWeight_cell = {};
        inputEncoderBiasFiles = dir([path_full runs{i} ...
           '/input_encoderBiases*' ]);
        inputEncoderBias_cell = {};
        for j = 1:length(inputEncoderWeightFiles)
           inputEncoderWeight_cell{end+1} = readmatrix(...
               [path_full runs{i} '/' ...
               inputEncoderWeightFiles(j).name ], 'Delimiter', ' ');
           inputEncoderBias_cell{end+1} = readmatrix(...
               [path_full runs{i} '/' ...
               inputEncoderBiasFiles(j).name ], 'Delimiter', ' ')';
        end

        inputDecoderWeightFiles = dir([path_full runs{i} ...
           '/input_decoderWeights*' ]);
        inputDecoderWeight_cell = {};
        inputDecoderBiasFiles = dir([path_full runs{i} ...
           '/input_decoderBiases*' ]);
        inputDecoderBias_cell = {};
        for j = 1:length(inputDecoderWeightFiles)
           inputDecoderWeight_cell{end+1} = readmatrix(...
               [path_full runs{i} '/' ...
               inputDecoderWeightFiles(j).name ], 'Delimiter', ' ');
           inputDecoderBias_cell{end+1} = readmatrix(...
               [path_full runs{i} '/' ...
               inputDecoderBiasFiles(j).name ], 'Delimiter', ' ')';
        end
        
        uapn = (uap - normData_u(1,:))./(3*sqrt(normData_u(2,:)));
        [Rlqr, xi] = ...
            getQz(uapn, ...
            inputEncoderWeight_cell, inputEncoderBias_cell, ...
            R, zeros(2), Kw, Ts);
        
        w0 = g_z(zeros(1, order_u), inputEncoderWeight_cell, inputEncoderBias_cell)';
        
    else
        hasNonlinearInput = false;
        Rlqr = R;
        w0 = zeros(order_u,1);
        
    end
    
    C = dlqr(Kz, Kw, Qlqr, Rlqr);
    
    x = zeros(order_x, length(t));
    x(:,1) = x0;
    uvec = zeros(order_u, length(t)-1,1);
    for k = 2:length(t)
        xn = (x(:,k-1) - normData_x(1,:)')./(3*sqrt(normData_x(2,:))');
        z = g_z(xn', stateEncoderWeight_cell, stateEncoderBias_cell)';
        if hasLinearState
            z = [xn;z];
        end
        
        uc = -C*(z - z0) + w0;
        
        if hasNonlinearInput
            un = g_z(uc', inputDecoderWeight_cell, inputDecoderBias_cell)';
        else
            un = uc;
        end
        
        u = un * 3 .* sqrt(normData_u(2,:))' + normData_u(1,:)';
        
        uvec(:,k-1) = u;
        
        fun = @(t, x) threeTankModel(t, x, [voltage2vdot_P1(u) voltage2vdot_P2(u)]);
        
        sol = rkSolver(fun, t(k-1:k), x(:, k-1));
        x(:, k) = sol(:,2);
    end
    
    
    
    figure(f1);
    for j = 1:order_x
        ax(j) = subplot(order_x+order_u+1,1,j);
        plot(t,x(j,:))
        hold on
        ylabel(['x_' num2str(j)])
    end
    
    for j = 1:order_u
        ax(order_x+j) = subplot(order_x+order_u+1,1,order_x+j);
        plot(t(1:end-1),uvec(j,:))
        hold on
        ylabel(['u_' num2str(j)])
    end
    
    ax(order_x+order_u+1) = subplot(order_x+order_u+1,1,order_x+order_u+1);
    Eu = controlEnergy(uvec' - normData_u(1,:), t);
    disp(['Eu = ' num2str(Eu(end-1))])
    plot(t(1:end-1),Eu)
    hold on
    ylabel('\int_0^t \Delta \mat{u}\left(\tau\right)^2 d \tau')
    xlabel('Time (s)')
    
    t5(i,1) = get_t_barrier(x(1,:), 0.188, t, 0.01);
    t5(i,2) = get_t_barrier(x(2,:), 0.09, t, 0.01);
    t5(i,3) = get_t_barrier(x(3,:), 0.188, t, 0.01);

    clear stateEncoderWeight_cell stateEncoderBias_cell ...
        stateDecoderWeight_cell stateDecoderBias_cell ...
        inputEncoderWeight_cell inputEncoderBias_cell ...
        inputDecoderWeight_cell inputDecoderBias_cell
    
end

figure(f1);
linkaxes(ax, 'x')

%% Linearization
xe = normData_x(1,:);
ue_voltage = normData_u(1,:);
ue = [voltage2vdot_P1(ue_voltage(1)), voltage2vdot_P2(ue_voltage(2))];

syms x1 x2 x3 u1 u2 real;
f = threeTankModel(0,[x1 x2 x3],[u1 u2]);
dfdxu = jacobian(f, [x1 x2 x3, u1 u2]);
A_lin = double(subs(dfdxu(:,1:3), [x1 x2 x3, u1 u2], [xe ue]));
B_lin = double(subs(dfdxu(:,4:end), [x1 x2 x3, u1 u2], [xe ue]));

%% DMD and linearization - tbd - [voltage2vdot_P1(u) voltage2vdot_P2(u)]
Q_lin = diag(Q(1,:));
R_lin = R;

vdot_ap = [voltage2vdot_P1(normData_u(1,1)), voltage2vdot_P1(normData_u(1,2))]';
saturation_vdot = [0, voltage2vdot_P1(0);0, voltage2vdot_P2(0)];
saturation_voltage = [-1, 0;-1 0];

Clin = lqr(A_lin, B_lin, Q_lin, R_lin);
    
x = zeros(order_x, length(t));
x(:,1) = x0;
uvec = zeros(order_u, length(t));
for k = 2:length(t)
    u = -Clin*(x(:,k-1) - normData_x(1,:)') + vdot_ap;
    if ~isempty(saturation_vdot)
        for i = 1:length(u)
            u(i) = sat(u(i), saturation_vdot(i, 1), saturation_vdot(i, 2));
        end
    end
    uvec(:,k-1) = u;
    fun = @(t, x) threeTankModel(t, x, u);

    sol = rkSolver(fun, t(k-1:k), x(:, k-1));
    x(:, k) = sol(:,2);
end

t5(end+1,1) = get_t_barrier(x(1,:), 0.188, t, 0.01);
t5(end,2) = get_t_barrier(x(2,:), 0.09, t, 0.01);
t5(end,3) = get_t_barrier(x(3,:), 0.188, t, 0.01);

% DMD
load('sys_dmd.mat')
Cdmd = dlqr(sys_dmd.A, sys_dmd.B, Q_lin, R_lin);

x_dmd = zeros(order_x, length(t));
x_dmd(:,1) = x0;
uvec_dmd = zeros(order_u, length(t));
for k = 2:length(t)
    xn = (x_dmd(:,k-1)' - normData_x(1,:))./(3*sqrt(normData_x(2,:)));
    xapn = (xap/1000 - normData_x(1,:))./(3*sqrt(normData_x(2,:)));
    
    uc = -Cdmd*xn';

    u = uc * 3 .* sqrt(normData_u(2,:))' + normData_u(1,:)';
    if ~isempty(saturation_voltage)
        for i = 1:length(u)
            u(i) = sat(u(i), saturation_voltage(i, 1), saturation_voltage(i, 2));
        end
    end
    
    uvec_dmd(:,k-1) = [voltage2vdot_P1(u(1)) voltage2vdot_P2(u(2))];

    fun = @(t, x) threeTankModel(t, x, [voltage2vdot_P1(u(1)) voltage2vdot_P2(u(2))]);

    sol = rkSolver(fun, t(k-1:k), x(:, k-1));
    x_dmd(:, k) = sol(:,2);
end

% Plotting
figure(f1);
for i = 1:order_x
    subplot(order_x+order_u+1,1,i)
    plot(t,x(i,:)')
    plot(t,x_dmd(i,:)')
    plot([t(1) t(end)], normData_x(1,i)*ones(1,2))
    hold off
    if i == 1
        title(['Q = diag(' num2str(Q(1,:)) ') - not optimized'])
    end
end
legend('S1','S2','S3','S4','Lin','DMD')

uvec = [vdot2voltage_P1(uvec(1,:)); vdot2voltage_P2(uvec(2,:))];
uvec_dmd = [vdot2voltage_P1(uvec_dmd(1,:)); vdot2voltage_P2(uvec_dmd(2,:))];

for j = 1:order_u
    subplot(order_x+order_u+1,1,order_x+j);
    plot(t,uvec(j,:))
    plot(t,uvec_dmd(j,:))
    hold on
    ylabel(['u_' num2str(j)])
end

subplot(order_x+order_u+1,1,order_x+order_u+1)
Eu = controlEnergy([vdot2voltage_P1(uvec(1,:)); vdot2voltage_P2(uvec(2,:))]'-normData_u(1,:), t);
plot(t(1:end),Eu)
Eu_dmd = controlEnergy([vdot2voltage_P1(uvec_dmd(1,:)); vdot2voltage_P2(uvec_dmd(2,:))]'-normData_u(1,:), t);
plot(t(1:end),Eu_dmd)


t5(end+1,1) = get_t_barrier(x_dmd(1,:), 0.188, t, 0.01);
t5(end,2) = get_t_barrier(x_dmd(2,:), 0.09, t, 0.01);
t5(end,3) = get_t_barrier(x_dmd(3,:), 0.188, t, 0.01)


