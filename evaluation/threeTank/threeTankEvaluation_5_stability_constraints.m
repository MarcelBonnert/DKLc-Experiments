%% Load necessary data
load('normalizationData.mat')
load('pumpCharacteristic.mat')

%% Pump charactersistic
voltage2vdot_P1 = @(x) m1*x+b1;
voltage2vdot_P2 = @(x) m2*x+b2;

vdot2voltage_P1 = @(x) x/m1-b1/m1;
vdot2voltage_P2 = @(x) x/m2-b2/m2;

%%
order_x = 3;
order_u = 2;

path_full = '../../experiments/threeTank/models/';

runs = {...
        'S3_n_7'};
    
data_path = '../../data/threeTank/threeTankTrainingStates.csv';

Q = diag([1 10 1]);
R = eye(2);

normData_x = [  mean_data_ds(3:end);
                var_data_ds(3:end)/1000]/1000;
normData_u = [  mean_data_ds(1:2);
                var_data_ds(1:2)];

normData_x = normData_x;
normData_u = normData_u;

vals = [];
xapn = [0 0 0];
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
    
    z0 = g_z([0 0 0], stateEncoderWeight_cell, stateEncoderBias_cell)';
    
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
        
       [Qlqr, xi] = getQz(xapn, ...
           stateEncoderWeight_cell, stateEncoderBias_cell, ...
           Q, Kz, Kw, 5)
        
    else
       hasLinearState = true;
       
       Qlqr = [ Q zeros(3, size(Kz,2)-3);
                zeros(size(Kz,1)-3, size(Kz,2))];
        z0 = [[0 0 0]'; z0];
    end

    if contains(runs{i}, 'nonlinearInput')
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
       
        [Rlqr, selectedInputs, identityError_u, identityErrorMatrix_u] = ...
            getQz([0 0], ...
            inputEncoderWeight_cell, inputEncoderBias_cell, ...
            inputDecoderWeight_cell, inputDecoderBias_cell, ...
            R, 1);
        
        w0 = g_z(0, inputEncoderWeight_cell, inputEncoderBias_cell)';
          
    else
        hasNonlinearInput = false;
        Rlqr = R;
        w0 = 0;
    end
    
    C = dlqr(Kz, Kw, Qlqr, Rlqr);
    
    if contains(runs{i}, 'nonlinearInput')
        cont = @(x) g_z(-C * g_z(x, stateEncoderWeight_cell, stateEncoderBias_cell)', ...
            inputDecoderWeight_cell, inputDecoderBias_cell);
    else
        cont = @(x) -C * g_z(x, stateEncoderWeight_cell, stateEncoderBias_cell)';
    end

    mat = readmatrix(data_path);

    % without controller
    [sub, alpha, gamma] = getStabilitConditions(Kz, eye(size(Kz,1)), ...
        mat, stateEncoderWeight_cell, stateEncoderBias_cell, 300);
    vals = [vals, [sub; alpha; gamma]];

    [val, val1, val2] = ...
        getKoopmanLyapunovAttractor_measured(Kz, stateEncoderWeight_cell, stateEncoderBias_cell, ...
        mat, alpha);
    
    dsfac = 1;
        
    % get grid map from data
    mat_t = mat(1:end-1,:)*3.*sqrt(var_data_ds(3:end))+mean_data_ds(3:end);
    index = 1:dsfac:size(mat_t(val,:),1);
    mat_t_s = mat_t;
    mat_t_s(~val,:) = -10;
    
    figure; 
    scatter3(mat_t(index,1),mat_t(index,2),mat_t(index,3),'x', 'linewidth',1.5)
    hold on
    
    scatter3(mat_t_s(index,1),mat_t_s(index,2),mat_t_s(index,3),'o', 'linewidth',1.5)

    [val, val1, val2] = ...
        getKoopmanLyapunovAttractor_measured(Kz-Kw*C, stateEncoderWeight_cell, stateEncoderBias_cell, ...
        mat, alpha);

    mat_t_s = mat_t;
    mat_t_s(~val,:) = -10;
    
    scatter3(mat_t_s(index,1),mat_t_s(index,2),mat_t_s(index,3),'^', 'linewidth',1.5)
    xlabel('x1')
    ylabel('x2')
    zlabel('x3')
    legend('Data', 'uncontrolled stable', 'controlled stable')
    
    xlim([0 400])
    ylim([0 200])
    zlim([0 400])
end







