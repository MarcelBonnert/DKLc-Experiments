% path to raw data
path = '';

files = dir(path);
files = files(5:end);
load('filters.mat', 'Cheby_Filter_hf_2_Ts_1000Hz')

%% Get experimental equilibirum
for i = 1:length(files)
    % Get dataset
    load([path '/', files(i).name])
    dataSet = eval([files(i).name(1:end-4), ';']);
    data = [dataSet.Y(1).Data', dataSet.Y(2).Data', ...
            dataSet.Y(3).Data', dataSet.Y(4).Data', dataSet.Y(5).Data'];
    eval(['clear ' files(i).name(1:end-4)])
    mi(i,:) = mean(data, 1);
    Li(i) = size(data, 1);
end

equi = sum(mi .* Li')/sum(Li);


%% 
chop_fac = 50;
Ts_target = 5;
Ts = 0.001;

for i = 1:length(files)
    % Get dataset
    load([path '/', files(i).name])
    dataSet = eval([files(i).name(1:end-4), ';']);
    
    % Get beginning
    d = diff(dataSet.Y(1).Data);
    start_indicies = find(d ~= 0);
    Start_index = start_indicies(1);

    % Get ending
    end_trigger = true;
    start_trigger = false;
    buffer_trigger = false;
    index = 1;
    while end_trigger 
        d = diff(dataSet.Y(1).Data(end-index:end));
        if start_trigger
            if d(1) == 0
                end_trigger = false;
            else
                buffer_trigger = false;
                start_trigger = false;
            end
        elseif buffer_trigger
            if d(1) == 0
                start_trigger = true;
            end
        else
            if any(d ~= 0)
                buffer_trigger = true;
            end
        end
        index = index + 1;
    end

    End_index = index - 3;

    data = [dataSet.Y(1).Data', dataSet.Y(2).Data', ...
            dataSet.Y(3).Data', dataSet.Y(4).Data', dataSet.Y(5).Data'];
    
    % Filter
    data_filt = filter(Cheby_Filter_hf_2_Ts_1000Hz, data-equi) + equi;
    clear data
    
    % Chop data
    data_final = data_filt(Start_index:end-End_index, :);
    clear data_filt
    
    % Down sampling
    t = dataSet.X.Data;

    t = t(Start_index:end-End_index)';
    t = t - t(1);

    index = 1:Ts_target/Ts:length(t);
    t_ds = t(index);
    data_ds = data_final(index,:);
    clear data_final
    
    % Chop at multiples of cop_fac
    number_traj = floor(size(data_ds,1)/chop_fac);
    data_chop = data_ds(2:number_traj * chop_fac + 1, :);
    clear data_ds
    
    final_data{i} = data_chop;
end

%% Normalize data and save it
data = [];
for i = 1:length(final_data)
    data = [data; final_data{i}];
end

lim_min = min(data);
lim_max = max(data);
var_data = var(data);

data_n = (data - (lim_min + lim_max)/2) ./ ...
    (3 * sqrt(var_data));

sep = [80, 15, 5];
L = size(data_n,1);
L_train = sep(1)/100*L - mod(sep(1)/100*L, chop_fac);
L_val = sep(2)/100*L - mod(sep(2)/100*L, chop_fac);
L_test = sep(3)/100*L - mod(sep(3)/100*L, chop_fac);

data_train = data_n(1:L_train,:);
data_val = data_n(1 + L_train:L_train + L_val,:);
data_test = data_n(1 + L_train + L_val:end,:);

number_traj = floor(size(data_n,1)/chop_fac);
T = [];
t = 0:Ts_target:(chop_fac - 1) * Ts_target;
for i = 1:number_traj
    T = [T;t'];
end

writematrix(data_train(:,3:end), 'preprocessedData/threeTankTrainingStates.csv')
writematrix(data_train(:,1:2), 'preprocessedData/threeTankTrainingInputs.csv')
writematrix(T(1:L_train), 'preprocessedData/threeTankTrainingTime.csv')

writematrix(data_val(:,3:end), 'preprocessedData/threeTankValidationStates.csv')
writematrix(data_val(:,1:2), 'preprocessedData/threeTankValidationInputs.csv')
writematrix(T(1:L_val), 'preprocessedData/threeTankValidationTime.csv')

writematrix(data_test(:,3:end), 'preprocessedData/threeTankTestStates.csv')
writematrix(data_test(:,1:2), 'preprocessedData/threeTankTestInputs.csv')
writematrix(T(1:L_test), 'preprocessedData/threeTankTestTime.csv')




