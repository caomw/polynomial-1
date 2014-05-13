energy_dc = cell(1, 100);
time_dc = cell(1, 100);
y_dc = cell(1, 100);
mean_energy_dc = 0;
mean_time_dc = 0;
for i = 1:100
    load (['res_',num2str(i),'_noise_2_admm_100_rule_0']);
    time_dc{i} = (1:length(primal))/length(primal)*runtime;
    energy_dc{i} = primal;
    y_dc{i} = y;
    mean_energy_dc = mean_energy_dc + primal(end);
    mean_time_dc = mean_time_dc + runtime;
end
mean_time_dc = mean_time_dc / 100
mean_energy_dc = mean_energy_dc / 100


energy_cccp = cell(1, 100);
time_cccp = cell(1, 100);
y_cccp = cell(1, 100);
mean_energy_cccp = 0;
mean_time_cccp = 0;
for i = 1 : 100
    fileID = fopen(['./150/outputCCCP', num2str(i),'.txt'],'r');
    formatSpec = '%f ';
    sizeA = [1 Inf];
    datalist = fscanf(fileID,formatSpec,sizeA);
    y_cccp{i} = datalist(1:9);
    datalist = reshape(datalist(10:end), [], 2)';
    time_cccp{i} = datalist(2,:);
    energy_cccp{i} = datalist(1,:);
    
    fileID = fopen(['./140/outputCCCP', num2str(i),'.txt'],'r');
    formatSpec = '%f ';
    sizeA = [1 Inf];
    datalist = fscanf(fileID,formatSpec,sizeA);
    datalist = reshape(datalist(10:end), [], 2)';
    temp_y = datalist(1:9);
    temp_time = datalist(2,:);
    temp_energy = datalist(1,:);
    if (temp_energy(end) < energy_cccp{i}(end))
        time_cccp{i} = temp_time;
        energy_cccp{i} = temp_energy;
        y_cccp{i} = temp_y;
    end
    mean_energy_cccp = mean_energy_cccp + energy_cccp{i}(end);
    mean_time_cccp = mean_time_cccp + time_cccp{i}(end);
    fclose(fileID);
end
mean_time_cccp = mean_time_cccp / 100
mean_energy_cccp = mean_energy_cccp / 100


energy_bfgs = cell(1, 100);
time_bfgs = cell(1, 100);
y_bfgs = cell(1, 100);
mean_energy_bfgs = 0;
mean_time_bfgs = 0;
for i = 1 : 100
    fileID = fopen(['./150/outputRegistration',num2str(i),'.txt'],'r');
    formatSpec = '%f ';
    sizeA = [1 1];
    energy_bfgs{i} = fscanf(fileID,formatSpec,sizeA);
    time_bfgs{i} = fscanf(fileID,formatSpec,sizeA);
    y_bfgs{i} = fscanf(fileID,formatSpec,[1 9]);
    mean_energy_bfgs = mean_energy_bfgs + energy_bfgs{i};
    mean_time_bfgs = mean_time_bfgs + time_bfgs{i};
    fclose(fileID);
end
mean_time_bfgs = mean_time_bfgs / 100
mean_energy_bfgs = mean_energy_bfgs / 100

for i = 1 : 100
    relative_energy(i) = min(energy_cccp{i}) - min(energy_dc{i});
    relative_energy_bfgs(i) = min(energy_bfgs{i}) - min(energy_dc{i});
end
subplot(1, 2, 1);plot(relative_energy);hold on;subplot(1, 2, 2);plot(relative_energy_bfgs, 'r');
