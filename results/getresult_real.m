energy_cccp = cell(1, 100);
time_cccp = cell(1, 100);
y_cccp = cell(1, 100);
mean_energy_cccp = 0;
mean_time_cccp = 0;
for i = 1 : 100
    fileID = fopen(['./150nowait/outputCCCP', num2str(i),'.txt'],'r');
    formatSpec = '%f ';
    sizeA = [1 Inf];
    datalist = fscanf(fileID,formatSpec,sizeA);
    y_cccp{i} = datalist(1:81);
    datalist = reshape(datalist(82:end), [], 2)';
    time_cccp{i} = datalist(2,:);
    energy_cccp{i} = datalist(1,:);
    mean_energy_cccp = mean_energy_cccp + energy_cccp{i}(end);
    mean_time_cccp = mean_time_cccp + time_cccp{i}(end);
    fclose(fileID);
end
mean_time_cccp = mean_time_cccp / 100
mean_energy_cccp = mean_energy_cccp / 100
for i = 1 : 100
    relative_energy(i) = min(energy_cccp{i});
end
plot(relative_energy);