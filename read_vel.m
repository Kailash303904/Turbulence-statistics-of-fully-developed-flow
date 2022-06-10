function vel = read_vel(filename)

startRow = 2;
formatSpec = '%3C%4C%4C%14f%13f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

fclose(fileID);
vel = table(dataArray{1:end-1}, 'VariableNames', {'i','j','k','U','V','W'});
clearvars filename startRow formatSpec fileID dataArray ans;

end

