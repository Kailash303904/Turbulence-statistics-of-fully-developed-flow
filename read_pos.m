function pos = read_pos(filename)

opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = 2;
opts.Delimiter = " ";
opts.VariableNames = ["i", "xh"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
pos = readtable(filename, opts);

end