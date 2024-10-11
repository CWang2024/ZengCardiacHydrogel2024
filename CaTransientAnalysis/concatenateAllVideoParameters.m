function concatenateAllVideoParameters(folderPath, outputFileName, sheetName)

    excelFiles = dir(fullfile(folderPath, '*.xlsx'));

    excelPath = fullfile(folderPath, outputFileName);
    excelFiles = excelFiles(~strcmp({excelFiles.name}, outputFileName));
    % Loop through each Excel file
    colNames = readcell(fullfile(folderPath, excelFiles(1).name), 'Sheet', sheetName, 'Range', 'A:A');
    
    writecell(colNames, excelPath, 'Sheet', sheetName, 'Range', 'A1');
    for i = 1:length(excelFiles)
        % Get the full path of the current Excel file
        filePath = fullfile(folderPath, excelFiles(i).name);
        
        % Read the data from Sheet 1, Column 2 (starting from B2)
        data = readmatrix(filePath, 'Sheet', sheetName, 'Range', 'B:B'); % Adjust range if necessary
        
        % Remove any NaN values (in case of empty cells)
        data = data(~isnan(data));
        
        % Determine the column letter for the current file
        colLetter = char('A' + i);  % Converts index to letter (A, B, C, ...)
        header = excelFiles(i).name(1:end-5);  % Remove .xlsx extension
        writecell({header}, excelPath, 'Sheet', sheetName, 'Range', [colLetter, '1']);
        % Write the data to the corresponding column in the new file
        writematrix(data, excelPath, 'Sheet', sheetName, 'Range', [colLetter, '2']);
    end
    




 