function DeleteXlsSheet1(excelFileName)
% Function to delete the 1st sheet of a spreadsheet. Allows for custom page
% names without the extra sheet.

% excelFileName = 'C:\Users\Tim\Documents\Academic - Research\Data\Emily_drops\Experiment_summary_auto.xlsx';
excelFilePath = pwd; % Current working directory.
sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
% Open Excel file.
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(fullfile(excelFilePath, excelFileName)); % Full path is necessary!
% objExcel.Workbooks.Open(excelFileName); % Full path is necessary!
% Delete sheets.
try
      % Throws an error if the sheets do not exist.
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
catch
      ; % Do nothing.
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;