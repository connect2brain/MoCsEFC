function [subjData] = readNeuroneFromTable(dataPath, row)
%READNEURONEFROMTABLE Reads the neurone-data specified by the dataPath (a
% path), and the table-row (with at least the fields neuroneFile [a cell
% array containing a string] and neuroneIndex [a natural number])
% returns the raw neurone-data-structure
session_file = row.neuroneFile{:};
neurone_idx  = row.neuroneIndex;
file = sprintf('%s\\%s\\NeurOne-%s.ses', dataPath, session_file, session_file);
warning('off')
subjData = module_read_neurone(file, neurone_idx);
warning('on')
end

