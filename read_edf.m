
%the current file read edf was made with the help of AI
% source PSG data could not be uploaded to github due to file size constraint
function [data, header] = read_edf(filename)
    % EDF reader for PSG data
    % Usage: [data, header] = read_edf('path/to/file.edf')
    
    fprintf('Opening file: %s\n', filename);
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Read header (256 bytes)
    header.version = deblank(fread(fid, 8, '*char')');
    header.patient = deblank(fread(fid, 80, '*char')');
    header.recording = deblank(fread(fid, 80, '*char')');
    header.startdate = deblank(fread(fid, 8, '*char')');
    header.starttime = deblank(fread(fid, 8, '*char')');
    header.bytes = str2double(fread(fid, 8, '*char')');
    fread(fid, 44, '*char');  % reserved
    header.records = str2double(fread(fid, 8, '*char')');
    header.duration = str2double(fread(fid, 8, '*char')');
    header.channels = str2double(fread(fid, 4, '*char')');
    
    fprintf('File contains %d channels, %d records\n', header.channels, header.records);
    
    % Read names of channels used
    for i = 1:header.channels
        header.label{i} = deblank(fread(fid, 16, '*char')');
    end
    
    % Skip transducer type- not requried for this analysis
    fread(fid, 80 * header.channels, '*char');
    
    % Read units
    for i = 1:header.channels
        header.units{i} = deblank(fread(fid, 8, '*char')');
    end
    
    % Read physical min/max
    for i = 1:header.channels
        header.physicalMin(i) = str2double(fread(fid, 8, '*char')');
    end
    for i = 1:header.channels
        header.physicalMax(i) = str2double(fread(fid, 8, '*char')');
    end
    
    % Read digital min/max
    for i = 1:header.channels
        header.digitalMin(i) = str2double(fread(fid, 8, '*char')');
    end
    for i = 1:header.channels
        header.digitalMax(i) = str2double(fread(fid, 8, '*char')');
    end
    
    % Skip prefiltering
    fread(fid, 80 * header.channels, '*char');
    
    % Read samples per record
    for i = 1:header.channels
        header.samples(i) = str2double(fread(fid, 8, '*char')');
    end
    
    % Skip reserved
    fread(fid, 32 * header.channels, '*char');
    
    % Calculate sampling rate
    header.fs = header.samples(1) / header.duration;
    
    fprintf('Sampling rate: %d Hz\n', header.fs);
    fprintf('Duration: %.1f seconds\n', header.records * header.duration);
    
    % Read data 
    records_to_read = min(100, header.records);
    fprintf('Reading first %d records (out of %d total)...\n', records_to_read, header.records);
    
    total_samples = records_to_read * header.samples(1);
    data = zeros(header.channels, total_samples);
    
    for rec = 1:records_to_read
        if mod(rec, 10) == 0
            fprintf('  Reading record %d/%d\n', rec, records_to_read);
        end
        
        for ch = 1:header.channels
            % Read raw values
            raw = fread(fid, header.samples(ch), 'int16');
            
            % Convert to real units
            gain = (header.physicalMax(ch) - header.physicalMin(ch)) / ...
                   (header.digitalMax(ch) - header.digitalMin(ch));
            offset = header.physicalMin(ch) - gain * header.digitalMin(ch);
            
            start_idx = (rec-1) * header.samples(ch) + 1;
            end_idx = rec * header.samples(ch);
            data(ch, start_idx:end_idx) = gain * double(raw) + offset;
        end
    end
    
    fclose(fid);
    
    fprintf('Successfully loaded: %d channels x %d samples\n', size(data,1), size(data,2));
    fprintf('Memory used: %.1f MB\n', (numel(data) * 8) / (1024^2));

end




