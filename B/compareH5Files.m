function compareH5Files(filename1, filename2)
% Function to compare if two HDF5 files are identical. filename1 is used as
% the reference.
% author: Bradley Treeby
% date: 27th August 2012

% get the contents of files
info1 = h5info(filename1);
info2 = h5info(filename2);

% get the variables names
var_names1 = {info1.Datasets.Name};
var_names2 = {info2.Datasets.Name};

% use the file with more datasets as reference
if length(var_names1) >= length(var_names2)
    var_names = var_names1;
else
    var_names = var_names2;
    disp('WARNING: The second file has more datasets than the first, using this file as reference...');
end

% set ok flag
ok_flag = 0;

disp('---------------------------------------------------------');
disp('Comparing Datasets');
disp('---------------------------------------------------------'); 

% loop through each of the variables, then check the contents
for var_index = 1:length(var_names)
    
    % load reference data
    var1 = h5read(filename1, ['/' var_names{var_index}]);
    
    try
        % load comparison data
        var2 = h5read(filename2, ['/' var_names{var_index}]);
        
        % check values
        if any(size(var1) ~= size(var2))
            disp(['WARNING: Sizes for ' var_names{var_index} ' are different']);
            ok_flag = ok_flag + 1;
        elseif sum(var1(:) - var2(:))
            mx_error = max(abs(var1(:) - var2(:)));
            if  mx_error < 10*eps
                disp(['Contents for ' var_names{var_index} ' are within 10^-15']);
            else
                disp(['WARNING: Contents for ' var_names{var_index} ' are different (L_inf = ' num2str(mx_error) ')']);
                ok_flag = ok_flag + 1;
            end
        else
            disp(['Contents for ' var_names{var_index} ' are equal']);
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:imagesci:h5read:datasetDoesNotExist')
            disp(['WARNING: Contents for ' var_names{var_index} ' is missing in second file']);
            ok_flag = ok_flag + 1;
        else
            rethrow(ME);
        end
        
    end
end

disp('---------------------------------------------------------');
disp('Comparing Attributes');
disp('---------------------------------------------------------'); 

% get the variables names
att_names = {info1.Attributes.Name};

% loop through each of the attributes, then check the contents
for att_index = 1:length(att_names)
    
    % load reference data
    att1 = h5readatt(filename1, '/', att_names{att_index});
    
    try
        % load comparison data
        att2 = h5readatt(filename2, '/', att_names{att_index});
        
        % check values
        if ~strcmp(att1, att2)
            disp(['WARNING: Contents for ' att_names{att_index} ' are different']);
            ok_flag = ok_flag + 1;
        else
            disp(['Contents for ' att_names{att_index} ' are equal']);
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:imagesci:h5readatt:cannotOpenAttribute')
            disp(['WARNING: Contents for ' att_names{att_index} ' is missing in second file']);
            ok_flag = ok_flag + 1;
        else
            rethrow(ME);
        end
        
    end
end

% display summary
disp('---------------------------------------------------------');
disp(['Total of ' num2str(ok_flag) ' differences between the files']);
disp('---------------------------------------------------------');    