% bemobil_check_config - Checks if all fields of the config are present. If fields are missing, they will be shown, and
% in case a default for this field exists it will be entered. If no default exists, the filed is considered mandatory
% and an error is thrown.
%
% Usage:
%   >>  [bemobil_config] = bemobil_check_config(bemobil_config)
% 
% Inputs:
%   bemobil_config              - current bemobil_config struct
%
% Outputs:
%   bemobil_config              - checked bemobil_config struct with default filled if missing
%
% See also:
%   example_bemobil_config
%
% Authors: Marius Klug, 2021

function bemobil_config = bemobil_check_config(bemobil_config)

return %for now!

if nargin == 0
    help bemobil_check_config
    return
end


%% load defaults

load default_config
load full_example_config

default_config_fields = fieldnames(default_config);
example_config_fields = fieldnames(full_example_config);

%% check fields

disp('Checking config fields...')

for i_field = 1:length(example_config_fields)
    
    if ~isfield(bemobil_config,example_config_fields{i_field})
        warning(['Field not present: ''' example_config_fields{i_field} ''''])
        
        if isfield(default_config,example_config_fields{i_field})
            disp('Using default:')
            
            default_idx = find(~cellfun(@isempty,strfind(default_config_fields,example_config_fields{i_field})));
            default_config.(default_config_fields{default_idx})
            bemobil_config.(default_config_fields{default_idx}) = default_config.(default_config_fields{default_idx});
        else
            error(['Mandatory field missing: bemobil_config.' example_config_fields{i_field}])
        end
    end
end

disp('All fields checked!')