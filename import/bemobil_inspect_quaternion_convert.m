% bemobil_inspect_quaternion_convert visualizes all 24 possible permutations when transforming a quaternion dataset into
% euler angles. It takes in a freshly loaded stream from load_xdf. If no other info is provided, it just displays all
% channel names so the user can determine the initial indices. If indices are provided, all permutations are plotted. If
% necessary, the angles can also be unwrapped. If they are unwrapped, the axes are also linked for inspection.
%
% example usage for the import:
%           s = load_xdf(config.filename);
%           stream = s{5} % or whatever index you want to inspect
%           bemobil_inspect_quaternion_convert(stream);
%           quaternionIndices = [4:7];
%           bemobil_inspect_quaternion_convert(stream, quaternionIndices);
%
% inputs: 
%           stream                  - one stream from the raw xdf file
%           quaternionIndices       - vector of 4 indices that are the quaternion channels    
%           do_unwrap               - boolean whether or not unwrap should be applied (default = 0)
% 
% 
% outputs: 
%           plothandles             - handles to all created plots, e.g. to save them on disc and close automatically
% 
% see also: util_quat2eul, unwrap

function plothandles = bemobil_inspect_quaternion_convert(stream, quaternionIndices, do_unwrap)

if nargin == 0
    help bemobil_inspect_quaternion_convert
    return
end

if ~exist('quaternionIndices','var')
    for i_chan = 1:length(stream.info.desc.channels.channel)
        disp(stream.info.desc.channels.channel{i_chan}.label)
    end
    return
end
if ~exist('do_unwrap','var')
    do_unwrap = 0;
end


allperms = perms(quaternionIndices);

for i_perm = 1:24
    
    orientationInEuler         = util_quat2eul(stream.time_series(allperms(i_perm,:),:)'); % the BeMoBIL util script
    
    if do_unwrap
        orientationInEuler  = transpose(unwrap(orientationInEuler', [], 2));
    end
    
    plothandles(i_perm) = figure('color','w','Position', [100 200 1400 800]);
    ax(1) = subplot(3,1,1);
    plot(orientationInEuler(:,1))
    title({['permutation of stream ' stream.info.name ' with quaternion channels: ' num2str(allperms(i_perm,:))]
        [stream.info.desc.channels.channel{allperms(i_perm,1)}.label ', ' stream.info.desc.channels.channel{allperms(i_perm,2)}.label ', '...
        stream.info.desc.channels.channel{allperms(i_perm,3)}.label ', ' stream.info.desc.channels.channel{allperms(i_perm,4)}.label]},...
        'interpreter','none')
    ax(2) = subplot(3,1,2);
    plot(orientationInEuler(:,2))
    ax(3) = subplot(3,1,3);
    plot(orientationInEuler(:,3))
    if ~do_unwrap;linkaxes(ax);end
    
end
drawnow