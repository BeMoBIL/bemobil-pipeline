function [X,ndx,dbg] = natsortfiles(X,rgx,varargin)
% Alphanumeric / Natural-Order sort a cell array of filenames/filepaths (1xN char).
%
% (c) 2014-2020 Stephen Cobeldick
%
% Alphanumeric sort of a cell array of filenames or foldernames, i.e. sort
% the names by character order and by the values of any numbers within the
% names. The file/folder names, file extensions, and path directories (if
% supplied) are sorted separately: this ensures that shorter names sort
% before longer ones. For names without file extensions (e.g. foldernames,
% filenames without extensions) the 'noext' option should be used.
%
%%% Example:
% P = 'C:\SomeDir\SubDir';
% S = dir(fullfile(P,'*.txt'));
% C = natsortfiles({S.name});
% for k = 1:numel(C)
%     fullfile(P,C{k})
% end
%
%%% Syntax:
%  Y = natsortfiles(X)
%  Y = natsortfiles(X,rgx)
%  Y = natsortfiles(X,rgx,<options>)
% [Y,ndx,dbg] = natsortfiles(X,...)
%
% To sort all of the strings in a cell array use NATSORT (File Exchange 34464).
% To sort the rows of a cell array of strings use NATSORTROWS (File Exchange 47433).
%
%% File Dependency %%
%
% NATSORTFILES requires the function NATSORT (File Exchange 34464). The optional
% arguments <options> are passed directly to NATSORT. See NATSORT for case
% sensitivity, sort direction, number substring matching, and other options.
%
%% Explanation %%
%
% Using SORT on filenames will sort any of char(0:45), including the printing
% characters ' !"#$%&''()*+,-', before the file extension separator character '.'.
% Therefore this function splits the name and extension and sorts them separately.
%
% Similarly the file separator character within filepaths can cause longer
% directory names to sort before shorter ones, as char(0:46)<'/' and
% char(0:91)<'\'. Check this example to see why this matters:
%
% >> X = {'A1\B', 'A+/B', 'A/B1', 'A=/B', 'A\B0'};
% >> sort(X)
% ans =   'A+/B'  'A/B1'  'A1\B'  'A=/B'  'A\B0'
% >> natsortfiles(X)
% ans =   'A\B0'  'A/B1'  'A1\B'  'A+/B'  'A=/B'
%
% NATSORTFILES splits filepaths at each file separator character and sorts
% every level of the directory hierarchy separately, ensuring that shorter
% directory names sort before longer, regardless of the characters in the names.
%
%% Examples %%
%
% >> A = {'a2.txt', 'a10.txt', 'a1.txt'};
% >> sort(A)
% ans = 'a1.txt'  'a10.txt'  'a2.txt'
% >> natsortfiles(A)
% ans = 'a1.txt'  'a2.txt'  'a10.txt'
%
% >> B = {'test_new.m'; 'test-old.m'; 'test.m'};
% >> sort(B) % Note '-' sorts before '.':
% ans =
%    'test-old.m'
%    'test.m'
%    'test_new.m'
% >> natsortfiles(B) % Shorter names before longer:
% ans =
%    'test.m'
%    'test-old.m'
%    'test_new.m'
%
% >> C = {'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'};
% >> sort(C) % Wrong number order:
% ans =
%    'test.m'
%    'test1.m'
%    'test10-old.m'
%    'test10.m'
%    'test2.m'
% >> natsortfiles(C) % Shorter names before longer:
% ans =
%    'test.m'
%    'test1.m'
%    'test2.m'
%    'test10.m'
%    'test10-old.m'
%
%%% Directory Names:
% >> D = {'A2-old\test.m';'A10\test.m';'A2\test.m';'A1archive.zip';'A1\test.m'};
% >> sort(D) % Wrong number order, and '-' sorts before '\':
% ans =
%    'A10\test.m'
%    'A1\test.m'
%    'A1archive.zip'
%    'A2-old\test.m'
%    'A2\test.m'
% >> natsortfiles(D) % Shorter names before longer:
% ans =
%    'A1archive.zip'
%    'A1\test.m'
%    'A2\test.m'
%    'A2-old\test.m'
%    'A10\test.m'
%
%% Input and Output Arguments %%
%
%%% Inputs (**=default):
% X   = CellArrayOfCharRowVectors, with filenames or filepaths to be sorted.
% rgx = Regular expression to match number substrings, '\d+'**
%     = [] uses the default regular expression, which matches integers.
% <options> can be supplied in any order:
%     = 'noext' should be specified for names without extensions.
%     = all remaining options are passed directly to NATSORT.
%
%%% Outputs:
% Y   = CellArrayOfCharRowVectors, filenames of <X> sorted alphanumerically.
% ndx = NumericMatrix, same size as <X>. Indices such that Y = X(ndx).
% dbg = CellVectorOfCellArrays, size 1xMAX(2+NumberOfDirectoryLevels).
%       Each cell contains the debug cell array for directory names, filenames,
%       and file extensions. Helps debug the regular expression. See NATSORT.
%
% See also SORT NATSORT NATSORTROWS DIR FILEPARTS FULLFILE NEXTNAME CELLSTR REGEXP IREGEXP SSCANF

%% Input Wrangling %%
%
assert(iscell(X),...
	'SC:natsortfiles:X:NotCellArray',...
	'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
assert(all(tmp(:)),...
	'SC:natsortfiles:X:ContentNotCharVectors',...
	'First input <X> must be a cell array of strings (1xN character).')
%
ide = cellfun(@(s)ischar(s)&&strcmpi(s,'noext'),varargin);
assert(nnz(ide)<2,...
	'SC:natsortfiles:options:noextOverspecified',...
	'File-extension handling is overspecified.')
varargin(ide) = [];
%
if nargin>1
	varargin = [{rgx},varargin];
end
%
%% Split and Sort File Names/Paths %%
%
% Split full filepaths into file [path,name,extension]:
[pth,fnm,ext] = cellfun(@fileparts,X(:),'UniformOutput',false);
% Split path into {dir,subdir,subsubdir,...}:
pth = regexp(pth,'[^/\\]+','match'); % either / or \ as filesep.
len = cellfun('length',pth);
num = max(len);
vec = cell(numel(len),1);
%
if any(ide)
	fnm = strcat(fnm,ext);
	ext(:) = {''};
end
%
% Alphanumeric sort of the file extensions and filenames:
if isempty(num)
	ndx = [];
	ids = [];
	dbg = {};
elseif nargout<3 % faster:
	[~,ndx] = natsort(ext,varargin{:});
	[~,ids] = natsort(fnm(ndx),varargin{:});
else % for debugging:
	[~,ndx,dbg{num+2}] = natsort(ext,varargin{:});
	[~,ids,tmp] = natsort(fnm(ndx),varargin{:});
	[~,idd] = sort(ndx);
	dbg{num+1} = tmp(idd,:);
end
ndx = ndx(ids);
%
% Alphanumeric sort of the directory names:
for k = num:-1:1
	idx = len>=k;
	vec(:) = {''};
	vec(idx) = cellfun(@(c)c(k),pth(idx));
	if nargout<3 % faster:
		[~,ids] = natsort(vec(ndx),varargin{:});
	else % for debugging:
		[~,ids,tmp] = natsort(vec(ndx),varargin{:});
		[~,idd] = sort(ndx);
		dbg{k} = tmp(idd,:);
	end
	ndx = ndx(ids);
end
%
% Return the sorted array and indices:
ndx = reshape(ndx,size(X));
X = X(ndx);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortfiles