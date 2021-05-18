function chk = testfun(fun)
% Test function for NATSORT, NATSORTFILES, and NATSORTROWS. Do NOT call!
%
% (c) 2012-2020 Stephen Cobeldick
%
% See also NATSORT_TEST NATSORTFILES_TEST NATSORTROWS_TEST
tmp = {'%s %3d','<a href="matlab:opentoline(''%1$s'',%2$3d)">%1$s %2$3d</a>'};
fmt = tmp{1+feature('hotlinks')};
itr = 0;
cnt = 0;
chk = @nestfun;
	function nestfun(varargin) % (function inputs..., fun, expected function outputs...)
		dbs = dbstack();
		%
		if ~nargin % post-processing
			fprintf('%s: %d of %d testcases failed.\n',dbs(2).file,cnt,itr)
			return
		end
		%
		itr = itr+1;
		%
		idx = find(cellfun(@(f)isequal(f,fun),varargin));
		assert(nnz(idx)==1,'SC:testfun:MissFun','Missing function handle.')
		xfa = varargin(idx+1:end);
		ofa = cell(size(xfa));
		%
		[ofa{:}] = fun(varargin{1:idx-1});
		%
		boo = false;
		for k = 1:numel(xfa)
			if ~isequal(xfa{k},[])
				if ~isequal(size(ofa{k}),size(xfa{k}))
					boo = true;
					otx = sprintf('x%d',size(ofa{k}));
					xtx = sprintf('x%d',size(xfa{k}));
				elseif ~isequaln(ofa{k},xfa{k})
					boo = true;
					oar = ofa{k}.';
					xar = xfa{k}.';
					switch k
						case 1 % cell array of char vectors
							otx = sprintf(', ''%s''',oar{:});
							xtx = sprintf(', ''%s''',xar{:});
						case 2 % indices
							otx = sprintf(', %d',oar);
							xtx = sprintf(', %d',xar);
						case 3 % debug cell array
							obo = cellfun(@ischar,oar);
							xbo = cellfun(@ischar,xar);
							ocs = cellfun(@num2str,oar,'uni',0);
							xcs = cellfun(@num2str,xar,'uni',0);
							ocs(~obo) = strcat('[',ocs(~obo),']');
							xcs(~xbo) = strcat('[',xcs(~xbo),']');
							ocs(obo) = strcat('''',ocs(obo),'''');
							xcs(xbo) = strcat('''',xcs(xbo),'''');
							otx = sprintf(', %s',ocs{:});
							xtx = sprintf(', %s',xcs{:});
					end
				end
			end
			if boo
				fprintf(fmt, dbs(2).file, dbs(2).line);
				fprintf(' (output argument %d)',k);
				fprintf('\noutput:%s\nexpect:%s\n', otx(2:end), xtx(2:end));
			end
		end
		cnt = cnt+boo;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%testfun