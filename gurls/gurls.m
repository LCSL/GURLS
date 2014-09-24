% The GURLS Package
% Copyright (c) 2011, P. Mallapragada, A. Tacchetti
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 
%     * Redistributions of source code must retain the above
%       copyright notice, this list of conditions and the following
%       disclaimer.
%     * Redistributions in binary form must reproduce the above
%       copyright notice, this list of conditions and the following
%       disclaimer in the documentation and/or other materials
%       provided with the distribution.
%     * Neither the name(s) of the copyright holders nor the names
%       of its contributors or of the Massacusetts Institute of Technology
%		may be used to endorse or promote products
%       derived from this software without specific prior written
%       permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [opt] = gurls(X, y, opt, jobid)
	


%% Enums
IGN = 0; % Ignore
CPT = 1; % Compute
CSV = 2; % Compute and Save
LDF = 3; % Load from file
DEL = 4; % Explicitly Delete
%% End Enums

if jobid == 0
	return;
end

if ~isa(opt, 'GurlsOptions')
    opt = GurlsOptions(opt);
end

for i = 1:numel(opt.process)
	if numel(opt.process{i}) ~= numel(opt.seq)
		error('Number of elements in process and sequence not same\n');
	end
end

%setup(); % make this more general

seq = opt.seq;
% Load and copy

if exist(opt.savefile) == 2
	t = load(opt.savefile);
	if isprop(t.opt,'time');
		opt.time = t.opt.time;
	end	
else
    if ~isequal(opt.name,'')
        fprintf('Could not load %s. Starting from scratch.\n', opt.savefile);
    end
end	
%try
%	t = load(opt.savefile);
%	if isfield(t.opt,'time')
%		opt.time = t.opt.time;
%	end
%catch
%	fprintf('Could not load %s. Starting from scratch.\n', opt.savefile);
%end

process = opt.process{jobid};


%for i = 1:numel(opt.process) % Go by the length of process.
opt.time{jobid} = struct;
%end

fprintf('\nNew job sequence...\n');
for i = 1:numel(process) % Go by the length of process.
	reg = regexp(seq{i},':','split');
	if length(reg) < 2
		error('Command format incorrect. Requires a ":" \n');
	end
	fprintf('[Job %d: %15s]: %15s ',jobid, reg{1}, reg{2});

	switch process(i)
	case IGN
		fprintf('\tignored\n');
		continue;

	case {CPT, CSV, ~isfield(opt,reg{1})}

		fName = [reg{1} '_' reg{2}];
		fun = str2func(fName);
		tic;
        opt.newprop(reg{1}, fun(X,y,opt));
		opt.time{jobid} = setfield(opt.time{jobid},reg{1}, toc);
		fprintf('\tdone\n');

	case LDF,
		if exist('t','var') && isprop(opt, reg{1})
			opt.reg{1} = t.opt.reg{1};
			fprintf('\tcopied\n');
		else
			fprintf('\tcopy failed\n');
		end

	otherwise
		fprintf('Unknown process statement\n');
	end	
end

if ~isequal(opt.name, '')
    fprintf('\nSave cycle...\n');
    % Delete whats not necessary
    for i = 1:numel(process)
        reg = regexp(seq{i},':','split');
        fprintf('[Job %d: %15s] %15s: ',jobid, reg{1}, reg{2});
        switch process(i)
            case {CSV, LDF}
                fprintf('\tsaving..\n');
            otherwise
                if isprop(opt, reg{1})
                    opt.(reg{1}) = [];
                    fprintf('\tremoving..\n');
                else
                    fprintf('\tnot found..\n');
                end
        end
    end
    save(opt.savefile, 'opt', '-v7.3');
    fprintf('Saving opt in %s\n', opt.savefile);
end
