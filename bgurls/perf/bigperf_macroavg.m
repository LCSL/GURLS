function [p] = bigperf_macroavg(X,y,opt)
% bigperf_macroavg(X,y,opt);
% Computes the average classification accuracy per class.
%
% INPUTS:
% -X: input data bigarray
% -y: labels bigarray
% -OPT: structure of options with the following fields (and subfields):
%   fields that need to be set through previous gurls tasks:
%       -pred (set by the pred_* routines)
%   fields that need to be set through the bigdefopt function:
%       -nb_pred
% 
% OUTPUT: struct with the following fields:
% -acc: array of prediction accuracy for each class
% -forho: ""
% -forplot: ""


		if isfield(opt, 'perf')
			p = opt.perf;
		end

		nb_pred = opt.nb_pred;
		transpose = opt.pred.Transpose(true);
		T = y.Sizes();
		T = T{1};

		n_class = zeros(1,T);
		flatcost = zeros(1,T);

		for i1 = 1:y.BlockSize : y.NumItems
			i2 = min(i1 + y.BlockSize -1 , y.NumItems);
			y_block = y(i1 : i2, : );
			ypred_block = opt.pred(i1 : i2, : );

            [dummy,IY]= max(y_block,[],2);
            [dummy,IYpred]= sort(ypred_block,2,'descend');
        
            for i=1:size(y_block,1)
                flatcost(IY(i)) = flatcost(IY(i)) + ~ismember(IY(i),IYpred(i,1:nb_pred));
                n_class(IY(i)) = n_class(IY(i)) + 1;
            end
    
        end

        p.acc = 1-(flatcost./n_class);
        p.forho = 1-(flatcost./n_class); %paramsel needs a performance and not a cost...
        p.forplot = 1-(flatcost./n_class);
		opt.pred.Transpose(transpose);
