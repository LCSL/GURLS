function [scores] = pred_bagging(X, y, opt)

t_opt = struct;
t_opt.rls = struct;


switch opt.predbagmethod
	case 'avg'
	%% 1. averaged score
		avg_scores = 0;
		for i = 1:numel(opt.rls.W) % for each of the bagging classifier
			t_opt.rls.W = opt.rls.W{i};
			avg_scores = avg_scores + pred_primal(X, y, t_opt); % nxT (n: # of data, T: # of class)
		end
		scores = avg_scores;
	case 'vote'

		%% 2. voting - MAR 25, 2011. - sangwan
		vote_scores=0;
		zero_vote=zeros(size(X,1),size(y,2)); % nxT
		for i = 1:numel(opt.rls.W) % for each of the bagging classifier
			t_opt.rls.W = opt.rls.W{i}; % i: i-th classifier
		    pred_vector=pred_primal(X, y, t_opt);
		    [tmp vote_ind]=max(pred_vector');
		    vote_mat=zero_vote;
		    for j=1:1:length(vote_ind)
		        vote_mat(j,vote_ind(j))=1;
		    end
			vote_scores = vote_scores + vote_mat; % nxT (n: # of data, T: # of class)
		end
		scores = vote_scores;
end 
