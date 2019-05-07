function [currentState, logP] = hmmviterbi(self,data,opt_varargin)


    n=0;            
    if isempty(opt_varargin)
        program_option;
        n=size(opt_varargin,2);
    elseif isstruct(opt_varargin{1})
        opt=opt_varargin{1};
    else
        program_option;
        n=size(opt_varargin,2);
    end
    ma=[];
    for j=1:2:n
        opt=setfield(opt,opt_varargin{j},opt_varargin{j+1});
    end

    nstates=self.nstates;
    ndata=self.emis_model.ndatafun(data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute emission probabilities of sequence given each state: e(state,1:L)
    emit = zeros(nstates,ndata);
    for state = 1:nstates
        [aux logE(state,:)] = self.emis_model.prob(opt.train,data,state);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A2 logTR]  = self.trans_model.expect(opt.train);
    % allocate space
    pTR = zeros(nstates,ndata);
    % assumption is that model is in state 1 at step 0
    v = -Inf(nstates,1);
    v(1,1) = 0;
    vOld = v;
    % loop through the model
    for count = 1:ndata
        for state = 1:nstates
            % for each state we calculate
            % v(state) = e(state,seq(count))* max_k(vOld(:)*tr(k,state));
            bestVal = -inf;
            bestPTR = 0;
            % use a loop to avoid lots of calls to max
            for inner = 1:nstates 
                val = vOld(inner) + logTR(inner,state);
                if val > bestVal
                    bestVal = val;
                    bestPTR = inner;
                end
            end
            % save the best transition information for later backtracking
            pTR(state,count) = bestPTR;
            % update v
            v(state) = logE(state,count) + bestVal;
        end
        vOld = v;
    end
    % decide which of the final states is post probable
    [logP, finalState] = max(v);

    % Now back trace through the model
    currentState(ndata) = finalState;
    for count = ndata-1:-1:1
        currentState(count) = pTR(currentState(count+1),count+1);
        if currentState(count) == 0
            error(message('stats:hmmviterbi:ZeroTransitionProbability', currentState( count + 1 )));
        end
    end
end