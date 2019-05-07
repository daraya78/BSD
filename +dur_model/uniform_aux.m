classdef uniform_aux < handle
    properties
        parsampl
        nstates
    end
    methods
        function self=uniform_aux(nstates,mean,vari)
            if exist('mean','var'), self.parsampl.mean=mean; else self.parsampl.mean=[]; end
            if exist('vari','var'), self.parsampl.var=vari; else self.parsampl.var=[]; end
            if exist('nstates','var'), self.nstates=nstates; else self.nstates=[]; end
        end
        function self=parsamplfun(self,option)
     
        end
        
        function re = sample(self,num,state)
            ma=self.parsampl.mean{state}+self.parsampl.var{state};
            mi=self.parsampl.mean{state}-self.parsampl.var{state};
            re=randi([mi ma],[num 1]);
            re = max(1,re);
        end
        
        
        
    end
end

