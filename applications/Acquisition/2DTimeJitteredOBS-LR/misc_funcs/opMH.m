classdef opMH < opSpot
    % source-receiver to midpoint-offset conversion.
    %
    % use:
    %   op = opMH(s,r)
    %
    %
    % You may use this code only under the conditions and terms of the
    % license contained in the file LICENSE provided with this source
    % code. If you do not agree to these terms you may not use this
    % software.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        s,r;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opMH(s,r)
            
            op = op@opSpot('opMH',(min(s,r)+floor(abs(r-s)/2))*(s+r-1),s*r);
            op.cflag     = 0;
            op.linear    = 1;
            op.children  = [];
            op.sweepflag = true;
            op.s         = s;
            op.r         = r;
        end
    end
    
    
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = multiply(op,x,mode)
            if mode==1
                out = SR2MH(x,op.s,op.r,1);
            else
                out = SR2MH(x,op.s,op.r,-1);
            end
        end %multiply
    end %protected methods
    
end %classdef


function [B] =SR2MH(A1,s,r,mode)

if mode==1
    A1 = reshape(A1,s,r);
    B=zeros(min(s,r)+floor(abs(r-s)/2),s+r-1);
    maxx=r-1;
    for k=1-s:1:maxx
        B(floor(abs(k)/2)+1:floor(abs(k)/2)+length(diag(A1,k)),s+k)=diag(A1,k).';
    end
    B = vec(B);
else
    
    s1 = min(s,r)+floor(abs(r-s)/2);
    r1 = s+r-1;
    A1 = reshape(A1,s1,r1);
    B=zeros(s,r);
    Blow=zeros(s+r-1,s+r-1);
    
    if (s>r)
        s1=r;
        r1=s;
        s=s1;
        r=r1;
        A1=fliplr(A1);
        B=B.';
        blow_size=s+r-1;
        Blow=zeros(blow_size,blow_size);
        maxx=r-1;
        max_trans=abs(s-r);
        min_rhs=max_trans+1;
        
        for k=1-s:1:0   %abs(s-r)
            Blow(abs(k)+1:2:2*min(r,s)-1+k,s+k)=A1(floor(abs(k)/2)+1:s+round(k/2),s+k);
        end
        
        
        for k=1:1:max_trans
            Blow(abs(k)+1:2:2*min(r,s)-1+k,s+k)=A1(floor(k/2)+1:floor(k/2)+s,s+k);
        end
        
        for k=min_rhs:1:maxx
            Blow(k+1:2:blow_size-abs(k-max_trans),k+s)=A1(1+floor(k/2):floor(k/2)+s-abs(k-max_trans),k+s);
        end
        
    elseif ((s<=r))
        s1=r;
        r1=s;
        
        blow_size=s+r-1;
        Blow=zeros(blow_size,blow_size);
        maxx=r-1;
        max_trans=abs(s-r);
        min_rhs=max_trans+1;
        
        for k=1-s:1:0   %abs(s-r)
            Blow(abs(k)+1:2:2*min(r,s)-1+k,s+k)=A1(floor(abs(k)/2)+1:s+round(k/2),s+k);
        end
        
        
        for k=1:1:max_trans
            Blow(abs(k)+1:2:2*min(r,s)-1+k,s+k)=A1(floor(k/2)+1:floor(k/2)+s,s+k);
        end
        
        for k=min_rhs:1:maxx
            Blow(k+1:2:blow_size-abs(k-max_trans),k+s)=A1(1+floor(k/2):floor(k/2)+s-abs(k-max_trans),k+s);
        end
        
    end
    
    for k = 1:s
        tmp= diag(Blow,s-1+2*(1-k));
        B(k,:) = tmp((s-1)/2-abs(k-1-(s-1)/2) + [1:r]);
    end
    
    
    if r1>s1
        B=transpose(B);
    end
    B = vec(B);
end %function
end


