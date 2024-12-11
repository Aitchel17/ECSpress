function [stack] = pre_pshiftcorrection(stack,pshift) 
   
    if pshift == 0
        disp('skip pshift correction, pshift =0')
    elseif pshift > 0
        stack(2:2:end, 1+abs(pshift):end, :) = stack(2:2:end, 1:end-abs(pshift), :);
            stack = stack (:,1+abs(pshift):end,:);
    elseif pshift < 0
        stack(1:2:end, 1+abs(pshift):end, :) = stack(1:2:end, 1:end-abs(pshift), :);
        stack = stack (:,1+abs(pshift):end,:);
    end
    
end

