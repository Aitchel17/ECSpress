function p = util_pbar(data, h)
persistent TOTAL COUNT H
if nargin == 2
    % initialisation mode
    H = h;
    TOTAL = data;
    COUNT = 0;
    
else
    % afterEach call, increment COUNT
    COUNT = 1 + COUNT;
    p = COUNT / TOTAL;
    waitbar(p, H, sprintf('Calculating correlation... %d %% (%d/%d)',round(p*100),COUNT,TOTAL));
end
end