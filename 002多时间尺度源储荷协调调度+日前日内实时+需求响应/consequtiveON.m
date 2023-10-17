function C = consequtiveON(x,minup)
if min(size(x))==1
    x = x(:)';
end
if size(x,1) ~= size(minup,1)
    error('MINUP should have as many rows as X');
end
Horizon = size(x,2);
C = [];
for k = 2:size(x,2)
    for unit = 1:size(x,1)
        % indicator will be 1 only when switched on
        indicator = x(unit,k)-x(unit,k-1);
        range = k:min(Horizon,k+minup(unit)-1);
        % Constraints will be redundant unless indicator = 1
        affected = x(unit,range);
        if strcmp(class(affected),'sdpvar')
        % ISA behaves buggy, hence we use class+strcmp
            C = [C, affected >= indicator];
        end
    end
end