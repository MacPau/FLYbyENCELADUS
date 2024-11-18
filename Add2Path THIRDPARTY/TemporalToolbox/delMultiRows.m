function M = delMultiRows(M)


%% DELMULTIROWS
% Function that deletes identical rows.
%
%% INPUTS
%   M       2D-Matrix
%
%% OUTPUTS
%   M       2D-Matrix M without identical rows
%
%% EXEMPLE
%   M = [1 2;
%        3 4;
%        5 6;
%        1 2;
%        2 1;
%        5 6;
%        5 6;
%        7 8;];
% 
%   M = delMultiRows(M)
%       returns M = 
%                      1     2
%                      3     4
%                      5     6
%                      2     1
%                      7     8
%
%
% (c)  Nicolas Croisard (2007)


k=1;
while k < size(M,1)
    
    ind = find(ismember(M,M(k,:),'rows'));
    M(ind(2:end),:) = [];
    k = k+1;

end