% Find conditional number of a sparse matrix
%
% matlab -nodesktop
% 1. copy matrixmarket -> ascii datafile
% 2. remove comments then 
% 3. change first line of `.mtx` to `_ _ 0`

matrs = [                   %        k
    'data/s_torsoL.dat',    %    2.5803e+67
    'data/s_tsopfL.dat'     %    2.5803e+67
]

for matr = matrs
    % http://matlab.izmiran.ru/help/techdoc/ref/load.html
    L = load('-ascii', matr);
    % https://www.mathworks.com/help/matlab/ref/spconvert.html
    S = spconvert(L);
    % https://www.mathworks.com/help/symbolic/cond.html
    k = cond(S)
end