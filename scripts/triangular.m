% compute solution of lower triangular solve
%
% 1. copy matrixmarket -> ascii datafile
% 2. remove comments and remove line for `(m n nnz)`
# 3. for `b` need to add additional line to give matrix size `ncols 1 0`

data_folder = './data/'

matrs = [
    's_small',
    's_medium',
    's_torso',
    's_tsopf'
]

for matr = ["s_tsopf"]
    % reads matrix market
    L = load('-ascii', strcat(data_folder, matr, 'L.dat'));
    b = load('-ascii', strcat(data_folder, matr, 'b.dat'));
    LL = spconvert(L);
    bb = spconvert(b);
    % solve
    xx = LL \ bb;
    [i,j,val] = find(xx);
    % outputs to matrix market
    [x_m ~] = size(LL);
    x_n = 1;
    [x_nnz ~] = size(i);
    fid = fopen(strcat(data_folder, matr, 'x.mtx'), 'w');
    fprintf(fid, '%s\n', "%%MatrixMarket matrix coordinate real general");      % banner
    fprintf(fid, '%d  %d  %20.d\n', x_m, x_n, x_nnz);                           % m n nnz
    fprintf(fid, '%d  %d  %20.16f\n', [i,j,val]');                              % row col val
    fclose(fid);
end


