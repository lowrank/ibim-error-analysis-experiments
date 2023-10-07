function U = RandUGroup(n, varargin)
% Q = RandUGroup(n, sz1, sz2, ...)
% Q = RandUGroup(n, sz)
% Q = RandUGroup(n, sz, 'GType', GType)
% U = RandUGroup(n, sz, 'GType', GType)
% ... = RandUGroup(n, sz, GType)
% ... = RandUGroup(n, GType) % or 
% ... = RandUGroup(n, 'GType', GType) % assume default sz = 1
% U = RandUGroup(..., Propertie1, value1, ...)
%
% Generate matrix of one of these four supported types of groups:
%    O(n), SO(n), U(n), SU(n)
% Inputs:
%   n:  space dimension, must be provided
%   sz: (optional argument) size of of (n x n) matrices,
%       see Output decription. If sz is not provided RandUGroup assumes 
%       sz is 1 (scalar output)
%   sz can be split as comma list in input arguments: sz1, sz2, ...
%   GType is among
%       'O':  orthogonal group (real), Q.'*Q = Q*Q.'= eye(n)
%       'SO': special orthogonal group (real), Q.'*Q = Q*Q.'= eye(n) and det(Q) = 1
%       'U':  unitary group (complex), U'*U = U*U'= eye(n)
%       'SU': special unitary group (complex), U'*U = U*U'= eye(n) and det(U) = 1
%   Optional Properties (put them after GType argument) are:
%       'DataType': 'double' (default) or 'single'
%       'UseQRFlag': scalar boolean orthogonalization engine
%           if TRUE (default) uses qr, otherwise eig.
%           It can also be among {'MultipleQR', 'qr'} to select qr engine
%           The MultipleQR engine requires the File Exchange
%           https://www.mathworks.com/matlabcentral/fileexchange/68976-multipleqr
%       'RandFun': among {'randn','randsym'}, default 'randn'
%           randomize engine before orthogonalization
%       'DetFun': among {'LUdet','Laplacedet'}. Applied only for special
%           group to compute the determinant. 'LUdet' uses MATLAB det()
%           'Laplacedet' uses recursive Laplacian formula.
%           Default 'LUdet' for n>=4 and 'Laplacedet' for n<=3.
% Output:
%   Q/U: array of size (n x n x sz), each U(:,:,k) satisfies the above
%       description
%
% EXAMPLES:
%
% >> Q=RandUGroup(3,'O')
% 
% Q =
% 
%    -0.6611    0.4898    0.5684
%    -0.2935    0.5284   -0.7966
%     0.6906    0.6934    0.2056
% 
% Q'*Q
% 
% ans =
% 
%     1.0000    0.0000    0.0000
%     0.0000    1.0000   -0.0000
%     0.0000   -0.0000    1.0000
% 
% det(Q)
% 
% ans =
% 
%    -1.0000
% 
% Q=RandUGroup(3,'SO');
% det(Q)
% 
% ans =
% 
%     1.0000
% 
% U=RandUGroup(3,'U');
% det(U)
% 
% ans =
% 
%    0.9172 - 0.3985i
% 
% U=RandUGroup(3,'SU');
% det(U)
% 
% ans =
% 
%    1.0000 + 0.0000i
% 
% U=RandUGroup(3,1,2,'O')
% 
% U(:,:,1,1) =
% 
%     0.2770    0.7586    0.5898
%     0.9447   -0.3270   -0.0230
%    -0.1754   -0.5635    0.8072
% 
% 
% U(:,:,1,2) =
% 
%     0.4886   -0.6677    0.5616
%     0.8617    0.2681   -0.4309
%    -0.1371   -0.6945   -0.7063
%
% Author: brunoluong@yahoo.com
% History: 23/Hov/2020 initial version
argin = varargin;
% parse sz
sz = 1;
if ~isempty(argin)
    ilast = find(cellfun(@isnumeric, argin),1,'last');
    if ~isempty(ilast)
        sz = cat(2,argin{1:ilast});
        argin = argin(ilast+1:end);
    end
end
if ~isempty(argin)
    arg1 = argin{1};
    % calling syntax RandUGroup(n, p, GType)
    if any(strcmpi(arg1, {'O' 'SO' 'U' 'SU'}))
        argin = [{'GType'} argin];
    end
end
options = struct(argin{:});
if nargin < 1 || isempty(n)
    n = 2;
end
if isscalar(sz)
    sz = [sz sz];
end
p = prod(sz);
GType = getoptions(options, 'GType', 'O');
GType = strtrim(GType);
if contains(GType, {'u' 'U'})
    DefComplexFlag = true;
else
    DefComplexFlag = false;
end
if n <= 3
    DefDetfun = 'Laplacedet';
else
    DefDetfun = 'LUdet';
end
DefSpecial      = upper(GType(1)) == 'S' ;
ComplexFlag     = getoptions(options, 'ComplexFlag',    DefComplexFlag);	% complex or real
UseQRFlag       = getoptions(options, 'UseQRFlag',      true);              % QR or EIG
RandFun         = getoptions(options, 'RandFun',        'randn');           % among {'randn','randsym','rand'}
DetFun          = getoptions(options, 'DetFun',         DefDetfun);         % among {'LUdet','Laplacedet'}
Special         = getoptions(options, 'Special',        DefSpecial);
DataType        = getoptions(options, 'DataType',       'double');
switch RandFun
    case 'randn'
        rfun = @(varargin) randn(varargin{:}, DataType);
    case 'randsym'
        rfun = @(varargin) 2*rand(varargin{:}, DataType)-1;  % only for EIG
        if UseQRFlag
            warning('randsym should no used with QR');
        end
    case 'rand'
        rfun = @(varargin) rand(varargin{:}, DataType);      % NOT this for any method
        if UseQRFlag
            warning('rand should no be used');
        end
end
if ComplexFlag
    U = rfun(n,n,p) + 1i*rfun(n,n,p);
else
    U = rfun(n,n,p);
end
detU = NaN;
if UseQRFlag
    if ischar(UseQRFlag)
        MQRFlag = strcmpi(UseQRFlag, 'MultipleQR');
    else
        MQRFlag = n <= 4;
    end
    MQRFlag = MQRFlag && exist('MultipleQR','file') == 2;
    if MQRFlag
        % FEX https://www.mathworks.com/matlabcentral/fileexchange/68976-multipleqr
        U = MultipleQR(U);
        % Beware: We assume MultipleQR that uses Householder reflection,
        % thus returns det that depends only on dimension as below
        detU = (-1)^n;
    else
        for k=1:p
            [U(:,:,k),~] = qr(U(:,:,k));
        end
        if ~ComplexFlag
            % Beware: We assume MATLAB qr that uses Householder reflection,
            % thus returns det that depends only on dimension as below
            % this is the case up to R2020b
            detU = (-1)^(n-1);
        end
    end
else
    U = U + permute(conj(U),[2 1 3]);
    for k=1:p
        [U(:,:,k),~] = eig(U(:,:,k));
    end
    % NOTE:
    % for complex, eig of hermitian matrix selects phase st U(end,:,:) are real
    %                  of non hermitian selects phase st U(imax,:,:) are
    %                  real >= 0, imax = argmax_i(abs(U(i,:,:))
end
% random phase
if ComplexFlag
    rsign = exp((2i*pi)*rand(1,n,p));
else
    rsign = 2*(rand(1,n,p)>0.5)-1;
end
if Special
    % force determinant to be 1 by multipling the first vector with a
    % scalar (of norm 1)
    if isfinite(detU)
        % known determinant from the orthogonlization algorithm
        d = detU;
    else
        switch DetFun
            case 'Laplacedet'
                d = Laplaciandet(U);
            otherwise % LUdet
                d = zeros(1,1,p);
                for k=1:p
                    d(k) = det(U(:,:,k));
                end
        end
        d = d./abs(d); % in case numerical error
    end
    rsign(1,1,:) = 1;
    rsign(1,1,:) = 1./(d.*prod(rsign,2));
end
    
U = U.*rsign;
% Put result to final shape
U = reshape(U, [n n sz]);
end % RandUGroup
% Compute determinant by Laplacian recursive formula
function d = Laplaciandet(A)
n = size(A,1);
rdet(n);
d = rdet(A,1:n,1:n);
end
% Compute determinant recursively with bookeeping
function d = rdet(A,r,c)
persistent RC D
if nargin == 1
    % Reset data-base
    q = A;
    RC = cell(1,q);
    for k=1:q
        RC{k} = zeros(0,2*k);
    end
    D = cell(1,q);
    return
end
q = length(r);
if q == 1
    d = A(r,c,:);
    return
end
[tf,loc] = ismember([r,c],RC{q},'rows');
if tf
    d = D{q}(loc,1,:);
    return
end
% Update the row-column database
RC{q}(end+1,:) = [r,c];
% Recursive call of determinant
A1 = A(r,c(1),:);
c(1) = [];
di = zeros(size(A1));
for i=1:q
    ri = r;
    ri(i) = [];
    di(i,1,:) = rdet(A,ri,c);
end
A1(2:2:end,1,:) = -A1(2:2:end,1,:);
d = sum(A1.*di,1);
% Store the result
D{q}(end+1,1,:) = d;
end
%%
function value = getoptions(options, name, defaultvalue)
% function value = getoptions(options, name, defaultvalue)
fields = fieldnames(options);
found = strcmpi(name,fields);
if any(found)
    value = options.(fields{found});
    if isempty(value)
        value = defaultvalue;
    end
else
    value = defaultvalue;
end
end % getoptions