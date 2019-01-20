% Copyright (c) 2013, Guido Albertin All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [YI] = lsq_lut_piecewise(x, y, XI)
% LSQ_LUT_PIECEWISE Piecewise linear interpolation for 1-D interpolation
% (table lookup)
%   YI = lsq_lut_piecewise( x, y, XI ) obtain optimal (least-square sense)
%   vector to be used with linear interpolation routine. The target is
%   finding Y given X the minimization of function
%           f = |y-interp1(XI,YI,x)|^2
%   
%   INPUT
%       x measured data vector
%       y measured data vector
%       XI break points of 1-D table
%
%   OUTPUT
%       YI interpolation points of 1-D table
%           y = interp1(XI,YI,x)
%


if size(x,2) ~= 1
    error('Vector x must have dimension n x 1.');   
elseif size(y,2) ~= 1
    error('Vector y must have dimension n x 1.');    
elseif size(x,1) ~= size(x,1)
    error('Vector x and y must have dimension n x 1.'); 
end

% matrix defined by x measurements
A = sparse([]); 

% vector for y measurements
Y = []; 

for j=2:length(XI)
    
    % get index of points in bin [XI(j-1) XI(j)]
    ix = x>=XI(j-1) & x<XI(j);
    
    % check if we have data points in bin
    if ~any(ix)
        warning(sprintf('Bin [%f %f] has no data points, check estimation. Please re-define X vector accordingly.',XI(j-1),XI(j)));
    end
    
    % get x and y data subset
    x_ = x(ix);
    y_ = y(ix);
    
    % create temporary matrix to be added to A
    tmp = [(( -x_+XI(j-1) ) / ( XI(j)-XI(j-1) ) + 1) (( x_-XI(j-1) ) / ( XI(j)-XI(j-1) ))];
    
    % build matrix of measurement with constraints
    [m1,n1]=size(A);
    [m2,n2]=size(tmp);
    A = [[A zeros(m1,n2-1)];[zeros(m2,n1-1) tmp]];
    
    % concatenate y measurements of bin
    Y = [Y; y_];
end

% obtain least-squares Y estimation
YI=A\Y;

