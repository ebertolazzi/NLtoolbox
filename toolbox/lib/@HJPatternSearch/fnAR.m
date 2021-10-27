
% =================================================================
% NRMG Algorithm
% =================================================================

function R = fnAR(~,X)
  %FNAR This method implements the N-dimensional Rotation Matrix
  %Generation Algorithm (NRMG), here used to align the first direction
  %of the search direction matrix to the descending gradient of
  %the value function. Reference to NRMG can be found in DOI:
  %10.5923/j.ajcam.20170702.04 .
    
  N    = length(X); %X have to be row vector (transposed)
  R    = eye(N);    %Initial rotation matrix = Identity matrix
  step = 1;         %Initial step
  % Loop to create matrices of stages
  while step < N
    A = eye(N);
    n = 1;
    while n <= N-step
      r2 = X(n)*X(n) + X(n+step)*X(n+step);
      if r2 > 0
        r = sqrt(r2);
        pcos = X(n)/r;
        psin = -X(n+step)/r;
        % Base 2-dimensional rotation
        A(n, n) = pcos;
        A(n, n+step) = -psin;
        A(n+step, n) = psin;
        A(n+step, n+step) =  pcos;
      end
      n = n+2*step;  % Move to the next base operation
    end
    step = step*2;
    X    = (A*X')';
    R    = A*R;  % Multiply R by current matrix of stage A
  end
end
