
classdef Potra3 < handle

  properties (SetAccess = private, Hidden = true)
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = Potra3()
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function delete( ~ )
      %% Destroy the C++ class instance
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out1 = F( self, varargin )
      % extract arguments
      x = varargin{1};
      y = varargin{2};
      z = varargin{3};
      % -----------------

      %F
      %    OUT1 = F(X,Y,Z)
      
      %    This function was generated by the Symbolic Math Toolbox version 9.0.
      %    09-Oct-2021 16:07:38
      
      out1 = [x.*3.0-cos(y.*z)-1.0./2.0;x.^2-y.^2.*6.25e+2;z.*2.0e+1+exp(-x.*y)+9.471975511965978];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out1 = J( self, varargin )
      % extract arguments
      x = varargin{1};
      y = varargin{2};
      z = varargin{3};
      % -----------------

      %J
      %    OUT1 = J(X,Y,Z)
      
      %    This function was generated by the Symbolic Math Toolbox version 9.0.
      %    09-Oct-2021 16:07:39
      
      t2 = x.*y;
      t3 = y.*z;
      t4 = sin(t3);
      t5 = -t2;
      t6 = exp(t5);
      out1 = reshape([3.0,x.*2.0,-t6.*y,t4.*z,y.*-1.25e+3,-t6.*x,t4.*y,0.0,2.0e+1],[3,3]);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [out1,out2] = FJ( self, varargin )
      % extract arguments
      x = varargin{1};
      y = varargin{2};
      z = varargin{3};
      % -----------------

      %FJ
      %    [OUT1,OUT2] = FJ(X,Y,Z)
      
      %    This function was generated by the Symbolic Math Toolbox version 9.0.
      %    09-Oct-2021 16:07:39
      
      t2 = x.*y;
      t3 = y.*z;
      t4 = sin(t3);
      t5 = -t2;
      t6 = exp(t5);
      out1 = [x.*3.0-cos(t3)-1.0./2.0;x.^2-y.^2.*6.25e+2;t6+z.*2.0e+1+9.471975511965978];
      if nargout > 1
          out2 = reshape([3.0,x.*2.0,-t6.*y,t4.*z,y.*-1.25e+3,-t6.*x,t4.*y,0.0,2.0e+1],[3,3]);
      end
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out1 = fun( self, varargin )
      % extract arguments
      x = varargin{1};
      y = varargin{2};
      z = varargin{3};
      % -----------------

      %FUN
      %    OUT1 = FUN(X,Y,Z)
      
      %    This function was generated by the Symbolic Math Toolbox version 9.0.
      %    09-Oct-2021 16:07:39
      
      out1 = (z.*2.0e+1+exp(-x.*y)+9.471975511965978).^2+(x.*-3.0+cos(y.*z)+1.0./2.0).^2+(x.^2-y.^2.*6.25e+2).^2;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out1 = grad( self, varargin )
      % extract arguments
      x = varargin{1};
      y = varargin{2};
      z = varargin{3};
      % -----------------

      %GRAD
      %    OUT1 = GRAD(X,Y,Z)
      
      %    This function was generated by the Symbolic Math Toolbox version 9.0.
      %    09-Oct-2021 16:07:39
      
      t2 = x.*y;
      t3 = y.*z;
      t4 = x.*3.0;
      t5 = x.^2;
      t6 = y.^2;
      t10 = z.*2.0e+1;
      t7 = cos(t3);
      t8 = sin(t3);
      t9 = -t4;
      t11 = -t2;
      t13 = t6.*6.25e+2;
      t12 = exp(t11);
      t14 = -t13;
      t15 = t7+t9+1.0./2.0;
      t16 = t5+t14;
      t17 = t10+t12+9.471975511965978;
      out1 = [t7.*-6.0+x.*1.8e+1+t16.*x.*4.0-t12.*t17.*y.*2.0-3.0,t16.*y.*-2.5e+3-t12.*t17.*x.*2.0-t8.*t15.*z.*2.0,t12.*4.0e+1+z.*8.0e+2-t8.*t15.*y.*2.0+3.788790204786391e+2];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function out1 = hess( self, varargin )
      % extract arguments
      x = varargin{1};
      y = varargin{2};
      z = varargin{3};
      % -----------------

      %HESS
      %    OUT1 = HESS(X,Y,Z)
      
      %    This function was generated by the Symbolic Math Toolbox version 9.0.
      %    09-Oct-2021 16:07:39
      
      t2 = x.*y;
      t3 = y.*z;
      t4 = x.*3.0;
      t5 = x.^2;
      t6 = y.^2;
      t7 = z.^2;
      t12 = z.*2.0e+1;
      t8 = t2.*2.0;
      t9 = cos(t3);
      t10 = sin(t3);
      t11 = -t4;
      t13 = -t2;
      t18 = t2.*5.0e+3;
      t14 = -t8;
      t15 = t10.^2;
      t16 = exp(t13);
      t19 = t10.*y.*6.0;
      t20 = t10.*z.*6.0;
      t21 = -t18;
      t28 = t9+t11+1.0./2.0;
      t17 = exp(t14);
      t23 = t3.*t15.*2.0;
      t24 = t16.*x.*4.0e+1;
      t25 = t16.*y.*4.0e+1;
      t29 = t10.*t28.*2.0;
      t32 = t3.*t9.*t28.*2.0;
      t34 = t12+t16+9.471975511965978;
      t22 = t8.*t17;
      t26 = -t24;
      t27 = -t25;
      t30 = -t29;
      t33 = -t32;
      t35 = t16.*t34.*2.0;
      t37 = t8.*t16.*t34;
      t31 = t19+t27;
      t36 = -t35;
      t38 = t23+t26+t30+t33;
      t39 = t20+t21+t22+t36+t37;
      out1 = reshape([t5.*1.2e+1-t6.*2.5e+3+t6.*t17.*2.0+t6.*t35+1.8e+1,t39,t31,t39,t5.*-2.5e+3+t6.*4.6875e+6+t5.*t17.*2.0+t7.*t15.*2.0+t5.*t35-t7.*t9.*t28.*2.0,t38,t31,t38,t6.*t15.*2.0-t6.*t9.*t28.*2.0+8.0e+2],[3,3]);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [out1,out2,out3] = fgH( self, varargin )
      % extract arguments
      x = varargin{1};
      y = varargin{2};
      z = varargin{3};
      % -----------------

      %fgH
      %    [OUT1,OUT2,OUT3] = fgH(X,Y,Z)
      
      %    This function was generated by the Symbolic Math Toolbox version 9.0.
      %    09-Oct-2021 16:07:39
      
      t2 = x.*y;
      t3 = y.*z;
      t4 = x.*3.0;
      t5 = x.^2;
      t6 = y.^2;
      t7 = z.^2;
      t12 = z.*2.0e+1;
      t8 = t2.*2.0;
      t9 = cos(t3);
      t10 = sin(t3);
      t11 = -t4;
      t13 = -t2;
      t18 = t2.*5.0e+3;
      t21 = t6.*6.25e+2;
      t14 = -t8;
      t15 = t10.^2;
      t16 = exp(t13);
      t19 = t10.*y.*6.0;
      t20 = t10.*z.*6.0;
      t22 = -t18;
      t23 = -t21;
      t30 = t9+t11+1.0./2.0;
      t17 = exp(t14);
      t25 = t3.*t15.*2.0;
      t26 = t16.*x.*4.0e+1;
      t27 = t16.*y.*4.0e+1;
      t31 = t5+t23;
      t32 = t10.*t30.*2.0;
      t35 = t3.*t9.*t30.*2.0;
      t37 = t12+t16+9.471975511965978;
      out1 = t30.^2+t31.^2+t37.^2;
      if nargout > 1
          out2 = [t9.*-6.0+x.*1.8e+1+t31.*x.*4.0-t16.*t37.*y.*2.0-3.0,t31.*y.*-2.5e+3-t16.*t37.*x.*2.0-t10.*t30.*z.*2.0,t16.*4.0e+1+z.*8.0e+2-t10.*t30.*y.*2.0+3.788790204786391e+2];
      end
      if nargout > 2
          t24 = t8.*t17;
          t28 = -t26;
          t29 = -t27;
          t33 = -t32;
          t36 = -t35;
          t38 = t16.*t37.*2.0;
          t40 = t8.*t16.*t37;
          t34 = t19+t29;
          t39 = -t38;
          t41 = t25+t28+t33+t36;
          t42 = t20+t22+t24+t39+t40;
          out3 = reshape([t5.*1.2e+1-t6.*2.5e+3+t6.*t17.*2.0+t6.*t38+1.8e+1,t42,t34,t42,t5.*-2.5e+3+t6.*4.6875e+6+t5.*t17.*2.0+t7.*t15.*2.0+t5.*t38-t7.*t9.*t30.*2.0,t41,t34,t41,t6.*t15.*2.0-t6.*t9.*t30.*2.0+8.0e+2],[3,3]);
      end
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
end
