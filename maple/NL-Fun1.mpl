Fun1 := module()

  description "fun1";

  # Module defined as a package (i.e.) collection of procedures
  option package, load = ModuleLoad, unload = ModuleUnLoad;

  export fun, grad, hess, F, JF;
  
  local ModuleLoad,
        ModuleUnLoad,
        f1, f2, f12,
        CST;

  uses LinearAlgebra;

  ModuleLoad := proc()
    CST := 53/27;
    NULL;
  end proc;

  ModuleUnLoad := proc()
    NULL;
  end proc:
  
  # Explicitly call ModuleLoad here so the type is registered when this
  # code is cut&pasted in.  ModuleLoad gets called when the module is
  # read from the repository, but the code used to define a module
  # (like the command below) is not called at library read time.
  ModuleLoad();

  f1  := unapply( sin(3*x+5*y)+x+exp(3*(x-y))-1, (x,y) );
  f2  := unapply( (CST*x^2+2*y^2+cos(x-y)^2)-1, (x,y) );
  f12 := unapply( f1(x,y)^2+f2(x,y)^2, (x,y) );

  fun := proc( X::Vector )
    local x, y;
    x := X[1];
    y := X[2];
    f12(x,y);
  end proc;

  grad := proc( X::Vector )
    local x, y;
    x := X[1];
    y := X[2];
    eval(<D[1](f12)(x,y),D[2](f12)(x,y)>);
  end proc;

  hess := proc( X::Vector )
    local x, y, d11, d12, d22;
    x   := X[1];
    y   := X[2];
    d11 := eval(D[1,1](f12)(x,y));
    d12 := eval(D[1,2](f12)(x,y));
    d22 := eval(D[2,2](f12)(x,y));
    <<d11,d12>|<d12,d22>>;
  end proc;

  F := proc( X::Vector )
    local x, y;
    x := X[1];
    y := X[2];
    <f1(x,y),f2(x,y)>;
  end proc;

  JF := proc( X::Vector )
    local x, y;
    x := X[1];
    y := X[2];
    <<D[1](f1)(x,y),D[1](f2)(x,y)>|
     <D[2](f1)(x,y),D[2](f2)(x,y)>>;
  end proc;

end module:
