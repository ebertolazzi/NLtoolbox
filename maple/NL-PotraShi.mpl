PotraShi := module()

  description "linesearch"
              "Efficient Line Search Algorithm for Unconstrained Optimization"
              "F. A. POTRA AND Y.SHI"
              "JOURNAL OF OPTIMIZATION THEORY AND APPLICATIONS"
              "Vol.85,No.3,pp.677-704,1995";

  # Module defined as a package (i.e.) collection of procedures
  option package, load = ModuleLoad, unload = ModuleUnLoad;

  export lines, linesearch_A2, esatta, AB_list, doplot, doplot2;
  
  local ModuleLoad,
        ModuleUnLoad,
        _debug,
        min_quadratic,
        min_cubic,
        sigma, rho, tau1, tau2, tau3, Jbig,
        max_iter, A_list, B_list;

  uses LinearAlgebra;

  ModuleUnLoad := proc()
    printf("PotraShi unload\n");
    NULL;
  end proc:

  ModuleLoad := proc()
    printf("PotraShi load\n");
    _debug   := true;
    rho      := 0.25;
    sigma    := 0.75;
    tau1     := 0.1;
    tau2     := 0.49;
    tau3     := 2.5;
    Jbig     := 2;
    max_iter := 40;
    NULL;
  end proc;
  
  # Explicitly call ModuleLoad here so the type is registered when this
  # code is cut&pasted in.  ModuleLoad gets called when the module is
  # read from the repository, but the code used to define a module
  # (like the command below) is not called at library read time.
  ModuleLoad();

  AB_list := proc()
    A_list, B_list;
  end proc;

  min_quadratic := proc ( a, fa, Dfa, b, fb )
    local t, delta; # A+B*(x-a)+C*(x-a)^2
    # p  = fa + Dfa*(x-a)+t*(x-a)^2
    # p' = Dfa+2*t*(x-a)
    ##t     := ((fb-fa)/(b-a)-Dfa)/(b-a);
    ##delta := -Dfa/(2*t);
    t     := (fb-fa)/(b-a)-Dfa;
    delta := max(tau1,min(1-tau1,-Dfa/(2*t)));
    a+(b-a)*delta;
  end proc;

  lines := proc( f0, df0, L )
    [ [0,f0], [L,f0+L*rho*df0] ], [ [0,f0], [L,f0+L*sigma*df0] ];
  end proc;

  linesearch_A2 := proc( f_in, df_in )
    local f, df, a, b, c, alpha,
          fa, dfa, fb, fc,
          r_df0, s_df0,
          f0, df0, df1, t1, t2, delta, k;

    f  := x->evalf(f_in(x));
    df := x->evalf(df_in(x));

    alpha  := 1;
    a      := 0;
    b      := 1;
    fa     := f(0);
    fb     := f(1);
    df0    := df(0);

    if df0 >= 0 then
      error "[linesearch_A2] df0 = %1 must be < 0!\n", df0;
    end;

    r_df0  := rho*df0;
    s_df0  := sigma*df0;
    f0     := fa;

    A_list := [];
    B_list := [];

    # step 1 (check alpha = 1)
    if fb <= f0 + r_df0 then
      if sigma > 0.5 then
        if fb >= f0 + s_df0 then
          if _debug then
            printf("[step1] return 1 (a case)\n");
          end;
          return 1;
        end;
      else
        df1 := df(1);
        if df1 >= s_df0 then
          if _debug then
            printf("[step1] return 1 (b case)\n");
          end;
          return 1;
        end;
      end;
      # step 2
      if _debug then
        if sigma > 0.5 then
          printf(
            "\n[step2] now\n"
            "f(1) [%g] <= f(0)+rho*f'(0) [%g]\n"
            "and f(1) [%g] < f(0)+sigma*f'(0) [%g]\n",
            fb, f0 + r_df0,
            fb, f0 + s_df0
          );
        else
          printf(
            "\n[step2] now\n"
            "f(1) [%g] <= f(0)+rho*f'(0) [%g]\n"
            "and f'(1) [%g] < sigma*f'(0) [%g]\n",
            fb, f0 + r_df0,
            df1, s_df0
          );
        end;
      end;
      a      := 1;
      b      := Jbig;
      A_list := [op(A_list), a];
      B_list := [op(B_list), b];
      fa := fb;
      fb := f(b);
      while fb <= fa + (b-a)*r_df0 do

        if fb >= fa + (b-a)*s_df0 then
          if _debug then
            printf( "[step2] return b = %g\n", b );
          end;
          return b;
        end;

        a  := b;
        b  := Jbig*b;
        fa := fb;
        fb := f(b);

        if _debug then
          printf( "[step2] a=%g b=%g f(a)=%g f(b)=%g\n", a, b, fa, fb );
        end;

        A_list := [op(A_list), a];
        B_list := [op(B_list), b];

        if a > 1000 then
          error "a too big\n";
        end;
      end;
    else
      printf(
        "\n[step2] case\n"
        "f(1) [%g] > f(0)+rho*f'(0) [%g]\n",
        fb, f0 + r_df0
      );
    end;
    
    # aggiunta ad algoritmo
    if a > 0 then
      dfa := df(a);
      if dfa >= 0 then
        if _debug then
          printf("[step3] return a = %g (added)\n", a);
        end;
        return a;
      end;
    else
      dfa := df0;
    end;

    # step 3
    if _debug then
      printf(
        "\n[step3] a=%g b=%g with f(0)=%g, f(a)=%g, f'(a)=%g\n"
        "f(a) [%g] <= f(0)+a*rho*f'(0) [%g]\n"
        "and f(b) [%g] > f(a)+(b-a)*sigma*f'(0) [%g]\n",
        a, b, f0, fa, dfa,
        fa, f0 + a*r_df0,
        fb, fa + (b-a)*s_df0
      );
    end;

    for k from 1 to max_iter do

      c := min_quadratic( a, fa, dfa, b, fb );
      #c  := (a+b)/2;
      fc := f(c);

      if fc <= fa + (b-a)*r_df0 and fc >= fa + (b-a)*s_df0 then
        if _debug then
          printf("[step3] return c = %g (case a)\n", c);
        end;
        return c;
      end;

      t1    := (fb-fc)/(b-c);
      t2    := (fc-fa)/(c-a);
      delta := abs( (t1-t2)/(b-a) );

      if fc <= fa + (c-a)*r_df0 then
        if (rho-sigma)*df0 >= tau3*(b-a)*delta then
          if _debug then
            printf("[step3] return c = %g (case b)\n", c);
          end;
          return c;
        end;
        a   := c;
        fa  := fc;
        dfa := df(c);
      else
        if (rho-sigma)*df0 >= tau3*(b-a)*delta and a > 0 then
          if _debug then
            printf("[step3] return a = %g\n", a);
          end;
          return a;
        end;
        b  := c;
        fb := fc;
      end;

      A_list := [op(A_list), a];
      B_list := [op(B_list), b];

      if _debug then
        printf(
          "[step3] a=%.4g b=%.4g f(a)=%.4g f'(a)=%.4g f(b)=%.4g f'(b)=%.4g\n",
          a, b, fa, dfa, fb, df(b)
        );
      end;

      if a > 1000 then
        error "a = %1\n", a;
      end;
    end;

    printf("[step3] search failed\n");
    return NULL;
    
  end proc;

  esatta := proc( fun )
    local res, alpha;
    #res := Optimization[Minimize](
    res := Optimization[NLPSolve](
      fun(alpha),
      {alpha >= 0, alpha <= 3},
      #initialpoint={alpha=1},
      iterationlimit=50
    ):
    subs(res[2],alpha);
  end proc:

  doplot := proc( f_in, df_in, alpha, amax)
    local f0, df0, x, AA, BB, CC, DD;
    f0  := evalf(f_in(0));
    df0 := evalf(df_in(0));

    AA := plot(f_in(x),x=0..amax);
    BB := plot([[0,f0],[amax,f0+amax*rho*df0]],color="LimeGreen");
    CC := plot([[0,f0],[amax,f0+amax*sigma*df0]],color="blue");
    DD := plot([[alpha,f_in(alpha)]],color="red",style = point);
    display(AA,BB,CC,DD);
  end;

  doplot2 := proc( f_in, df_in, alpha, amax, A_list, B_list )
    local i, f0, df0, x, AA, BB, CC, DD, EE;
    f0  := evalf(f_in(0));
    df0 := evalf(df_in(0));

    AA := plot(f_in(x),x=0..amax);
    BB := plot([[0,f0],[amax,f0+amax*rho*df0]],color="LimeGreen");
    CC := plot([[0,f0],[amax,f0+amax*sigma*df0]],color="blue");
    DD := plot([[alpha,f_in(alpha)]],color="red",style = point);
    EE := plot([seq([[A_list[i],-i/5],[B_list[i],-i/5]],i=1..nops(A_list))],color="black");
    display(AA,BB,CC,DD,EE);
  end;

end module:
