StrongWolfe := module()

  description "StrongWolfe linesearch";

  # Module defined as a package (i.e.) collection of procedures
  option package, load = ModuleLoad, unload = ModuleUnLoad;

  export lines, linesearch, interpolate, zoom, esatta, AB_list, doplot, doplot2;
  
  local ModuleLoad,
        ModuleUnLoad,
        _debug,
        c_1_vec, c_2_vec, rho, max_iterations,
        Dphi_0, phi_0, A_list, B_list;

  uses LinearAlgebra;

  ModuleUnLoad := proc()
    printf("StrongWolfe unload\n");
    NULL;
  end proc:

  ModuleLoad := proc()
    printf("StrongWolfe load\n");
    _debug         := true;
    c_1_vec        := [0.1,0.01,1e-4]; # Armijo condition
    c_2_vec        := [0.1,0.5,0.9];   # second (strong) Wolfe condition
    rho            := 2;  # bracket growth
    max_iterations := [5,5,10];
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

  lines := proc( f0, df0, L )
    [ [0,f0], [L,f0+L*rho*df0] ], [ [0,f0], [L,f0+L*sigma*df0] ];
  end proc;

  # a_lo = a_{i - 1}
  # a_hi = a_{i}
  interpolate := proc ( a_lo, a_hi, phi_lo, phi_hi, Dphi_lo, Dphi_hi )
    local d1, d2;
    d1 := Dphi_lo + Dphi_hi - 3 * (phi_lo - phi_hi) / (a_lo - a_hi);
    d2 := sqrt(d1 * d1 - Dphi_lo * Dphi_hi);
    return a_hi - (a_hi - a_lo) * ((Dphi_hi + d2 - d1) / (Dphi_hi - Dphi_lo + 2 * d2))
  end:

  zoom := proc( a_lo_in, a_hi_in, phi, Dphi )
    local kkk, iteration,
          c1_Dphi_0, c2_Dphi_0,
          a_lo, a_hi, a_j, 
          phi_lo, phi_hi, phi_j, 
          Dphi_lo, Dphi_hi, Dphi_j;

    a_lo := a_lo_in; phi_lo := phi(a_lo); Dphi_lo := Dphi(a_lo);
    a_hi := a_hi_in; phi_hi := phi(a_hi); Dphi_hi := Dphi(a_hi);

    # Shrink bracket
    for kkk from 1 to 3 do
      c1_Dphi_0 := c_1_vec[kkk]*Dphi_0;
      c2_Dphi_0 := c_2_vec[kkk]*Dphi_0;
      for iteration from 1 to max_iterations[kkk] do

        A_list := [ op(A_list), a_lo ];
        B_list := [ op(B_list), a_hi ];

        # Interpolate a_j
        if a_lo < a_hi then
          a_j := interpolate( a_lo, a_hi, phi_lo, phi_hi, Dphi_lo, Dphi_hi );
        else
          # TODO: Check if this is needed
          a_j := interpolate( a_hi, a_lo, phi_hi, phi_lo, Dphi_hi, Dphi_lo );
        end;

        # Evaluate phi(a_j)
        phi_j  := phi(a_j);
        Dphi_j := Dphi(a_j);

        # Check Armijo
        if phi_j > phi_0 + a_j * c1_Dphi_0 or phi_j > phi_lo then
          a_hi    := a_j;
          phi_hi  := phi_j;
          Dphi_hi := Dphi_j;
        else

          if abs(Dphi_j) <= -c2_Dphi_0 then
            return a_j;
          end;

          if Dphi_j * (a_hi - a_lo) >= 0 then
            a_hi    := a_lo;
            phi_hi  := phi_lo;
            Dphi_hi := Dphi_lo; 
          end;

          a_lo    := a_j;
          phi_lo  := phi_j;
          Dphi_lo := Dphi_j; 
        end;
      end;
    end;

    # Quasi-error response
    return a_j;
  end;

  # `StrongWolfe`: This linesearch algorithm guarantees that 
  # the step length satisfies the (strong) Wolfe conditions.
  # See Nocedal and Wright - Algorithms 3.5 and 3.6
  # This algorithm is mostly of theoretical interest,
  # users should most likely
  # `MoreThuente`, `HagerZhang` or `BackTracking`.

  linesearch := proc( phi, Dphi, alpha0 )
    local kkk, iteration, c1_Dphi_0, c2_Dphi_0,
          i, a_0, a_i, a_im1, phi_i, Dphi_i, phi_im1, a_star;

    phi_0  := phi(0);
    Dphi_0 := Dphi(0);

    A_list := [];
    B_list := [];

    # Step-sizes
    a_0     := 0;
    a_im1   := a_0;
    a_i     := alpha0;
    phi_im1 := phi_0;

    for kkk from 1 to 3 do
      c1_Dphi_0 := c_1_vec[kkk]*Dphi_0;
      c2_Dphi_0 := c_2_vec[kkk]*Dphi_0;
      for iteration from 1 to max_iterations[kkk] do

        phi_i := phi(a_i);

        # Test Wolfe conditions
        if (phi_i > phi_0 + a_i * c1_Dphi_0) or
           (phi_i >= phi_im1 and (iteration > 1 or kkk > 1) ) then
          a_star := zoom( a_im1, a_i, phi, Dphi );
          return a_star, phi(a_star);
        end;

        Dphi_i := Dphi(a_i);

        # Check condition 2
        if abs(Dphi_i) <= -c2_Dphi_0 then
          return a_i, phi_i;
        end;

        # Check condition 3
        if Dphi_i >= 0 then # FIXME untested!
          a_star := zoom( a_i, a_im1, phi, Dphi );
          return a_star, phi(a_star);
        end;

        A_list := [ op(A_list), a_im1 ];
        B_list := [ op(B_list), a_i ];

        # Choose a_iplus1 from the interval (a_i, a_max)
        a_im1 := a_i;
        a_i   := a_i * rho;

        # Update phi_im1
        phi_im1 := phi_i;
      end;
    end;

    # Quasi-error response TODO make this error instead
    printf("\n\nlinesearch failed!\n\n");
    return a_im1, phi(a_im1);
  end:

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
