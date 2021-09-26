# Plotting

Plotting := module()

  description "plotting";

  # Module defined as a package (i.e.) collection of procedures
  option package, load = ModuleLoad, unload = ModuleUnLoad;

  export plot_log, plot_iter;
  
  local ModuleLoad, ModuleUnLoad,
        EPSI, LLL, COLORS, CONTOURPARS;

  uses plots;

  ModuleUnLoad := proc()
    NULL;
  end proc:

  ModuleLoad := proc()
    NULL;
    EPSI        := 1e-50;
    LLL         := thickness=2,style=pointline;
    COLORS      := [ "blue", "ForestGreen", "black", "DeepPink", "Goldenrod" ];
    CONTOURPARS := filledregions = true, contours=50,
                   coloring = ["LimeGreen", "Goldenrod"],
                   transparency = 0.5;
  end proc;
  
  ModuleLoad();
  
  plot_log := proc( V_list )
    local i, iii, V, A, A_list;
    A_list := [];
    iii    := 1;
    for V in V_list do
      A      := logplot( [seq([i,abs(V[i])+EPSI],i=1..nops(V))],
                         LLL, color=COLORS[iii] );
      A_list := [op(A_list), A];
      iii    := iii+1;
    end;
    display(op(A_list));
  end proc:

  plot_iter := proc( epxr, XRANGE, YRANGE, X_list )
    local i, iii, X, AA, A, A_list:
    AA := contourplot( epxr, XRANGE, YRANGE, CONTOURPARS ):
    iii := 1;
    A_list := [];
    for X in X_list do
      A := plot( [seq([X[i][1],X[i][2]],i=1..nops(X))],
                  LLL, color=COLORS[iii] ):
      A_list := [op(A_list), A];
      iii    := iii+1;
      if iii > 4 then
        iii := 1;
      end;
    end;
    display( AA, op(A_list), scaling=constrained );
  end proc:

end module:
