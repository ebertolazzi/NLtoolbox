/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Function15 : public nonlinearSystem {
  typedef pair<int_type,int_type> INDEX;
  mutable map<INDEX,real_type> jac_idx_vals;
public:

  Function15( int_type neq )
  : nonlinearSystem(
      "Function 15",
      "@article{LaCruz:2003,\n"
      "  author    = { William {La Cruz}  and  Marcos Raydan},\n"
      "  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems},\n"
      "  journal   = {Optimization Methods and Software},\n"
      "  year      = {2003},\n"
      "  volume    = {18},\n"
      "  number    = {5},\n"
      "  pages     = {583--599},\n"
      "  publisher = {Taylor & Francis},\n"
      "  doi       = {10.1080/10556780310001610493},\n"
      "}\n",
      neq
    )
  {
    checkMinEquations(n,2);
    jac_idx_vals.clear();
    jac_idx_vals[INDEX(0,0)]     = 1;
    jac_idx_vals[INDEX(n-1,n-1)] = 1;
    jac_idx_vals[INDEX(n-1,n-2)] = 1;
    for ( int_type i = 1; i < n-1; ++i ) {
      jac_idx_vals[INDEX(i,i-1)] = 1;
      jac_idx_vals[INDEX(i,i)]   = 1;
      jac_idx_vals[INDEX(i,i+1)] = 1;
    }
    for ( int_type i = 0; i < n; ++i ) {
      jac_idx_vals[INDEX(i,n-5)] = 1;
      jac_idx_vals[INDEX(i,n-4)] = 1;
      jac_idx_vals[INDEX(i,n-3)] = 1;
      jac_idx_vals[INDEX(i,n-2)] = 1;
      jac_idx_vals[INDEX(i,n-1)] = 1;
    }
  }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type bf = 3*x(n-5) - x(n-4) - x(n-3) + 0.5 * x(n-2) - x(n-1) +1;
    if ( i == 0   ) return -2*x(0)*x(0) + 3*x(0) + bf;
    if ( i == n-1 ) return -2*x(n-1)*x(n-1) + 3*x(n-1) - x(n-2) + bf;
    return -2*x(i)*x(i) + 3*x(i) - x(i-1) - 2*x(i+1) + bf;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type bf = 3*x(n-5) - x(n-4) - x(n-3) + 0.5 * x(n-2) - x(n-1) +1;
    f(0)   = -2*x(0)*x(0)     + 3*x(0)            + bf;
    f(n-1) = -2*x(n-1)*x(n-1) + 3*x(n-1) - x(n-2) + bf;
    for ( int_type i = 1; i < n-1; ++i )
      f(i) = -2*x(i)*x(i) + 3*x(i) - x(i-1) - 2*x(i+1) + bf;
  }

  virtual
  int_type
  jacobianNnz() const override
  { return int_type(jac_idx_vals.size()); }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    map<pair<int_type,int_type>,real_type>::const_iterator it = jac_idx_vals.begin();
    int_type kk = 0;
    for (; it != jac_idx_vals.end(); ++it ) {
      ii(kk) = it->first.first;
      jj(kk) = it->first.second;
      ++kk;
    }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    map<pair<int_type,int_type>,real_type>::iterator it = jac_idx_vals.begin();
    for (; it != jac_idx_vals.end(); ++it ) it->second = 0;

    jac_idx_vals[INDEX(0,0)]     = -4*x(0)   + 3;
    jac_idx_vals[INDEX(n-1,n-1)] = -4*x(n-1) + 3;
    jac_idx_vals[INDEX(n-1,n-2)] = -1;
    for ( int_type i = 1; i < n-1; ++i ) {
      jac_idx_vals[INDEX(i,i-1)] = -1;
      jac_idx_vals[INDEX(i,i)]   =  -4*x(i) + 3;
      jac_idx_vals[INDEX(i,i+1)] = -2;
    }
    for ( int_type i = 1; i < n-1; ++i ) {
      jac_idx_vals[INDEX(i,n-5)] += 3;
      jac_idx_vals[INDEX(i,n-4)] -= 1;
      jac_idx_vals[INDEX(i,n-3)] -= 1;
      jac_idx_vals[INDEX(i,n-2)] += 0.5;
      jac_idx_vals[INDEX(i,n-1)] -= 1;
    }
    int_type kk = 0;
    for ( it = jac_idx_vals.begin(); it != jac_idx_vals.end(); ++it )
      jac(kk++) = it->second;
  }

  virtual
  int_type
  numExactSolution() const override { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(-1);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
