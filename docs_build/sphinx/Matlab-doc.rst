MATLAB toolbox
==============

The usage of `NLtoolbox` is easy:

Intantiate the class and select the test:

.. code-block:: matlab

    ref = NLtest();

    ref.select( 1 ); % select the first test.

    ref.select('Discrete boundary value function neq = 50');

the list of available names can be obtained as

.. code-block:: matlab

    res = ref.listall();

After the selection of the nonlinear system you can
evaluate it:

.. code-block:: matlab

    F = ref.evalF( x ); % evaluate F(x)

    JF = ref.evalJF( x ); % evaluate jacobian of F(x)

    P = ref.pattern(); % evaluate pasttern of the jacobian of F(x)

Some aux methods make easy to check you numerical solver
for nonlinear system


.. code-block:: matlab

    N     = ref.neq()
    NG    = ref.numGuess();
    G     = ref.guess( n );
    Ne    = ref.numExact();
    E     = ref.exact( n );
    [L,U] = ref.bbox();

See the API manual for details.

.. toctree::

  api-matlab/root.rst
