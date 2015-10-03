from sage.numerical.interactive_simplex_method import *

# ASK: change file name to abstract_backend_dictionary?
class LPAbstractBackendDictionary(LPAbstractDictionary):
    r"""
    Construct an abstract dictionary for an LP problem from an backend.

    INPUT:

        - ``backend`` -- the backend where the dictionary is
            constructed from

    OUTPUT:

       - a :class:`backend dictionary for an LP problem <LPAbstractBackendDictionary>`

    EXAMPLES:

    One needs an instance of :class:`MixedIntegerLinearProgram` to initialize
    this class::

        sage: from sage.numerical.backends.general_backend_dictionary \
              import LPAbstractBackendDictionary
        sage: p = MixedIntegerLinearProgram(maximization=True)
        sage: x = p.new_variable(nonnegative=True)
        sage: p.add_constraint(-x[0] + x[1] <= 2)
        sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
        sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
        sage: b = p.get_backend()
        sage: d = LPAbstractBackendDictionary(b)
        sage: d
        LP problem dictionary (use typeset mode to see details)
    """
    # ASK: name it like "LP Abstract Backend Dictionary"?
    # TODO: doc strings using itself
    def __init__(self, backend):
        super(LPAbstractBackendDictionary, self).__init__()
        self._backend = backend
        for i in range(self._backend.nrows()):
            if self._backend.row_bounds(i)[0] != None \
               or self._backend.row_bounds(i)[1] == None:
                raise AttributeError("Problem constraints "
                                     "not in standard form.")

        for i in range(self._backend.ncols()):
            if self._backend.variable_lower_bound(i) == None:
                raise AttributeError("Problem variables "
                                     "not in standard form.")

        col_vars = tuple(
            LPAbstractBackendDictionary._format_(
                name=self._backend.col_name(i), 
                symbol='x', index=i
            )
            for i in range(self._backend.ncols())
        )
        row_vars = tuple(
            LPAbstractBackendDictionary._format_(
                name=self._backend.row_name(i), 
                symbol='w', index=i
            )
            for i in range(self._backend.nrows())
        )
        self._names = ", ".join(col_vars + row_vars)
        self._R = PolynomialRing(self._backend.base_ring(),
                                 self._names, order="neglex")
        self._x = self._R.gens() # Note: made this a tuple by default

    def __eq__(self, other):
        r"""
        Check if two LP problem dictionaries have the same
        reference.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``True`` if ``other`` is an :class:`LPDictionary` with its
          reference the same as ``self``, ``False`` otherwise.

        TESTS:

        Setting up the problem::

            sage: from sage.numerical.backends.general_backend_dictionary \
                import LPAbstractBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: d = LPAbstractBackendDictionary(b)

        Test when two problems have the same reference:

            sage: d2 = d
            sage: d2 == d
            True

        Test when two problems have the same constrct:

            sage: d3 = LPAbstractBackendDictionary(copy(p).get_backend())
            sage: d3 == d
            False
        """
        # ASK: should check different thing if not the same class?
        return (isinstance(other, LPAbstractBackendDictionary) and
                self._backend == other._backend)

    @staticmethod # NOTE: made this static
    def _format_(name='', symbol='x', infix='_', index='0'):
        r"""
        Returns a proper name for a given parameters.

        INPUT::

        - ``name`` -- (defualt: '') the original name of the variable

        -``symbol`` -- (defualt: 'x') the symbol for the new name

        -``infix`` -- (default: '_') the character separate the symbol
        and its index for the new name

        -``index`` -- (default: '0') the index following the infix for the new name

        TESTS:

        If a name is given, then returns the name itself::

            sage: from sage.numerical.backends.general_backend_dictionary \
                  import LPAbstractBackendDictionary
            sage: LPAbstractBackendDictionary._format_('Name')
            'Name'

        However, if the name given is in the form 'symbol[index]', then it will be
        converted to 'symbol+infix+index' i.e. 'symbol_index' by default:

            sage: LPAbstractBackendDictionary._format_('x[3]')
            'x_3'
            sage: LPAbstractBackendDictionary._format_('x[2]', infix='~')
            'x~2'

        If no name is given, then a newly create name will be returned in the
        form of symbol+infix+index:

            sage: LPAbstractBackendDictionary._format_(symbol='w', index='7')
            'w_7'
        """
        if name:
            return name.replace('[', infix).strip(']')
        else:
            return symbol + infix + str(index)

    def objective_value(self):
        r"""
        Return the value of the objective value.

        OUTPUT:

        - a number

        EXAMPLES::

            sage: from sage.numerical.backends.general_backend_dictionary \
                  import LPAbstractBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPAbstractBackendDictionary(b)
            sage: d.objective_value()
            14.08
        """
        return self._backend.get_objective_value()

    def get_backend(self):
        r"""
        Return the backend used to create the dictionary.

        OUTPUT:

        - The corresponding backend associated with self

        EXAMPLES::

            sage: from sage.numerical.backends.general_backend_dictionary \
                  import LPAbstractBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPAbstractBackendDictionary(b)
            sage: d.get_backend()
            <sage.numerical.backends.coin_backend.CoinBackend object at ...>
        """
        return self._backend

    def dictionary(self):
        r"""
        Return a regular LP dictionary matching ``self``.

        OUTPUT:

        - an :class:`LP dictionary <LPDictionary>`

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
            sage: from sage.numerical.mip import MixedIntegerLinearProgram
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0] + x[1] - 7*x[2] + x[3] <= 22)
            sage: p.add_constraint(x[1] + 2*x[2] - x[3] <= 13)
            sage: p.add_constraint(5*x[0] + x[2] <= 11)
            sage: p.set_objective(2*x[0] + 3*x[1] + 4*x[2] + 13*x[3])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPCoinBackendDictionary(b)
            sage: view(d.dictionary()) # not tested

        # ASK: greek sign fine here?
        ζ is used as default problem name, and it can be changed:

            sage: b.problem_name("beta")
            sage: view(d.dictionary()) # not tested
        """
        rows = []
        for i in range(self._backend.nrows()):
            self._leaving = self.basic_variables()[i]
            rows.append(self.leaving_coefficients())
        import sage.matrix.constructor as matrix
        m = matrix.matrix(rows)
        name = self._backend.problem_name()
        if not name:
            name = 'zeta'
        D = LPDictionary(m,
                         self.constant_terms(),
                         self.objective_coefficients(),
                         self.objective_value(),
                         self.basic_variables(),
                         self.nonbasic_variables(),
                         name) # ASK: default name of the problem?
        D._entering = self._entering
        D._leaving = self._leaving
        return D

# ASK: general functions in backends?
# just find common backend functions in coin and glpk?