from sage.numerical.interactive_simplex_method import *
from sage.numerical.backends.coin_backend import *
from sage.numerical.backends.general_backend_dictionary import *


class LPCoinBackendDictionary(LPAbstractBackendDictionary):
    r"""
    Construct a dictionary for an LP problem from an backend.

    INPUT:

        - ``backend`` -- the backend where the dictionary is
            constructed from

    OUTPUT:

       - a :class:`backend dictionary for an LP problem <LPCoinBackendDictionary>`

    EXAMPLES:

    # ASK: wierd description?
    One needs an instance of :class:`MixedIntegerLinearProgram` to initialize
    this class::

        sage: from sage.numerical.backends.coin_backend_dictionary \
              import LPCoinBackendDictionary
        sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
        sage: x = p.new_variable(nonnegative=True)
        sage: p.add_constraint(-x[0] + x[1] <= 2)
        sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
        sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
        sage: b = p.get_backend()
        sage: d = LPCoinBackendDictionary(b)
        sage: d
        LP problem dictionary (use typeset mode to see details)
    """
    def __init__(self, backend):
        r"""
        See :class:`LPCoinBackendDictionary` for documentation.

        TESTS::
            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: d = LPCoinBackendDictionary(b)
            sage: TestSuite(d).run(skip=['_test_pickling'])

        An exception will be raised if the problem is not in standard form
        i.e. with <= constraints and >= 0 variable bounds::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(8 * x[0] + 2 * x[1], min=17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: d = LPCoinBackendDictionary(b)
            Traceback (most recent call last):
            ...
            AttributeError: Problem constraints not in standard form.
        """
        super(LPCoinBackendDictionary, self).__init__(backend)

    # def __eq__(self, other):
    #     r"""
    #     Check if two LP problem dictionaries have the same
    #     reference.

    #     INPUT:

    #     - ``other`` -- anything

    #     OUTPUT:

    #     - ``True`` if ``other`` is an :class:`LPDictionary` with all
    #       details the same as ``self``, ``False`` otherwise.

    #     TESTS:

    #     Setting up the problem::

    #         sage: from sage.numerical.backends.coin_backend_dictionary \
    #               import LPCoinBackendDictionary
    #         sage: p = MixedIntegerLinearProgram(maximization=True,\
    #                                             solver="Coin")
    #         sage: x = p.new_variable(nonnegative=True)
    #         sage: p.add_constraint(-x[0] + x[1] <= 2)
    #         sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
    #         sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
    #         sage: b = p.get_backend()
    #         sage: d = LPCoinBackendDictionary(b)

    #     Test when two problems have the same reference:

    #         sage: d2 = d
    #         sage: d2 == d
    #         True

    #     Test when two problems have the same constrct:

    #         sage: d3 = LPCoinBackendDictionary(copy(p).get_backend())
    #         sage: d3 == d
    #         False
    #     """
    #     return (isinstance(other, LPCoinBackendDictionary) and
    #             self._backend == other._backend)

    def basic_variables(self):
        r"""
        Return the basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPCoinBackendDictionary(b)
            sage: d.basic_variables()
            (x_0, x_1)
        """
        col_stat, row_stat = self._backend.get_basis_status()
        col_basics = tuple(
            self._x[i]
            for i in range(self._backend.ncols())
            if col_stat[i] == 1
        )
        row_basics = tuple(
            self._x[i+self._backend.ncols()]
            for i in range(self._backend.nrows())
            if row_stat[i] == 1
        )
        return vector(col_basics + row_basics)

    def constant_terms(self):
        r"""
        Return the constant terms of relations of ``self``.

        OUTPUT:

        - a vector.

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPCoinBackendDictionary(b)
            sage: d.constant_terms()
            (1.3, 3.3)
        """
        col_stat, row_stat = self._backend.get_basis_status()
        col_const = tuple(
            self._backend.get_variable_value(i)
            for i in range(self._backend.ncols())
            if col_stat[i] == 1
        )
        row_const = tuple(
            self._backend.row_bounds(i)[1]
            - self._backend.get_variable_value(self._backend.ncols()+i)
            for i in range(self._backend.nrows())
            if row_stat[i] == 1
        )
        return vector(col_const + row_const)

    def entering_coefficients(self):
        r"""
        Return coefficients of the entering variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
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
            sage: vars = d.nonbasic_variables()
            sage: vars
            (x_0, x_1, w_0, w_2)
            sage: d.enter(vars[0])
            sage: d.entering_coefficients()
            (36.0, 26.0, 5.0)
            sage: d.enter(vars[1])
            sage: d.entering_coefficients()
            (1.0, 2.0, 0.0)
        """
        if self._entering is None:
            raise ValueError("entering variable must be chosen to compute "
                             "its coefficients")

        index = self._x.index(self._entering)

        return vector(self._backend.get_binva_col(index))

    def leaving_coefficients(self):
        r"""
        Return coefficients of the leaving variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
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
            sage: vars = d.basic_variables()
            sage: vars
            (x_2, x_3, w_1)
            sage: d.leave(vars[0])
            sage: d.leaving_coefficients()
            (5.0, 0.0, 0.0, 1.0)
            sage: d.leave(vars[1])
            sage: d.leaving_coefficients()
            (36.0, 1.0, 1.0, 7.0)
        """
        if self._leaving is None:
            raise ValueError("leaving variable must be chosen to compute "
                             "its coefficients")

        var_index = self._x.index(self._leaving)
        row_indices = self._backend.get_basics()
        row_index = tuple(row_indices).index(var_index)
        from sage.misc.flatten import flatten
        row = flatten(self._backend.get_binva_row(row_index))
        nonbasic_indicies = [self._x.index(v) for v in self.nonbasic_variables()]
        return vector([row[i] for i in nonbasic_indicies])

    def nonbasic_variables(self):
        r"""
        Return non-basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPCoinBackendDictionary(b)
            sage: d.nonbasic_variables()
            (w_0, w_1)
        """
        col_stat, row_stat = self._backend.get_basis_status()
        col_nonbasics = tuple(
            self._x[i]
            for i in range(self._backend.ncols())
            if col_stat[i] != 1
        )
        row_nonbasics = tuple(
            self._x[i+self._backend.ncols()]
            for i in range(self._backend.nrows())
            if row_stat[i] != 1
        )
        return vector(col_nonbasics + row_nonbasics)

    def objective_coefficients(self):
        r"""
        Return coefficients of the objective of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
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
            sage: d.objective_coefficients()
            (486.0, 10.0, 13.0, 95.0)
        """
        col_stat, row_stat = self._backend.get_basis_status()
        cost = self._backend.get_reduced_cost()
        price = self._backend.get_row_price()
        col_coefs = tuple(
            cost[i]
            for i in range(self._backend.ncols())
            if col_stat[i] != 1
        )
        row_coefs = tuple(
            -price[i]
            for i in range(self._backend.nrows())
            if row_stat[i] != 1
        )
        return vector(col_coefs + row_coefs)

    # def objective_value(self):
    #     r"""
    #     Return the value of the objective value.

    #     OUTPUT:

    #     - a number

    #     EXAMPLES::

    #         sage: from sage.numerical.backends.coin_backend_dictionary \
    #               import LPCoinBackendDictionary
    #         sage: p = MixedIntegerLinearProgram(maximization=True,\
    #                                             solver="Coin")
    #         sage: x = p.new_variable(nonnegative=True)
    #         sage: p.add_constraint(-x[0] + x[1] <= 2)
    #         sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
    #         sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
    #         sage: b = p.get_backend()
    #         sage: b.solve()
    #         0
    #         sage: d = LPCoinBackendDictionary(b)
    #         sage: d.objective_value()
    #         14.08
    #     """
    #     return self._backend.get_objective_value()

    # def get_backend(self):
    #     r"""
    #     Return the backend used to create the dictionary.

    #     OUTPUT:

    #     - The corresponding dictionary

    #     EXAMPLES::

    #         sage: from sage.numerical.backends.coin_backend_dictionary \
    #               import LPCoinBackendDictionary
    #         sage: p = MixedIntegerLinearProgram(maximization=True,\
    #                                             solver="Coin")
    #         sage: x = p.new_variable(nonnegative=True)
    #         sage: p.add_constraint(-x[0] + x[1] <= 2)
    #         sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
    #         sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
    #         sage: b = p.get_backend()
    #         sage: b.solve()
    #         0
    #         sage: d = LPCoinBackendDictionary(b)
    #         sage: d.get_backend()
    #         <sage.numerical.backends.coin_backend.CoinBackend object at ...>
    #     """
    #     return self._backend

    def update(self):
        r"""
        Update ``self`` using previously set entering and leaving variables.

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
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
            sage: d.objective_value()
            1331.0
            sage: d.nonbasic_variables()
            (x_0, x_1, w_0, w_2)
            sage: d.enter(d.nonbasic_variables()[0])
            sage: d.basic_variables()
            (x_2, x_3, w_1)
            sage: d.leave(d.basic_variables()[0])
            sage: d.objective_value()
            1331.0
            sage: d.update()
            sage: d.basic_variables()
            (x_0, x_3, w_1)
            sage: d.nonbasic_variables()
            (x_1, x_2, w_0, w_2)
            sage: d.objective_value()
            261.8

        TESTS:

        An error will be raised if the pivot selected is zero::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
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
            sage: d.leave(d.basic_variables()[0])
            sage: d.leaving_coefficients()
            (5.0, 0.0, 0.0, 1.0)
            sage: d.enter(d.nonbasic_variables()[1])
            sage: d.update()
            Traceback (most recent call last):
            ...
            ValueError: incompatible choice of entering and leaving variables
        """
        entering = self._entering
        if entering is None:
            raise ValueError("entering variable must be set before updating")

        leaving = self._leaving
        if leaving is None:
            raise ValueError("leaving variable must be set before updating")

        matching_index = tuple(self.nonbasic_variables()).index(entering)
        coef = self.leaving_coefficients()[matching_index]
        if coef == 0:
            raise ValueError("incompatible choice of entering and leaving "
                             "variables")

        col_stat, row_stat = self._backend.get_basis_status()
        entering_index = self._x.index(entering)
        leaving_index = self._x.index(leaving)

        if entering_index < self._backend.ncols():
            col_stat[entering_index] = 1
        else:
            row_stat[entering_index-self._backend.ncols()] = 1

        if leaving_index < self._backend.ncols():
            col_stat[leaving_index] = 3
        else:
            row_stat[leaving_index-self._backend.ncols()] = 3

        self._backend.set_basis_status(col_stat, row_stat)

    def add_row(self, nonbasic_coef, constant, slack_variable, integer_slack_variable=False):
        r"""
        Update a dictionary with an additional row based on a given dictionary.

        INPUT:

        - ``nonbasic_coef``-- a list of nonbasic coefficients for the new row

        - ``constant``-- a number of the constant term for the new row

        - ``slack_variable``-- a string of the name for the new slack variable

        - ``integer_slack_variable``-- (default: false) a boolean value
        indicating if the new slack variable is integer or not.

        EXAMPLES::

            sage: from sage.numerical.backends.coin_backend_dictionary \
                  import LPCoinBackendDictionary
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
            sage: d.basic_variables()
            (x_2, x_3, w_1)
            sage: d.nonbasic_variables()
            (x_0, x_1, w_0, w_2)
            sage: d.objective_coefficients()
            (486.0, 10.0, 13.0, 95.0)
            sage: d.add_row(range(3,7), 2, 'z_0')
            sage: d.objective_coefficients()
            (486.0, 10.0, 13.0, 95.0)
            sage: d.basic_variables()
            (x_2, x_3, w_1, z_0)
            sage: d.leave(d.basic_variables()[3])
            sage: d.leaving_coefficients()
            (3.0, 4.0, 5.0, 6.0)
            sage: b.solve()
            0
            sage: d.basic_variables()
            (x_2, x_3, w_1, z_0)
            sage: d.nonbasic_variables()
            (x_0, x_1, w_0, w_2)

        Variables have 0 as their coefficient will not show up in the
        tableau:

            sage: d.add_row(range(-2, 2), 5, 'z_1')
            sage: d.get_backend().row(4)
            ([0, 1, 2], [-7.0, -1.0, -1.0])
        """
        if len(nonbasic_coef) != self._backend.ncols():
            raise ValueError("Length of nonbasic coefficients incompatible")

        # Convert to problem variable coefficients
        coefs = [0] * self._backend.ncols()
        for i, var in enumerate(self.nonbasic_variables()):
            index = self._x.index(var)
            if index < self._backend.ncols():
                coefs[index] += nonbasic_coef[i]
            else:
                row_pos = index-self._backend.ncols()
                row = self._backend.row(row_pos)
                for j, v in zip(*row):
                    coefs[j] -= nonbasic_coef[i] * v
                upper_bound = self._backend.row_bounds(row_pos)[1]
                constant -= nonbasic_coef[i] * upper_bound

        coef_pairs = [(i, coefs[i]) for i in range(self._backend.ncols())
                                    if coefs[i] != 0]        
        self._backend.add_linear_constraint(
            coef_pairs, None, constant, slack_variable)

        # print self._backend.get_reduced_cost(), self._backend.get_row_price()
        # Update buffered variables
        # TODO: extract small functions here
        self._names += ', '
        self._names += LPCoinBackendDictionary._format_(
                            name=self._backend.row_name(self._backend.nrows()-1),
                            symbol='w', index=self._backend.nrows()-1
                        )
        self._R = PolynomialRing(self._backend.base_ring(),
                                 self._names, order="neglex")
        self._x = tuple(self._R.gens())
        # print self._backend.get_reduced_cost(), self._backend.get_row_price()

        # Update basis status in the backend
        curr_basis = self._backend.get_basis_status()
        curr_basis[1].append(1)
        self._backend.set_basis_status(*curr_basis)
        # print self._backend.get_reduced_cost(), self._backend.get_row_price()
        # TODO: trim print statements

    # ASK: make this function general as well?
    # the problem is abstract backend dictionary will not be able to 
    # be viewed, so not consistent in the doc strings
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

# from sage.numerical.backends.coin_backend_dictionary \
#       import LPCoinBackendDictionary
# from sage.numerical.mip import MixedIntegerLinearProgram
# p = MixedIntegerLinearProgram(maximization=True,\
#                                     solver="Coin")
# x = p.new_variable(nonnegative=True)
# p.add_constraint(x[0] + x[1] - 7*x[2] + x[3] <= 22)
# p.add_constraint(x[1] + 2*x[2] - x[3] <= 13)
# p.add_constraint(5*x[0] + x[2] <= 11)
# p.set_objective(2*x[0] + 3*x[1] + 4*x[2] + 13*x[3])
# b = p.get_backend()
# b.solve()
# d = LPCoinBackendDictionary(b)
# view(d.dictionary())
# d.add_row(range(3,7), 2, 'z_0')
# view(d.dictionary())
# # b.solve()
# d.leave(d.basic_variables()[3])
# print d.leaving_coefficients()
# from sage.numerical.backends.generic_backend import get_solver
# p = get_solver(solver = "Coin")
# p.add_variables(3)
# p.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48)
# p.add_linear_constraint(zip([0, 1, 2], [4, 2, 1.5]), None, 20)
# p.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8)
# p.set_objective([60, 30, 20])
# p.set_basis_status([3, 3, 3], [1, 1, 1])
# p.get_row_price(), p.get_reduced_cost()
# p.set_basis_status([1, 3, 3], [1, 1, 3])
# p.get_row_price(), p.get_reduced_cost()
# p.add_linear_constraint(zip([0, 1, 2], [7, 3, 2]), None, 9)
# p.get_row_price(), p.get_reduced_cost()
# p.get_basis_status()
# p.set_basis_status([1, 3, 3], [1, 1, 3, 1])
# p.get_row_price(), p.get_reduced_cost()
