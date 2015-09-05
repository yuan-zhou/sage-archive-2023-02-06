# from sage.modules.all import random_vector, vector
import sage.numerical.backends.glpk_backend as backend
from sage.numerical.interactive_simplex_method import *
from sage.numerical.backends.glpk_backend import *
# from sage.symbolic.all import SR


def format(name, prefix, index):
    if name:
        return name.replace('[', '_').strip(']')
    else:
        return prefix + '_' + str(index)


class LPBackendDictionary(LPAbstractDictionary):
    r"""
    Construct a dictionary for an LP problem from an backend.

    INPUT:

        - ``backend`` -- the backend where the dictionary is
            constructed from

    OUTPUT:

       - a :class:`backend dictionary for an LP problem <LPBackendDictionary>`

    EXAMPLES:

    One needs an instance of :class:`MixedIntegerLinearProgram` to initialize
    this class::

        sage: p = MixedIntegerLinearProgram(maximization=True)
        sage: x = p.new_variable(integer=True, nonnegative=True)
        sage: p.add_constraint(-x[0] + x[1] <= 2)
        sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
        sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
        sage: b = p.get_backend()
        sage: d = LPBackendDictionary(b)
        sage: d
        LP problem dictionary (use typeset mode to see details)

    """
    def __init__(self, backend):
        r"""
        See :class:`LPBackendDictionary` for documentation.

        TESTS::

        sage: p = MixedIntegerLinearProgram(maximization=True)
        sage: x = p.new_variable(integer=True, nonnegative=True)
        sage: p.add_constraint(-x[0] + x[1] <= 2)
        sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
        sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
        sage: b = p.get_backend()
        sage: d = LPBackendDictionary(b)
        sage: TestSuite(d).run()

        """
        super(LPBackendDictionary, self).__init__()
        self._backend = backend
        col_vars = tuple(
            format(self._backend.col_name(i), 'x', i)
            for i in range(self._backend.ncols())
        )
        row_vars = tuple(
            format(self._backend.row_name(i), 'w', i)
            for i in range(self._backend.nrows())
        )
        names = ", ".join(col_vars + row_vars)
        self._R = PolynomialRing(self._backend.base_ring(),
                                 names, order="neglex")
        self._x = vector(self._R, self._R.gens())

    def __eq__(self, other):
        r"""
        Check if two LP problem dictionaries have the same
        reference.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``True`` if ``other`` is an :class:`LPDictionary` with all
          details the same as ``self``, ``False`` otherwise.

        TESTS:

        Setting up the problem::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: d = LPBackendDictionary(b)

        Test when two problems have the same reference:

            sage: d2 = d
            sage: d2 == d
            True

        Test when two problems have the same constrct:

            sage: d3 = copy(d)
            sage: d3 == d
            False

        """
        return (isinstance(other, LPBackendDictionary) and
                self._backend == other._backend)

    def basic_variables(self):
        r"""
        Return the basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES:

        Setting up the problem::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solver_parameter(
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()

        Use function in :class:`LPBackendDictionary`:

            sage: d = LPBackendDictionary(b)

        Use function in :class:`InteractiveLPProblem`:

            sage: lp, basis = p.interactive_linear_program()
            sage: lpd = lp.dictionary(*basis)

        Compare results:

            sage: d.basic_variables()
            (x_0, x_1)
            sage: lpd.basic_variables()
            (x_0, x_1)

        """
        col_basics = tuple(
            self._x[i]
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) == glp_bs
        )
        row_basics = tuple(
            self._x[i + self._backend.ncols()]
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) == glp_bs
        )
        return vector(col_basics + row_basics)

    def constant_terms(self):
        r"""
        Return the constant terms of relations of ``self``.

        OUTPUT:

        - a vector.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solver_parameter(
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            sage: d = LPBackendDictionary(b)
            sage: d.constant_terms()
            (1.3, 3.3)

        """
        col_const = tuple(
            self._backend.get_variable_value(i)
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) == glp_bs
        )
        row_const = tuple(
            self._backend.get_slack_value(i)
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) == glp_bs
        )
        return col_const + row_const

    def entering_coefficients(self):
        r"""
        Return coefficients of the entering variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: A = ([1, 1], [3, 1])
            sage: b = (1000, 1500)
            sage: c = (10, 5)
            sage: P = InteractiveLPProblemStandardForm(A, b, c)
            sage: D = P.initial_dictionary()
            sage: D.enter(1)
            sage: D.entering_coefficients()
            (1, 3)

        """
        if self._entering is None:
            raise ValueError("entering variable must be chosen to compute "
                             "its coefficients")
        k = tuple(self.nonbasic_variables()).index(self._entering)
        return k

    def leaving_coefficients(self):
        r"""
        Return coefficients of the leaving variable.

        OUTPUT:

        - a vector

        EXAMPLES::

        """
        if self._leaving is None:
            raise ValueError("leaving variable must be chosen to compute "
                             "its coefficients")
        k = tuple(self.nonbasic_variables()).index(self._leaving)
        return k

    def nonbasic_variables(self):
        r"""
        Return non-basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

        Return the basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES:

        Setting up the problem::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solver_parameter(
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()

        Use function in :class:`LPBackendDictionary`:

            sage: d = LPBackendDictionary(b)

        Use function in :class:`InteractiveLPProblem`:

            sage: lp, basis = p.interactive_linear_program()
            sage: lpd = lp.dictionary(*basis)

        Compare results:

            sage: d.nonbasic_variables()
            (w_0, w_1)
            sage: lpd.nonbasic_variables()
            (w_0, w_1)

        """
        col_nonbasics = tuple(
            self._x[i]
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) != glp_bs
        )
        row_nonbasics = tuple(
            self._x[i + self._backend.ncols()]
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) != glp_bs
        )
        return vector(col_nonbasics + row_nonbasics)

    def objective_coefficients(self):
        r"""
        Return coefficients of the objective of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES:

        Setting up the problem::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solver_parameter(
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()

        Use function in :class:`LPBackendDictionary`:

            sage: d = LPBackendDictionary(b)

        Use function in :class:`InteractiveLPProblem`:

            sage: lp, basis = p.interactive_linear_program()
            sage: lpd = lp.dictionary(*basis)

        Compare results:

            sage: d.objective_coefficients()
            (-0.58, -0.76)
            sage: lpd.objective_coefficients()
            (-0.5800000000000001, -0.76)

        """
        col_coefs = tuple(
            -self._backend.get_col_dual(i)
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) != glp_bs
        )
        row_coefs = tuple(
            -self._backend.get_row_dual(i)
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) != glp_bs
        )
        return col_coefs + row_coefs

    def objective_value(self):
        r"""
        Return the value of the objective value.

        OUTPUT:

        - a number

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solver_parameter(
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            sage: d = LPBackendDictionary(b)
            sage: d.objective_value()
            14.08

        """
        return self._backend.get_objective_value()

    def get_backend(self):
        r"""
        Return the backend used to create the dictionary.

        OUTPUT:

        - The corresponding dictionary

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solver_parameter(
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            sage: d = LPBackendDictionary(b)
            sage: d.get_backend()
            <sage.numerical.backends.glpk_backend.GLPKBackend object at \
            0x11a108e30>

        """
        return self._backend

    def update(self):
        r"""
        Update ``self`` using previously set entering and leaving variables.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solver_parameter(
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            sage: d = LPBackendDictionary(b)
            sage: d.objective_value()
            sage: d.enter(d.nonbasic_variables()[0])
            sage: d.leave(d.basic_variables()[0])
            sage: d.objective_value()
            sage: d.update()

        """
        entering = self._entering
        if entering is None:
            raise ValueError("entering variable must be set before updating")
        leaving = self._leaving
        if leaving is None:
            raise ValueError("leaving variable must be set before updating")


p = MixedIntegerLinearProgram(maximization=True)
x = p.new_variable(integer=True, nonnegative=True)

p.add_constraint(-x[0] + x[1] <= 2)
p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
p.set_objective(5.5 * x[0] + 2.1 * x[1])

print
print

print 'Through LPBackendDictionary()'
b = p.get_backend()
d = LPBackendDictionary(b)
print "Solving ......"
b.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
b.solve()
print 'basic vars:', d.basic_variables()
print 'nonbasic vars:', d.nonbasic_variables()
print 'constant terms:', d.constant_terms()
print 'obj coefs:', d.objective_coefficients()
print 'obj values:', d.objective_value()
print 'backend:', d.get_backend()

print
print

print 'Through interactive_linear_program()'
lp, basis = p.interactive_linear_program()
lpd = lp.dictionary(*basis)
print 'basic vars:', lpd.basic_variables()
print 'nonbasic vars:', lpd.nonbasic_variables()
print 'constant terms:', lpd.constant_terms()
print 'obj coefs:', lpd.objective_coefficients()
print 'obj values:', lpd.objective_value()

