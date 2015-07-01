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
    def __init__(self, backend):
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
        return (isinstance(other, LPBackendDictionary) and
                self._backend == other._backend)

    def basic_variables(self):
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
        if self._entering is None:
            raise ValueError("entering variable must be chosen to compute "
                             "its coefficients")
        k = tuple(self.nonbasic_variables()).index(self._entering)
        return k

    def leaving_coefficients(self):
        if self._leaving is None:
            raise ValueError("leaving variable must be chosen to compute "
                             "its coefficients")
        k = tuple(self.nonbasic_variables()).index(self._leaving)
        return k

    def nonbasic_variables(self):
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
        return col_nonbasics + row_nonbasics

    def objective_coefficients(self):
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
        return self._backend.get_objective_value()

    def get_backend(self):
        return self._backend

    def update(self):
        # ASK: todo?
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

