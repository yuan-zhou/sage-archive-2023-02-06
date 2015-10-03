from sage.numerical.interactive_simplex_method import *


class LPAbstractBackendDictionary(LPAbstractDictionary):
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
            self._format_(self._backend.col_name(i), 'x', i)
            for i in range(self._backend.ncols())
        )
        row_vars = tuple(
            self._format_(self._backend.row_name(i), 'w', i)
            for i in range(self._backend.nrows())
        )
        self._names = ", ".join(col_vars + row_vars)
        self._R = PolynomialRing(self._backend.base_ring(),
                                 self._names, order="neglex")
        self._x = tuple(self._R.gens()) # TODO: make this a tuple by default

    def __eq__(self, other):
        return (isinstance(other, LPBackendDictionary) and
                self._backend == other._backend)

    def _format_(self, name, prefix, index):
        if name:
            return name.replace('[', '_').strip(']')
        else:
            return prefix + '_' + str(index)

    def objective_value(self):
        return self._backend.get_objective_value()

    def get_backend(self):
        return self._backend

# ASK: general functions in backends?
# just find common backend functions in coin and glpk?