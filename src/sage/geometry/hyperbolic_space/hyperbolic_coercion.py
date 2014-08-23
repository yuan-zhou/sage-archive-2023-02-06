"""
Coercion Maps Between Hyperbolic Plane Models

This module implements the coercion maps between different hyperbolic
plane models.

AUTHORS:

- Travis Scrimshaw (2014): initial version
"""

#***********************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPoint
from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometry
from sage.categories.morphism import Morphism
from sage.symbolic.pynac import I
from sage.matrix.constructor import matrix
from sage.rings.all import RR
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.misc.lazy_import import lazy_import

class HyperbolicModelCoercion(Morphism):
    """
    Abstract base class for morphisms between the hyperbolic models.
    """
    def _repr_type(self):
        """
        Return the type of morphism.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi._repr_type()
            'Coercion Isometry'
        """
        return "Coercion Isometry"

    def _call_(self, x):
        """
        Return the image of ``x`` under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi(PD.get_point(0.5+0.5*I))
            Point in UHP 2.00000000000000 + 1.00000000000000*I

        It is an error to try to convert a boundary point to a model
        that doesn't support boundary points::

            sage: HM = HyperbolicPlane().HM()
            sage: psi = HM.coerce_map_from(UHP)
            sage: psi(UHP(infinity)
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented for the Hyperbolid Model
        """
        C = self.codomain()
        bdry = False
        if C.is_bounded():
            if self.domain().is_bounded():
                bdry = x.is_boundary()
            else:
                bdry = C.bdry_point_in_model(x)
        elif self.domain().is_bounded() and x.is_boundary():
            raise NotImplementedError("boundary points are not implemented for"
                                      " the {0}".format(C.name()))
        return C.element_class(C, self.image_point(x.coordinates()), bdry,
                               check=False, **x.graphics_options())

    def convert_geodesic(self, x):
        """
        Convert the geodesic ``x`` of the domain into a geodesic of
        the codomain.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi.convert_geodesic(PD.get_geodesic(0.5+0.5*I, -I))
            Geodesic in UHP from 2.00000000000000 + 1.00000000000000*I to 0
        """
        return self.codomain().get_geodesic(self(x.start()), self(x.end()),
                                            **x.graphics_options())

    def convert_isometry(self, x):
        """
        Convert the hyperbolic isometry ``x`` of the domain into an
        isometry of the codomain.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi(PD.get_point(0.5+0.5*I))
            Point in UHP 2.00000000000000 + 1.00000000000000*I
        """
        return self.codomain().get_isometry(self.image_isometry(x._matrix))

    def __invert__(self):
        """
        Return the inverse coercion of ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = UHP.coerce_map_from(PD)
            sage: ~phi
            Coercion Isometry morphism:
              From: Hyperbolic plane in the Upper Half Plane Model model
              To:   Hyperbolic plane in the Poincare Disk Model model
        """
        return self.domain().coerce_map_from(self.codomain())

############
# From UHP #
############

class CoercionUHPtoPD(HyperbolicModelCoercion):
    """
    Coercion from the UHP to PD model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: phi = PD.coerce_map_from(UHP)
            sage: phi.image_point(I)
            Point in PD 0
        """
        if x == infinity:
            return I
        return (x - I)/(Integer(1) - I*x)

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return matrix(2,[1,-I,-I,1]) * x * matrix(2,[1,I,I,1])/Integer(2)

class CoercionUHPtoKM(HyperbolicModelCoercion):
    """
    Coercion from the UHP to KM model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.point_to_model(3 + I, 'KM')
            (6/11, 9/11)
        """
        if p == infinity:
            return (0, 1)
        return ((2*real(p))/(real(p)**2 + imag(p)**2 + 1),
                (real(p)**2 + imag(p)**2 - 1)/(real(p)**2 + imag(p)**2 + 1))

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return SL2R_to_SO21(x)

class CoercionUHPtoHM(HyperbolicModelCoercion):
    """
    Coercion from the UHP to HM model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.point_to_model(3 + I, 'HM')
            (3, 9/2, 11/2)
        """
        return vector((real(x)/imag(x),
                      (real(x)**2 + imag(x)**2 - 1)/(2*imag(x)),
                      (real(x)**2 + imag(x)**2 + 1)/(2*imag(x))))

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return SL2R_to_SO21(x)

###########
# From PD #
###########

class CoercionPDtoUHP(HyperbolicModelCoercion):
    """
    Coercion from the PD to UHP model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: UHP = HyperbolicPlane().UHP()
            sage: phi = UHP.coerce_map_from(PD)
            sage: phi.image_point(0.5+0.5*I)
            2.00000000000000 + 1.00000000000000*I
            sage: phi.image_point(0)
            I
        """
        if x == I:
            return infinity
        return (x + I)/(Integer(1) + I*x)

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return matrix(2,[1,I,I,1]) * x * matrix(2,[1,-I,-I,1])/Integer(2)

class CoercionPDtoKM(HyperbolicModelCoercion):
    """
    Coercion from the PD to KM model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::
        """
        return (2*real(x)/(Integer(1) + real(x)**2 +imag(x)**2),
                2*imag(x)/(Integer(1) + real(x)**2 + imag(x)**2))

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return SL2R_to_SO21( matrix(2,[1,I,I,1]) * x *
                             matrix(2,[1,-I,-I,1])/Integer(2) )

class CoercionPDtoHM(HyperbolicModelCoercion):
    """
    Coercion from the PD to HM model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::
        """
        return vector((
            2*real(x)/(1 - real(x)**2 - imag(x)**2),
            2*imag(x)/(1 - real(x)**2 - imag(x)**2),
            (real(x)**2 + imag(x)**2 + 1)/(1 - real(x)**2 - imag(x)**2)
            ))

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return SL2R_to_SO21( matrix(2,[1,I,I,1]) * x *
                             matrix(2,[1,-I,-I,1])/Integer(2) )

###########
# From KM #
###########

class CoercionKMtoUHP(HyperbolicModelCoercion):
    """
    Coercion from the KM to UHP model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HyperbolicPlane.KM.point_to_model((0, 0), 'UHP')
            I
        """
        if tuple(x) == (0, 1):
            return infinity
        return ( -x[0]/(x[1] - 1)
                 + I*(-(sqrt(-x[0]**2 -x[1]**2 + 1) - x[0]**2 - x[1]**2 + 1)
                      / ((x[1] - 1)*sqrt(-x[0]**2 - x[1]**2 + 1) + x[1] - 1)) )

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return SO21_to_SL2R(x)

class CoercionKMtoPD(HyperbolicModelCoercion):
    """
    Coercion from the KM to PD model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::
        """
        return ( x[0]/(1 + (1 - x[0]**2 - x[1]**2).sqrt())
                 + I*x[1]/(1 + (1 - x[0]**2 - x[1]**2).sqrt()) )

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return (matrix(2,[1,-I,-I,1]) * SO21_to_SL2R(x) *
                matrix(2,[1,I,I,1])/Integer(2))

class CoercionKMtoHM(HyperbolicModelCoercion):
    """
    Coercion from the KM to HM model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HyperbolicPlane.KM.point_to_model((0, 0), 'HM')
            (0, 0, 1)
        """
        return (vector((2*x[0], 2*x[1], 1 + x[0]**2 + x[1]**2))
                / (1 - x[0]**2 - x[1]**2))

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return x

###########
# From HM #
###########

class CoercionHMtoUHP(HyperbolicModelCoercion):
    """
    Coercion from the HM to UHP model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HyperbolicPlane.HM.point_to_model(vector((0,0,1)), 'UHP')
            I
        """
        return -((x[0]*x[2] + x[0]) + I*(x[2] + 1)) / ((x[1] - 1)*x[2]
                                        - x[0]**2 - x[1]**2 + x[1] - 1)

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return SO21_to_SL2R(x)

class CoercionHMtoPD(HyperbolicModelCoercion):
    """
    Coercion from the HM to PD model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::
        """
        return x[0]/(1 + x[2]) + I*(x[1]/(1 + x[2]))

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return (matrix(2,[1,-I,-I,1]) * SO21_to_SL2R(x) *
                matrix(2,[1,I,I,1])/Integer(2))

class CoercionHMtoKM(HyperbolicModelCoercion):
    """
    Coercion from the HM to KM model.
    """
    def image_point(self, x):
        """
        Return the image of the coordinates of the hyperbolic point ``x``
        under ``self``.

        EXAMPLES::

            sage: HyperbolicPlane.HM.point_to_model(vector((0,0,1)), 'KM')
            (0, 0)
        """
        return (x[0]/(1 + x[2]), x[1]/(1 + x[2]))

    def image_isometry(self, x):
        """
        Return the image of the matrix of the hyperbolic isometry ``x``
        under ``self``.

        EXAMPLES::
        """
        return x

#####################################################################
## Helper functions

def SL2R_to_SO21(A):
    r"""
    Given a matrix in `SL(2, \RR)` return its irreducible representation in
    `O(2,1)`.

    Note that this is not the only homomorphism, but it is the only one
    that works in the context of the implemented 2D hyperbolic geometry
    models.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_coercion import SL2R_to_SO21
        sage: A = SL2R_to_SO21(identity_matrix(2))
        sage: J =  matrix([[1,0,0],[0,1,0],[0,0,-1]]) #Lorentzian Gram matrix
        sage: norm(A.transpose()*J*A - J) < 10**-4
        True
    """
    a,b,c,d = (A/A.det().sqrt()).list()
    B = matrix(3, [a*d + b*c, a*c - b*d, a*c + b*d, a*b - c*d,
                   Integer(1)/Integer(2)*a**2 - Integer(1)/Integer(2)*b**2 -
                   Integer(1)/Integer(2)*c**2 + Integer(1)/Integer(2)*d**2,
                   Integer(1)/Integer(2)*a**2 + Integer(1)/Integer(2)*b**2 -
                   Integer(1)/Integer(2)*c**2 - Integer(1)/Integer(2)*d**2,
                   a*b + c*d, Integer(1)/Integer(2)*a**2 -
                   Integer(1)/Integer(2)*b**2 + Integer(1)/Integer(2)*c**2 -
                   Integer(1)/Integer(2)*d**2, Integer(1)/Integer(2)*a**2 +
                   Integer(1)/Integer(2)*b**2 + Integer(1)/Integer(2)*c**2 +
                   Integer(1)/Integer(2)*d**2])
    B = B.apply_map(attrcall('real')) # Kill ~0 imaginary parts
    if A.det() > 0:
        return B
    else:
        # Orientation-reversing isometries swap the nappes of
        #  the lightcone.  This fixes that issue.
        return -B

def SO21_to_SL2R(M):
    r"""
    A homomorphism from `SO(2, 1)` to `SL(2, \RR)`.

    Note that this is not the only homomorphism, but it is the only one
    that works in the context of the implemented 2D hyperbolic geometry
    models.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_coercion import SO21_to_SL2R
        sage: (SO21_to_SL2R(identity_matrix(3)) - identity_matrix(2)).norm() < 10**-4
        True
    """
    ####################################################################
    # SL(2,R) is the double cover of SO (2,1)^+, so we need to choose  #
    # a lift.  I have formulas for the absolute values of each entry   #
    # a,b ,c,d of the lift matrix(2,[a,b,c,d]), but we need to choose  #
    # one entry to be positive.  I choose d for no particular reason,  #
    # unless d = 0, then we choose c > 0.  The basic strategy for this #
    # function is to find the linear map induced by the SO(2,1)        #
    # element on the Lie algebra sl(2, R).  This corresponds to the    #
    # Adjoint action by a matrix A or -A in SL(2,R).  To find which    #
    # matrix let X,Y,Z be a basis for sl(2,R) and look at the images   #
    # of X,Y,Z as well as the second and third standard basis vectors  #
    # for 2x2 matrices (these are traceless, so are in the Lie         #
    # algebra).  These corresponds to AXA^-1 etc and give formulas     #
    # for the entries of A.                                            #
    ####################################################################
    (m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9) = M.list()
    d = sqrt(Integer(1)/Integer(2)*m_5 - Integer(1)/Integer(2)*m_6 -
             Integer(1)/Integer(2)*m_8 + Integer(1)/Integer(2)*m_9)
    if M.det() > 0: #EPSILON?
        det_sign = 1
    elif M.det() < 0: #EPSILON?
        det_sign = -1
    if d > 0: #EPSILON?
        c = (-Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/d
        b = (-Integer(1)/Integer(2)*m_2 + Integer(1)/Integer(2)*m_3)/d
        ad = det_sign*1 + b*c # ad - bc = pm 1
        a = ad/d
    else: # d is 0, so we make c > 0
        c = sqrt(-Integer(1)/Integer(2)*m_5 - Integer(1)/Integer(2)*m_6 +
                  Integer(1)/Integer(2)*m_8 + Integer(1)/Integer(2)*m_9)
        d = (-Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/c
            #d = 0, so ad - bc = -bc = pm 1.
        b = - (det_sign*1)/c
        a = (Integer(1)/Integer(2)*m_4 + Integer(1)/Integer(2)*m_7)/b
    A = matrix(2,[a,b,c,d])
    return A

