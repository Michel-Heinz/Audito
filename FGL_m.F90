! Copyright (C) 2023 Michel Heinz

#  define _PASTE(X) X
#  define _CAT2(X,Y) _PASTE(X)Y
#  define _CAT3(X,Y,Z) _CAT2(_CAT2(X,Y),Z)


module FGL_m
    use, intrinsic :: ISO_FORTRAN_ENV, only: REAL64, REAL32, REAL128, INT8, INT16, INT32, INT64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_NEGATIVE_INF, IEEE_QUIET_NAN

    implicit none

    private
    public :: FGL_t, operator(+), operator(-), operator(*), operator(**), operator(/), ABS, SQRT, EXP, LOG, COS, SIN
    public :: TAN, INT, NINT, FRACTION, REAL, FLOOR
    public :: i1, i2, i4, i8, r4, r8, r16


    integer(INT32), parameter :: i1 = INT8
    integer(INT32), parameter :: i2 = INT16
    integer(INT32), parameter :: i4 = INT32
    integer(INT32), parameter :: i8 = INT64
    integer(INT32), parameter :: r4 = REAL32
    integer(INT32), parameter :: r8 = REAL64
    integer(INT32), parameter :: r16 = REAL128

    type FGL_t
        real(r8) :: f
        real(r8), allocatable :: g(:)
        real(r8) :: l
    contains
        procedure :: print => FGL_Print
        procedure, private :: FGL_Neg
        generic :: operator(-) => FGL_Neg
        procedure, private :: FGL_Add
        generic :: operator(+) => FGL_Add
        procedure, private :: FGL_Sub
        generic :: operator(-) => FGL_Sub
        procedure, private :: FGL_MulSame
        generic :: operator(*) => FGL_MulSame
        procedure, private :: FGL_DivSame
        generic :: operator(/) => FGL_DivSame
        procedure, private :: FGL_Pow
        generic :: operator(**) => FGL_Pow
        procedure, private :: FGL_Equal
        generic :: operator(==) => FGL_Equal
        procedure, private :: FGL_NotEqual
        generic :: operator(.ne.) => FGL_NotEqual
        procedure, private :: FGL_LessThan
        generic :: operator(<) => FGL_LessThan
        procedure, private :: FGL_GreaterThan
        generic :: operator(>) => FGL_GreaterThan
        procedure, private :: FGL_LessThanOrEqual
        generic :: operator(<=) => FGL_LessThanOrEqual
        procedure, private :: FGL_GreaterThanOrEqual
        generic :: operator(>=) => FGL_GreaterThanOrEqual
        procedure, private :: FGL_Assign
        generic :: assignment(=) => FGL_Assign
        procedure, private :: FGL_Mul
        generic :: operator(.mul.) => FGL_Mul
        procedure, private :: FGL_Div
        generic :: operator(.div.) => FGL_Div
    end type FGL_t

    interface operator(+)
        module procedure :: FGL_AddScalarObj_r8, FGL_AddObjScalar_r8, FGL_AddScalarObj_r4, FGL_AddObjScalar_r4
        module procedure :: FGL_AddScalarObj_r16, FGL_AddObjScalar_r16, FGL_AddScalarObj_i1, FGL_AddObjScalar_i1
        module procedure :: FGL_AddScalarObj_i2, FGL_AddObjScalar_i2, FGL_AddScalarObj_i4, FGL_AddObjScalar_i4
        module procedure :: FGL_AddScalarObj_i8, FGL_AddObjScalar_i8
    end interface

    interface operator(-)
        module procedure :: FGL_SubScalarObj_r8, FGL_SubObjScalar_r8, FGL_SubScalarObj_r4, FGL_SubObjScalar_r4
        module procedure :: FGL_SubScalarObj_r16, FGL_SubObjScalar_r16, FGL_SubScalarObj_i1, FGL_SubObjScalar_i1
        module procedure :: FGL_SubScalarObj_i2, FGL_SubObjScalar_i2, FGL_SubScalarObj_i4, FGL_SubObjScalar_i4
        module procedure :: FGL_SubScalarObj_i8, FGL_SubObjScalar_i8
    end interface

    interface operator(*)
        module procedure :: FGL_MulScalarObj_r8, FGL_MulObjScalar_r8, FGL_MulScalarObj_r4, FGL_MulObjScalar_r4
        module procedure :: FGL_MulScalarObj_r16, FGL_MulObjScalar_r16, FGL_MulScalarObj_i1, FGL_MulObjScalar_i1
        module procedure :: FGL_MulScalarObj_i2, FGL_MulObjScalar_i2, FGL_MulScalarObj_i4, FGL_MulObjScalar_i4
        module procedure :: FGL_MulScalarObj_i8, FGL_MulObjScalar_i8
    end interface

    interface operator(/)
        module procedure :: FGL_DivScalarObj_r8, FGL_DivObjScalar_r8, FGL_DivScalarObj_r4, FGL_DivObjScalar_r4
        module procedure :: FGL_DivScalarObj_r16, FGL_DivObjScalar_r16, FGL_DivScalarObj_i1, FGL_DivObjScalar_i1
        module procedure :: FGL_DivScalarObj_i2, FGL_DivObjScalar_i2, FGL_DivScalarObj_i4, FGL_DivObjScalar_i4
        module procedure :: FGL_DivScalarObj_i8, FGL_DivObjScalar_i8
    end interface

    interface operator(**)
        module procedure :: FGL_PowObjScalar_r4, FGL_PowObjScalar_r8, FGL_PowObjScalar_r16, FGL_PowObjScalar_i1
        module procedure :: FGL_PowObjScalar_i2, FGL_PowObjScalar_i4, FGL_PowObjScalar_i8
        module procedure :: FGL_PowScalarObj_r4, FGL_PowScalarObj_r8, FGL_PowScalarObj_r16, FGL_PowScalarObj_i1
        module procedure :: FGL_PowScalarObj_i2, FGL_PowScalarObj_i4, FGL_PowScalarObj_i8
    end interface

    interface ABS
        module procedure :: FGL_Abs_real
    end interface

    interface SQRT
        module procedure :: FGL_sqrt_real
    end interface

    interface EXP
        module procedure :: FGL_exp
    end interface

    interface LOG
        module procedure :: FGL_log_real
    end interface

    interface COS
        module procedure :: FGL_cos_real
    end interface

    interface SIN
        module procedure :: FGL_sin_real
    end interface

    interface TAN
        module procedure :: FGL_tan_real
    end interface

    interface INT
        module procedure :: FGL_int_real
    end interface

    interface NINT
        module procedure :: FGL_nint_real
    end interface

    interface FLOOR
        module procedure :: FGL_floor_real
    end interface

    interface FRACTION
        module procedure :: FGL_fraction_real
    end interface

    interface REAL
        module procedure :: FGL_real_real
    end interface

contains

    subroutine FGL_Print(this, iu)
        class(FGL_t), intent(in)   :: this
        integer, optional, intent(in) :: iu
        integer :: n, i
        character(len=10) :: n_str

        n = SIZE(this%g)
        write(n_str, '(i10)') n

        if (PRESENT(iu)) then
            write(iu,'(a)') 'value:'
            write(iu,'(e13.5)') this%f
            write(iu,'(a)') 'gradient:'
            write(iu,'('//TRIM(n_str)//'e11.3)') this%g
            write(iu,'(a)') 'Laplacian:'
            write(iu,'('//TRIM(n_str)//'e13.5)') this%l
        else
            print'(a)', 'value:'
            print'(e13.5)', this%f
            print'(a)', 'gradient:'
            print'('//TRIM(n_str)//'e11.3)', this%g
            print'(a)', 'Laplacian:'
            print'('//TRIM(n_str)//'e13.5)', this%l
        end if
    end subroutine FGL_Print


    elemental type(FGL_t) function FGL_DivSame(this, that) result(new)
        class(FGL_t), intent(in) :: this, that
!TODO needs testing!

!        new%f = this%f / that%f
!        new%g = ((this%g * that%f) - (this%f * that%g)) / that%f**2
!        new%l = (this%l*that%f**2 - 2*DOT_PRODUCT(this%g, that%g) + 2*DOT_PRODUCT(that%g, that%g)*this%f&
!        - that%l*this%f*that%f) / that%f**3

        new = this%FGL_MulSame(FGL_PowObjScalar_i1(that, -1_i1))
    end function FGL_DivSame


    elemental type(FGL_t) function FGL_Abs_real(this) result(new)
        class(FGL_t), intent(in) :: this

        if (this%f < 0) then
            new%f = -this%f
            new%g = -this%g
            new%l = -this%l
        else
            new%f = this%f
            new%g = this%g
            new%l = this%l
        end if
    end function FGL_Abs_real


    elemental type(FGL_t) function FGL_sqrt_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new%f = SQRT(this%f)
        new%g = this%g / (2 * new%f)
        new%l = (this%l * this%f - DOT_PRODUCT(this%g, this%g)) / (4*this%f*new%f)

    end function FGL_sqrt_real


    elemental type(FGL_t) function FGL_exp(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new

        new%f = EXP(this%f)
        new%g = new%f * this%g
        new%l = new%f * (DOT_PRODUCT(this%g, this%g) + this%l)

    end function FGL_exp


    elemental type(FGL_t) function FGL_log_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new
        real(r8)    :: dummy1(SIZE(this%g))

        if (this%f > 0) then
            new%f = LOG(this%f)
            new%g = this%g / this%f
            new%l = -1 / (this%f ** 2) * DOT_PRODUCT(this%g, this%g) + 1 / this%f * this%l
        else
            dummy1 = IEEE_VALUE(1._r8, IEEE_QUIET_NAN)
            new = FGL_t(IEEE_VALUE(1._r8, IEEE_NEGATIVE_INF), dummy1, IEEE_VALUE(1._r8, IEEE_QUIET_NAN))
        end if

    end function FGL_log_real


    elemental type(FGL_t) function FGL_cos_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new

        new%f = COS(this%f)
        new%g = - SIN(this%f) * this%g
        new%l = -new%f * DOT_PRODUCT(this%g, this%g) - SIN(this%f) * this%l

    end function FGL_cos_real


    elemental type(FGL_t) function FGL_sin_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new

        new%f = SIN(this%f)
        new%g = this%g * COS(this%f)
        new%l =  -new%f * DOT_PRODUCT(this%g, this%g) + COS(this%f) * this%l

    end function FGL_sin_real


    elemental type(FGL_t) function FGL_tan_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new, sec

        new%f = TAN(this%f)
        sec = 1 / cos(this%f)**2
        new%g = this%g * sec
        new%l = sec * (2 * new%f * DOT_PRODUCT(this%g, this%g) + this%l)

    end function FGL_tan_real


    elemental type(FGL_t) function FGL_int_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new%f = INT(this%f)
        new%g = this%g
        new%l = this%l

    end function FGL_int_real


    elemental type(FGL_t) function FGL_nint_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new%f = NINT(this%f)
        new%g = this%g
        new%l = this%l

    end function FGL_nint_real


    elemental type(FGL_t) function FGL_floor_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new%f = FLOOR(this%f)
        new%g = this%g
        new%l = this%l

    end function FGL_floor_real


    elemental type(FGL_t) function FGL_fraction_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new%f = FRACTION(this%f)
        new%g = this%g
        new%l = this%l

    end function FGL_fraction_real


    elemental type(FGL_t) function FGL_real_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new%f = REAL(this%f)
        new%g = this%g
        new%l = this%l

    end function FGL_real_real


    elemental type(FGL_t) function FGL_Div(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new = this%FGL_Mul(FGL_PowObjScalar_i1(that, -1_i1))
    end function FGL_Div


    elemental subroutine FGL_Assign(this, that)
        class(FGL_t), intent(inout) :: this
        class(FGL_t), intent(in)   :: that

        this%f = that%f
        this%g = that%g
        this%l = that%l
    end subroutine FGL_Assign


    elemental logical function FGL_Equal(this, that) result(equal)
        class(FGL_t), intent(in) :: this, that

        equal = .false.
        if (this%f == that%f .and. ALL(this%g == that%g) .and. this%l == that%l) &
                equal = .true.
    end function FGL_Equal


    elemental logical function FGL_NotEqual(this, that) result(notEqual)
        class(FGL_t), intent(in) :: this, that

        notEqual = .false.
        if (this%f /= that%f .and. ALL(this%g /= that%g) .and. this%l /= that%l) &
                notEqual = .true.
    end function FGL_NotEqual


    elemental logical function FGL_LessThan(this, that) result(LessThan)
        class(FGL_t), intent(in) :: this, that

        LessThan = .false.
        if (this%f < that%f) LessThan = .true.
    end function FGL_LessThan


    elemental logical function FGL_GreaterThan(this, that) result(GreaterThan)
        class(FGL_t), intent(in) :: this, that

        GreaterThan = .false.
        if (this%f > that%f) GreaterThan = .true.
    end function FGL_GreaterThan


    elemental logical function FGL_LessThanOrEqual(this, that) result(LessThanOrEqual)
        class(FGL_t), intent(in) :: this, that

        LessThanOrEqual = .false.
        if (this%f <= that%f) LessThanOrEqual = .true.
    end function FGL_LessThanOrEqual


    elemental logical function FGL_GreaterThanOrEqual(this, that) result(GreaterThanOrEqual)
        class(FGL_t), intent(in) :: this, that

        GreaterThanOrEqual = .false.
        if (this%f >= that%f) GreaterThanOrEqual = .true.
    end function FGL_GreaterThanOrEqual


    elemental type(FGL_t) function FGL_Neg(this) result(new)
        class(FGL_t), intent(in) :: this

        new%f = -this%f
        new%g = -this%g
        new%l = -this%l

    end function FGL_Neg


    elemental type(FGL_t) function FGL_Add(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new%f = this%f + that%f
        new%g = this%g + that%g
        new%l = this%l + that%l

    end function FGL_Add


    elemental type(FGL_t) function FGL_Sub(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new%f = this%f - that%f
        new%g = this%g - that%g
        new%l = this%l - that%l

    end function FGL_Sub


    elemental type(FGL_t) function FGL_Mul(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new = FGL_t(this%f * that%f, [that%f * this%g, this%f * that%g], &
        2 * DOT_PRODUCT(this%g, that%g) + this%f * that%l + that%f * this%l)
    end function FGL_Mul


    elemental type(FGL_t) function FGL_MulSame(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new%f = this%f * that%f
        new%g = this%f * that%g + that%f * this%g
        new%l = this%l * that%f + that%l * this%f + 2 * DOT_PRODUCT(this%g, that%g)

    end function FGL_MulSame
    
    
    elemental type(FGL_t) function FGL_Pow(this, that) result(new)
        class(FGL_t), intent(in) :: this, that
        real(r8) :: dotxy, dotxx
        real(r8) :: dotyy, logval, frac

        dotxy = DOT_PRODUCT(this%g, that%g)
        dotxx = DOT_PRODUCT(this%g, this%g)
        dotyy = DOT_PRODUCT(that%g, that%g)

        logval = LOG(this%f)
        frac = (that%f/this%f)
        new%f = this%f ** that%f
        new%g = new%f * (that%g * LOG(this%f) + (that%f / this%f) * this%g)
        new%l = new%f * (dotyy * logval**2 + dotxx * frac**2 + 2 * frac * logval * dotxy +&
                that%l * logval + 2 * dotxy / this%f - dotxx / this%f**2 + this%l * frac)

    end function FGL_Pow

#define _FGL true
    ! single precision
#define _F_TYPE FGL_t
#define _F FGL
#define _KIND r4
#define _TYPE real(r4)
#include "F_functions.inc"



    ! double precision
#define _F_TYPE FGL_t
#define _F FGL
#define _KIND r8
#define _TYPE real(r8)
#include "F_functions.inc"



    ! quad precision
#define _F_TYPE FGL_t
#define _F FGL
#define _KIND r16
#define _TYPE real(r16)
#include "F_functions.inc"



    ! 8bit int
#define _INT true
#define _F_TYPE FGL_t
#define _F FGL
#define _KIND i1
#define _TYPE integer(i1)
#include "F_functions.inc"



    ! 16bit int
#define _INT true
#define _F_TYPE FGL_t
#define _F FGL
#define _KIND i2
#define _TYPE integer(i2)
#include "F_functions.inc"



    ! 32bit int
#define _INT true
#define _F_TYPE FGL_t
#define _F FGL
#define _KIND i4
#define _TYPE integer(i4)
#include "F_functions.inc"



    ! 64bit int
#define _INT true
#define _F_TYPE FGL_t
#define _F FGL
#define _KIND i8
#define _TYPE integer(i8)
#include "F_functions.inc"


#undef _FGL
end module FGL_m