! Copyright (C) 2023 Michel Heinz

#  define _PASTE(X) X
#  define _CAT2(X,Y) _PASTE(X)Y
#  define _CAT3(X,Y,Z) _CAT2(_CAT2(X,Y),Z)


module FGH_m
    use, intrinsic :: ISO_FORTRAN_ENV, only: REAL64, REAL32, REAL128, INT8, INT16, INT32, INT64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_NEGATIVE_INF, IEEE_QUIET_NAN

    implicit none

    private
    public :: FGH_t, operator(+), operator(-), operator(*), operator(**), operator(/), ABS, SQRT, EXP, LOG, COS, SIN
    public :: TAN, INT, NINT, FRACTION, REAL, FLOOR
    public :: i1, i2, i4, i8, r4, r8, r16


    integer(INT32), parameter :: i1 = INT8
    integer(INT32), parameter :: i2 = INT16
    integer(INT32), parameter :: i4 = INT32
    integer(INT32), parameter :: i8 = INT64
    integer(INT32), parameter :: r4 = REAL32
    integer(INT32), parameter :: r8 = REAL64
    integer(INT32), parameter :: r16 = REAL128

    type FGH_t
        real(r8) :: f
        real(r8), allocatable :: g(:)
        real(r8), allocatable :: h(:, :)
    contains
        procedure :: print => FGH_Print
        procedure, private :: FGH_Neg
        generic :: operator(-) => FGH_Neg
        procedure, private :: FGH_Add
        generic :: operator(+) => FGH_Add
        procedure, private :: FGH_Sub
        generic :: operator(-) => FGH_Sub
        procedure, private :: FGH_MulSame
        generic :: operator(*) => FGH_MulSame
        procedure, private :: FGH_DivSame
        generic :: operator(/) => FGH_DivSame
        procedure, private :: FGH_Pow
        generic :: operator(**) => FGH_Pow
        procedure, private :: FGH_Equal
        generic :: operator(==) => FGH_Equal
        procedure, private :: FGH_NotEqual
        generic :: operator(.ne.) => FGH_NotEqual
        procedure, private :: FGH_LessThan
        generic :: operator(<) => FGH_LessThan
        procedure, private :: FGH_GreaterThan
        generic :: operator(>) => FGH_GreaterThan
        procedure, private :: FGH_LessThanOrEqual
        generic :: operator(<=) => FGH_LessThanOrEqual
        procedure, private :: FGH_GreaterThanOrEqual
        generic :: operator(>=) => FGH_GreaterThanOrEqual
        procedure, private :: FGH_Assign
        generic :: assignment(=) => FGH_Assign
        procedure, private :: FGH_Mul
        generic :: operator(.mul.) => FGH_Mul
        procedure, private :: FGH_Div
        generic :: operator(.div.) => FGH_Div
    end type FGH_t

    interface operator(+)
        module procedure :: FGH_AddScalarObj_r8, FGH_AddObjScalar_r8, FGH_AddScalarObj_r4, FGH_AddObjScalar_r4
        module procedure :: FGH_AddScalarObj_r16, FGH_AddObjScalar_r16, FGH_AddScalarObj_i1, FGH_AddObjScalar_i1
        module procedure :: FGH_AddScalarObj_i2, FGH_AddObjScalar_i2, FGH_AddScalarObj_i4, FGH_AddObjScalar_i4
        module procedure :: FGH_AddScalarObj_i8, FGH_AddObjScalar_i8
    end interface

    interface operator(-)
        module procedure :: FGH_SubScalarObj_r8, FGH_SubObjScalar_r8, FGH_SubScalarObj_r4, FGH_SubObjScalar_r4
        module procedure :: FGH_SubScalarObj_r16, FGH_SubObjScalar_r16, FGH_SubScalarObj_i1, FGH_SubObjScalar_i1
        module procedure :: FGH_SubScalarObj_i2, FGH_SubObjScalar_i2, FGH_SubScalarObj_i4, FGH_SubObjScalar_i4
        module procedure :: FGH_SubScalarObj_i8, FGH_SubObjScalar_i8
    end interface

    interface operator(*)
        module procedure :: FGH_MulScalarObj_r8, FGH_MulObjScalar_r8, FGH_MulScalarObj_r4, FGH_MulObjScalar_r4
        module procedure :: FGH_MulScalarObj_r16, FGH_MulObjScalar_r16, FGH_MulScalarObj_i1, FGH_MulObjScalar_i1
        module procedure :: FGH_MulScalarObj_i2, FGH_MulObjScalar_i2, FGH_MulScalarObj_i4, FGH_MulObjScalar_i4
        module procedure :: FGH_MulScalarObj_i8, FGH_MulObjScalar_i8
    end interface

    interface operator(/)
        module procedure :: FGH_DivScalarObj_r8, FGH_DivObjScalar_r8, FGH_DivScalarObj_r4, FGH_DivObjScalar_r4
        module procedure :: FGH_DivScalarObj_r16, FGH_DivObjScalar_r16, FGH_DivScalarObj_i1, FGH_DivObjScalar_i1
        module procedure :: FGH_DivScalarObj_i2, FGH_DivObjScalar_i2, FGH_DivScalarObj_i4, FGH_DivObjScalar_i4
        module procedure :: FGH_DivScalarObj_i8, FGH_DivObjScalar_i8
    end interface

    interface operator(**)
        module procedure :: FGH_PowObjScalar_r4, FGH_PowObjScalar_r8, FGH_PowObjScalar_r16, FGH_PowObjScalar_i1
        module procedure :: FGH_PowObjScalar_i2, FGH_PowObjScalar_i4, FGH_PowObjScalar_i8
        module procedure :: FGH_PowScalarObj_r4, FGH_PowScalarObj_r8, FGH_PowScalarObj_r16, FGH_PowScalarObj_i1
        module procedure :: FGH_PowScalarObj_i2, FGH_PowScalarObj_i4, FGH_PowScalarObj_i8
    end interface

    interface ABS
        module procedure :: FGH_Abs_real
    end interface

    interface SQRT
        module procedure :: FGH_sqrt_real
    end interface

    interface EXP
        module procedure :: FGH_exp_real
    end interface

    interface LOG
        module procedure :: FGH_log_real
    end interface

    interface COS
        module procedure :: FGH_cos_real
    end interface

    interface SIN
        module procedure :: FGH_sin_real
    end interface

    interface TAN
        module procedure :: FGH_tan_real
    end interface

    interface INT
        module procedure :: FGH_int_real
    end interface

    interface NINT
        module procedure :: FGH_nint_real
    end interface

    interface FLOOR
        module procedure :: FGH_floor_real
    end interface

    interface FRACTION
        module procedure :: FGH_fraction_real
    end interface

    interface REAL
        module procedure :: FGH_real_real
    end interface

    interface FGH_t
        module procedure :: constructor1
    end interface FGH_t

contains

    pure type(FGH_t) function constructor1(f, g, h) result(new)
        real(r8), intent(in) :: f, g(:), h(:,:)

        new%f = f
        if (SIZE(h(1,:)) /= SIZE(h(:,1))) ERROR STOP 'FGH: Hessian must be a square matrix!'
        if (SIZE(g) /= SIZE(h(:,1)) .or. SIZE(g) /= SIZE(h(1,:))) ERROR STOP 'FGH: gradient and Hessian have different&
                sizes!'
        new%g = g
        new%h = h
    end function constructor1


    subroutine FGH_Print(this, iu)
        class(FGH_t), intent(in)   :: this
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
            write(iu,'(a)') 'Hessian:'
            do i = 1, n
                write(iu,'('//TRIM(n_str)//'e11.3)') this%h(i,:)
            end do
        else
            print'(a)', 'value:'
            print'(e13.5)', this%f
            print'(a)', 'gradient:'
            print'('//TRIM(n_str)//'e11.3)', this%g
            print'(a)', 'Hessian:'
            do i = 1, n
                print'('//TRIM(n_str)//'e11.3)', this%h(i,:)
            end do
        end if
    end subroutine FGH_Print


    elemental type(FGH_t) function FGH_DivSame(this, that) result(new)
        class(FGH_t), intent(in) :: this, that
        !TODO needs testing!
        real(r8)    :: outxy(SIZE(this%g), SIZE(this%g)), outyy(SIZE(this%g), SIZE(this%g))
        integer     :: a


        a = SIZE(this%g)
        outxy = SPREAD(this%g, 2, a) * SPREAD(that%g, 1, a)
        outyy = SPREAD(that%g, 2, a) * SPREAD(that%g, 1, a)
        new%f = this%f / that%f
        new%g = ((this%g * that%f) - (this%f * that%g)) / that%f**2
        new%h = (this%h*that%f**2 - TRANSPOSE(outxy)*that%f - outxy*that%f + 2*outyy*this%f - that%h*this%f*that%f) / &
                that%f**3

    end function FGH_DivSame


    elemental type(FGH_t) function FGH_Abs_real(this) result(new)
        class(FGH_t), intent(in) :: this

        if (this%f < 0) then
            new%f = -this%f
            new%g = -this%g
            new%h = -this%h
        else
            new%f = this%f
            new%g = this%g
            new%h = this%h
        end if
    end function FGH_Abs_real


    elemental type(FGH_t) function FGH_sqrt_real(this) result(new)
        class(FGH_t), intent(in) :: this
        real(r8)    :: out(SIZE(this%g), SIZE(this%g))
        integer     :: a

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        new%f = SQRT(this%f)
        new%g = this%g / (2 * new%f)
        new%h = (this%h * this%f - out) / (4*this%f*new%f)

    end function FGH_sqrt_real


    elemental type(FGH_t) function FGH_exp_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g))

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        new%f = EXP(this%f)
        new%g = new%f * this%g
        new%h = new%f * (out + this%h)

    end function FGH_exp_real


    elemental type(FGH_t) function FGH_log_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g))
        real(r8)    :: dummy1(SIZE(this%g)), dummy2(SIZE(this%g), SIZE(this%g))

        a = SIZE(this%g)
        if (this%f > 0) then
            out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
            new%f = LOG(this%f)
            new%g = this%g / this%f
            new%h = (this%f * this%h - out) / this%f**2
        else
            dummy1 = IEEE_VALUE(1._r8, IEEE_QUIET_NAN)
            dummy2 = IEEE_VALUE(1._r8, IEEE_QUIET_NAN)
            new = FGH_t(IEEE_VALUE(1._r8, IEEE_NEGATIVE_INF), dummy1, dummy2)
        end if

    end function FGH_log_real


    elemental type(FGH_t) function FGH_cos_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g))

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        new%f = COS(this%f)
        new%g = - SIN(this%f) * this%g
        new%h = -new%f * out - SIN(this%f) * this%h

    end function FGH_cos_real


    elemental type(FGH_t) function FGH_sin_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g))

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        new%f = SIN(this%f)
        new%g = this%g * COS(this%f)
        new%h =  -new%f * out + COS(this%f) * this%h

    end function FGH_sin_real


    elemental type(FGH_t) function FGH_tan_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g)), sec

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        new%f = TAN(this%f)
        sec = 1 / cos(this%f)**2
        new%g = this%g * sec
        new%h = sec * (2 * new%f * out + this%h)

    end function FGH_tan_real


    elemental type(FGH_t) function FGH_int_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new%f = INT(this%f)
        new%g = this%g
        new%h = this%h

    end function FGH_int_real


    elemental type(FGH_t) function FGH_nint_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new%f = NINT(this%f)
        new%g = this%g
        new%h = this%h

    end function FGH_nint_real


    elemental type(FGH_t) function FGH_floor_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new%f = FLOOR(this%f)
        new%g = this%g
        new%h = this%h

    end function FGH_floor_real


    elemental type(FGH_t) function FGH_fraction_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new%f = FRACTION(this%f)
        new%g = this%g
        new%h = this%h

    end function FGH_fraction_real


    elemental type(FGH_t) function FGH_real_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new%f = REAL(this%f)
        new%g = this%g
        new%h = this%h

    end function FGH_real_real


    elemental type(FGH_t) function FGH_Div(this, that) result(new)
        class(FGH_t), intent(in) :: this, that

        new = this%FGH_Mul(FGH_PowObjScalar_i1(that, -1_i1))

    end function FGH_Div


    elemental subroutine FGH_Assign(this, that)
        class(FGH_t), intent(inout) :: this
        class(FGH_t), intent(in)   :: that

        this%f = that%f
        this%g = that%g
        this%h = that%h

    end subroutine FGH_Assign


    elemental logical function FGH_Equal(this, that) result(equal)
        class(FGH_t), intent(in) :: this, that

        equal = .false.
        if (this%f == that%f .and. ALL(this%g == that%g) .and. ALL(this%h == that%h)) &
            equal = .true.
    end function FGH_Equal


    elemental logical function FGH_NotEqual(this, that) result(notEqual)
        class(FGH_t), intent(in) :: this, that

        notEqual = .false.
        if (this%f /= that%f .and. ALL(this%g /= that%g) .and. ALL(this%h /= that%h)) &
                notEqual = .true.
    end function FGH_NotEqual


    elemental logical function FGH_LessThan(this, that) result(LessThan)
        class(FGH_t), intent(in) :: this, that

        LessThan = .false.
        if (this%f < that%f) LessThan = .true.
    end function FGH_LessThan


    elemental logical function FGH_GreaterThan(this, that) result(GreaterThan)
        class(FGH_t), intent(in) :: this, that

        GreaterThan = .false.
        if (this%f > that%f) GreaterThan = .true.
    end function FGH_GreaterThan


    elemental logical function FGH_LessThanOrEqual(this, that) result(LessThanOrEqual)
        class(FGH_t), intent(in) :: this, that

        LessThanOrEqual = .false.
        if (this%f <= that%f) LessThanOrEqual = .true.
    end function FGH_LessThanOrEqual


    elemental logical function FGH_GreaterThanOrEqual(this, that) result(GreaterThanOrEqual)
        class(FGH_t), intent(in) :: this, that

        GreaterThanOrEqual = .false.
        if (this%f >= that%f) GreaterThanOrEqual = .true.
    end function FGH_GreaterThanOrEqual


    elemental type(FGH_t) function FGH_Neg(this) result(new)
        class(FGH_t), intent(in) :: this

        new%f = -this%f
        new%g = -this%g
        new%h = -this%h
    end function FGH_Neg


    elemental type(FGH_t) function FGH_Add(this, that) result(new)
        class(FGH_t), intent(in) :: this, that

        new%f = this%f + that%f
        new%g = this%g + that%g
        new%h = this%h + that%h

    end function FGH_Add


    elemental type(FGH_t) function FGH_Sub(this, that) result(new)
        class(FGH_t), intent(in) :: this, that

        new%f = this%f - that%f
        new%g = this%g - that%g
        new%h = this%h - that%h

    end function FGH_Sub


    elemental type(FGH_t) function FGH_Mul(this, that) result(new)
        class(FGH_t), intent(in) :: this, that

        integer :: a, b
        real(r8) :: out(SIZE(this%g), SIZE(that%g)), &
                h(SIZE(this%g) + SIZE(that%g), SIZE(this%g) + SIZE(that%g))

        a = SIZE(this%g)
        b = SIZE(that%g)
        out = SPREAD(this%g, 2, b) * SPREAD(that%g, 1, a)
        h(:a, :a) = this%h * that%f
        h(a+1:, a+1:) = that%h * this%f
        h(:a, a+1:) = out
        h(a+1:, :a) = TRANSPOSE(out)

        new = FGH_t(this%f * that%f, [that%f * this%g, this%f * that%g], h)
    end function FGH_Mul


    elemental type(FGH_t) function FGH_MulSame(this, that) result(new)
        class(FGH_t), intent(in) :: this, that
        real(r8) :: out(SIZE(this%g), SIZE(that%g))
        integer :: a

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(that%g, 1, a)

        new%f = this%f * that%f
        new%g = this%f * that%g + that%f * this%g
        new%h = this%h * that%f + that%h * this%f + out + TRANSPOSE(out)

    end function FGH_MulSame


    elemental type(FGH_t) function FGH_Pow(this, that) result(new)
        class(FGH_t), intent(in) :: this, that
        real(r8) :: outxy(SIZE(this%g), SIZE(that%g)), outxx(SIZE(this%g), SIZE(that%g))
        real(r8) :: outyy(SIZE(this%g), SIZE(that%g)), logval, frac
        integer :: a

        a = SIZE(this%g)
        outxy = SPREAD(this%g, 2, a) * SPREAD(that%g, 1, a)
        outxx = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        outyy = SPREAD(that%g, 2, a) * SPREAD(that%g, 1, a)

        logval = LOG(this%f)
        frac = (that%f/this%f)
        new%f = this%f ** that%f
        new%g = new%f * (that%g * LOG(this%f) + (that%f / this%f) * this%g)
        new%h = new%f * (outyy * logval**2 + outxx * frac**2 + frac * logval * (TRANSPOSE(outxy) + outxy) +&
                that%h * logval + (outxy + TRANSPOSE(outxy)) / this%f - outxx / this%f**2 + this%h * frac)

    end function FGH_Pow

#define _FGH true
    ! single precision
#define _F_TYPE FGH_t
#define _F FGH
#define _KIND r4
#define _TYPE real(r4)
#include "F_functions.inc"



    ! double precision
#define _F_TYPE FGH_t
#define _F FGH
#define _KIND r8
#define _TYPE real(r8)
#include "F_functions.inc"



    ! quad precision
#define _F_TYPE FGH_t
#define _F FGH
#define _KIND r16
#define _TYPE real(r16)
#include "F_functions.inc"



    ! 8bit int
#define _INT true
#define _F_TYPE FGH_t
#define _F FGH
#define _KIND i1
#define _TYPE integer(i1)
#include "F_functions.inc"



    ! 16bit int
#define _INT true
#define _F_TYPE FGH_t
#define _F FGH
#define _KIND i2
#define _TYPE integer(i2)
#include "F_functions.inc"



    ! 32bit int
#define _INT true
#define _F_TYPE FGH_t
#define _F FGH
#define _KIND i4
#define _TYPE integer(i4)
#include "F_functions.inc"



    ! 64bit int
#define _INT true
#define _F_TYPE FGH_t
#define _F FGH
#define _KIND i8
#define _TYPE integer(i8)
#include "F_functions.inc"


#undef _FGH

end module FGH_m