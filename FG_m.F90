! Copyright (C) 2023 Michel Heinz

#  define _PASTE(X) X
#  define _CAT2(X,Y) _PASTE(X)Y
#  define _CAT3(X,Y,Z) _CAT2(_CAT2(X,Y),Z)


module FG_m
    use, intrinsic :: ISO_FORTRAN_ENV, only: REAL64, REAL32, REAL128, INT8, INT16, INT32, INT64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_NEGATIVE_INF, IEEE_QUIET_NAN

    implicit none

    private
    public :: FG_t, operator(+), operator(-), operator(*), operator(**), operator(/), ABS, SQRT, EXP, LOG, COS, SIN
    public :: TAN, INT, NINT, FRACTION, REAL, FLOOR
    public :: i1, i2, i4, i8, r4, r8, r16

    integer(INT32), parameter :: i1 = INT8
    integer(INT32), parameter :: i2 = INT16
    integer(INT32), parameter :: i4 = INT32
    integer(INT32), parameter :: i8 = INT64
    integer(INT32), parameter :: r4 = REAL32
    integer(INT32), parameter :: r8 = REAL64
    integer(INT32), parameter :: r16 = REAL128

    type FG_t
        real(r8) :: f
        real(r8), allocatable :: g(:)
    contains
        procedure :: print => FG_Print
        procedure, private :: FG_Neg
        generic :: operator(-) => FG_Neg
        procedure, private :: FG_Add
        generic :: operator(+) => FG_Add
        procedure, private :: FG_Sub
        generic :: operator(-) => FG_Sub
        procedure, private :: FG_MulSame
        generic :: operator(*) => FG_MulSame
        procedure, private :: FG_DivSame
        generic :: operator(/) => FG_DivSame
        procedure, private :: FG_Pow
        generic :: operator(**) => FG_Pow
        procedure, private :: FG_Equal
        generic :: operator(==) => FG_Equal
        procedure, private :: FG_NotEqual
        generic :: operator(.ne.) => FG_NotEqual
        procedure, private :: FG_LessThan
        generic :: operator(<) => FG_LessThan
        procedure, private :: FG_GreaterThan
        generic :: operator(>) => FG_GreaterThan
        procedure, private :: FG_LessThanOrEqual
        generic :: operator(<=) => FG_LessThanOrEqual
        procedure, private :: FG_GreaterThanOrEqual
        generic :: operator(>=) => FG_GreaterThanOrEqual
        procedure, private :: FG_Assign
        generic :: assignment(=) => FG_Assign
        procedure, private :: FG_Mul
        generic :: operator(.mul.) => FG_Mul
        procedure, private :: FG_Div
        generic :: operator(.div.) => FG_Div
        procedure          :: destroy
    end type FG_t

    interface operator(+)
        module procedure :: FG_AddScalarObj_r8, FG_AddObjScalar_r8, FG_AddScalarObj_r4, FG_AddObjScalar_r4
        module procedure :: FG_AddScalarObj_r16, FG_AddObjScalar_r16, FG_AddScalarObj_i1, FG_AddObjScalar_i1
        module procedure :: FG_AddScalarObj_i2, FG_AddObjScalar_i2, FG_AddScalarObj_i4, FG_AddObjScalar_i4
        module procedure :: FG_AddScalarObj_i8, FG_AddObjScalar_i8
    end interface

    interface operator(-)
        module procedure :: FG_SubScalarObj_r8, FG_SubObjScalar_r8, FG_SubScalarObj_r4, FG_SubObjScalar_r4
        module procedure :: FG_SubScalarObj_r16, FG_SubObjScalar_r16, FG_SubScalarObj_i1, FG_SubObjScalar_i1
        module procedure :: FG_SubScalarObj_i2, FG_SubObjScalar_i2, FG_SubScalarObj_i4, FG_SubObjScalar_i4
        module procedure :: FG_SubScalarObj_i8, FG_SubObjScalar_i8
    end interface

    interface operator(*)
        module procedure :: FG_MulScalarObj_r8, FG_MulObjScalar_r8, FG_MulScalarObj_r4, FG_MulObjScalar_r4
        module procedure :: FG_MulScalarObj_r16, FG_MulObjScalar_r16, FG_MulScalarObj_i1, FG_MulObjScalar_i1
        module procedure :: FG_MulScalarObj_i2, FG_MulObjScalar_i2, FG_MulScalarObj_i4, FG_MulObjScalar_i4
        module procedure :: FG_MulScalarObj_i8, FG_MulObjScalar_i8
    end interface

    interface operator(/)
        module procedure :: FG_DivScalarObj_r8, FG_DivObjScalar_r8, FG_DivScalarObj_r4, FG_DivObjScalar_r4
        module procedure :: FG_DivScalarObj_r16, FG_DivObjScalar_r16, FG_DivScalarObj_i1, FG_DivObjScalar_i1
        module procedure :: FG_DivScalarObj_i2, FG_DivObjScalar_i2, FG_DivScalarObj_i4, FG_DivObjScalar_i4
        module procedure :: FG_DivScalarObj_i8, FG_DivObjScalar_i8
    end interface

    interface operator(**)
        module procedure :: FG_PowObjScalar_r4, FG_PowObjScalar_r8, FG_PowObjScalar_r16, FG_PowObjScalar_i1
        module procedure :: FG_PowObjScalar_i2, FG_PowObjScalar_i4, FG_PowObjScalar_i8
        module procedure :: FG_PowScalarObj_r4, FG_PowScalarObj_r8, FG_PowScalarObj_r16, FG_PowScalarObj_i1
        module procedure :: FG_PowScalarObj_i2, FG_PowScalarObj_i4, FG_PowScalarObj_i8

    end interface

    interface ABS
        module procedure :: FG_Abs_real
    end interface

    interface SQRT
        module procedure :: FG_sqrt_real
    end interface

    interface EXP
        module procedure :: FG_exp_real
    end interface

    interface LOG
        module procedure :: FG_log_real
    end interface

    interface COS
        module procedure :: FG_cos_real
    end interface

    interface SIN
        module procedure :: FG_sin_real
    end interface

    interface TAN
        module procedure :: FG_tan_real
    end interface

    interface INT
        module procedure :: FG_int_real
    end interface

    interface NINT
        module procedure :: FG_nint_real
    end interface

    interface FLOOR
        module procedure :: FG_floor_real
    end interface

    interface FRACTION
        module procedure :: FG_fraction_real
    end interface

    interface REAL
        module procedure :: FG_real_real
    end interface

!    interface FG_t
!        module procedure :: constructor1
!    end interface FG_t
contains
!    pure type(FG_t) function constructor1(f, g) result(new)
!        real(r8), intent(in) :: f, g(:)
!
!        new%f = f
!        new%g = g
!    end function constructor1

    subroutine destroy(this)
        class(FG_t), intent(inout) :: this

        if (allocated(this%g)) deallocate(this%g)
    end subroutine

    subroutine FG_Print(this, iu)
        class(FG_t), intent(in)   :: this
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
        else
            print'(a)', 'value:'
            print'(e13.5)', this%f
            print'(a)', 'gradient:'
            print'('//TRIM(n_str)//'e11.3)', this%g
        end if
    end subroutine FG_Print


    elemental type(FG_t) function FG_DivSame(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new%f = this%f / that%f
        new%g = (this%g * that%f - that%g * this%f) / that%f ** 2
    end function FG_DivSame


    elemental type(FG_t) function FG_Abs_real(this) result(new)
        class(FG_t), intent(in) :: this

        if (this%f < 0) then
            new%f = -this%f
            new%g = -this%g
        else
            new%f = this%f
            new%g = this%g
        end if
    end function FG_Abs_real


    elemental type(FG_t) function FG_sqrt_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = SQRT(this%f)
        new%g = this%g / (2 * SQRT(this%f))
    end function FG_sqrt_real


    elemental type(FG_t) function FG_exp_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = EXP(this%f)
        new%g = new%f * this%g
    end function FG_exp_real


    elemental type(FG_t) function FG_log_real(this) result(new)
        class(FG_t), intent(in) :: this
        real(r8)    :: dummy1(SIZE(this%g))

        if (this%f > 0) then
            new%f = LOG(this%f)
            new%g = this%g / this%f
        else
            dummy1 = IEEE_VALUE(1._r8, IEEE_QUIET_NAN)
            new = FG_t(IEEE_VALUE(1._r8, IEEE_NEGATIVE_INF), dummy1)
        end if

    end function FG_log_real


    elemental type(FG_t) function FG_cos_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = COS(this%f)
        new%g = - SIN(this%f) * this%g

    end function FG_cos_real


    elemental type(FG_t) function FG_sin_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = SIN(this%f)
        new%g = this%g * COS(this%f)

    end function FG_sin_real


    elemental type(FG_t) function FG_tan_real(this) result(new)
        class(FG_t), intent(in) :: this
        real(r8)    :: sec

        new%f = TAN(this%f)
        sec = 1 / cos(this%f)**2
        new%g = this%g * sec

    end function FG_tan_real


    elemental type(FG_t) function FG_int_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = INT(this%f)
        new%g = this%g

    end function FG_int_real


    elemental type(FG_t) function FG_nint_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = NINT(this%f)
        new%g = this%g

    end function FG_nint_real


    elemental type(FG_t) function FG_floor_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = FLOOR(this%f)
        new%g = this%g

    end function FG_floor_real


    elemental type(FG_t) function FG_fraction_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = FRACTION(this%f)
        new%g = this%g

    end function FG_fraction_real


    elemental type(FG_t) function FG_real_real(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = REAL(this%f)
        new%g = this%g

    end function FG_real_real


    elemental type(FG_t) function FG_Div(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new = this%FG_Mul(FG_PowObjScalar_i1(that, -1_i1))
    end function FG_Div


    elemental subroutine FG_Assign(this, that)
        class(FG_t), intent(inout) :: this
        class(FG_t), intent(in)   :: that

        this%f = that%f
        this%g = that%g

    end subroutine FG_Assign


    elemental logical function FG_Equal(this, that) result(equal)
        class(FG_t), intent(in) :: this, that

        equal = .false.
        if (this%f == that%f .and. ALL(this%g == that%g)) &
                equal = .true.
    end function FG_Equal


    elemental logical function FG_NotEqual(this, that) result(notEqual)
        class(FG_t), intent(in) :: this, that

        notEqual = .false.
        if (this%f /= that%f .and. ALL(this%g /= that%g)) &
                notEqual = .true.
    end function FG_NotEqual


    elemental logical function FG_LessThan(this, that) result(LessThan)
        class(FG_t), intent(in) :: this, that

        LessThan = .false.
        if (this%f < that%f) LessThan = .true.
    end function FG_LessThan


    elemental logical function FG_GreaterThan(this, that) result(GreaterThan)
        class(FG_t), intent(in) :: this, that

        GreaterThan = .false.
        if (this%f > that%f) GreaterThan = .true.
    end function FG_GreaterThan


    elemental logical function FG_LessThanOrEqual(this, that) result(LessThanOrEqual)
        class(FG_t), intent(in) :: this, that

        LessThanOrEqual = .false.
        if (this%f <= that%f) LessThanOrEqual = .true.
    end function FG_LessThanOrEqual


    elemental logical function FG_GreaterThanOrEqual(this, that) result(GreaterThanOrEqual)
        class(FG_t), intent(in) :: this, that

        GreaterThanOrEqual = .false.
        if (this%f >= that%f) GreaterThanOrEqual = .true.
    end function FG_GreaterThanOrEqual


    elemental type(FG_t) function FG_Neg(this) result(new)
        class(FG_t), intent(in) :: this

        new%f = -this%f
        new%g = -this%g
    end function FG_Neg


    elemental type(FG_t) function FG_Add(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new%f = this%f + that%f
        new%g = this%g + that%g
    end function FG_Add


    elemental type(FG_t) function FG_Sub(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new%f = this%f - that%f
        new%g = this%g - that%g
    end function FG_Sub


    elemental type(FG_t) function FG_Mul(this, that) result(new)
        class(FG_t), intent(in) :: this, that
        integer :: a, b

        new = FG_t(this%f * that%f, [that%f * this%g, this%f * that%g])
    end function FG_Mul


    elemental type(FG_t) function FG_MulSame(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new%f = this%f * that%f
        new%g = this%g * that%f + that%g * this%f
    end function FG_MulSame


    elemental type(FG_t) function FG_Pow(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new%f = this%f ** that%f
        new%g = new%f * (that%g * LOG(this%f) + (that%f / this%f) * this%g)

    end function FG_Pow


    ! single precision
#define _F_TYPE FG_t
#define _F FG
#define _KIND r4
#define _TYPE real(r4)
#include "F_functions.inc"



    ! double precision
#define _F_TYPE FG_t
#define _F FG
#define _KIND r8
#define _TYPE real(r8)
#include "F_functions.inc"



    ! quad precision
#define _F_TYPE FG_t
#define _F FG
#define _KIND r16
#define _TYPE real(r16)
#include "F_functions.inc"



    ! 8bit int
#define _INT true
#define _F_TYPE FG_t
#define _F FG
#define _KIND i1
#define _TYPE integer(i1)
#include "F_functions.inc"



    ! 16bit int
#define _INT true
#define _F_TYPE FG_t
#define _F FG
#define _KIND i2
#define _TYPE integer(i2)
#include "F_functions.inc"



    ! 32bit int
#define _INT true
#define _F_TYPE FG_t
#define _F FG
#define _KIND i4
#define _TYPE integer(i4)
#include "F_functions.inc"



    ! 64bit int
#define _INT true
#define _F_TYPE FG_t
#define _F FG
#define _KIND i8
#define _TYPE integer(i8)
#include "F_functions.inc"


end module FG_m