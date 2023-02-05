module FG_m
    use, intrinsic :: ISO_FORTRAN_ENV, only: r8 => REAL64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_NEGATIVE_INF, IEEE_QUIET_NAN

    implicit none

    private
    public :: FG_t, operator(+), operator(-), operator(*), r8, operator(**), operator(/), ABS, SQRT, EXP, LOG, COS, SIN
    public :: TAN, INT, NINT, FRACTION, REAL, FLOOR

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
        module procedure :: FG_AddScalarObj, FG_AddObjScalar
    end interface

    interface operator(-)
        module procedure :: FG_SubScalarObj, FG_SubObjScalar
    end interface

    interface operator(*)
        module procedure :: FG_MulScalarObj, FG_MulObjScalar
    end interface

    interface operator(/)
        module procedure :: FG_DivScalarObj, FG_DivObjScalar
    end interface

    interface operator(**)
        module procedure :: FG_Pow_int, FG_Pow_real
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

contains

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

        new = this%FG_MulSame(FG_Pow_int(that, -1))
    end function FG_DivSame


    elemental type(FG_t) function FG_Abs_real(this) result(new)
        class(FG_t), intent(in) :: this

        if (this%f < 0) then
            new = -this
        else
            new = this
        end if
    end function FG_Abs_real


    elemental type(FG_t) function FG_sqrt_real(this) result(new)
        class(FG_t), intent(in) :: this

        new = FG_pow_real(this, 0.5_r8)

    end function FG_sqrt_real


    elemental type(FG_t) function FG_exp_real(this) result(new)
        class(FG_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new

        f_new = EXP(this%f)
        new = FG_t(f_new, f_new * this%g)

    end function FG_exp_real


    elemental type(FG_t) function FG_log_real(this) result(new)
        class(FG_t), intent(in) :: this
        real(r8)    :: f_new
        real(r8)    :: dummy1(SIZE(this%g))

        if (this%f > 0) then
            f_new = LOG(this%f)
            new = FG_t(f_new, this%g / this%f)
        else
            dummy1 = IEEE_VALUE(1._r8, IEEE_QUIET_NAN)
            new = FG_t(IEEE_VALUE(1._r8, IEEE_NEGATIVE_INF), dummy1)
        end if

    end function FG_log_real


    elemental type(FG_t) function FG_cos_real(this) result(new)
        class(FG_t), intent(in) :: this
        real(r8)    :: f_new

        f_new = COS(this%f)
        new = FG_t(f_new, - SIN(this%f) * this%g)

    end function FG_cos_real


    elemental type(FG_t) function FG_sin_real(this) result(new)
        class(FG_t), intent(in) :: this
        real(r8)    :: f_new

        f_new = SIN(this%f)
        new = FG_t(f_new, this%g * COS(this%f))

    end function FG_sin_real


    elemental type(FG_t) function FG_tan_real(this) result(new)
        class(FG_t), intent(in) :: this
        real(r8)    :: f_new, sec

        f_new = TAN(this%f)
        sec = 1 / cos(this%f)**2
        new = FG_t(f_new, this%g * sec)

    end function FG_tan_real


    elemental type(FG_t) function FG_int_real(this) result(new)
        class(FG_t), intent(in) :: this

        new = FG_t(INT(this%f), this%g)

    end function FG_int_real


    elemental type(FG_t) function FG_nint_real(this) result(new)
        class(FG_t), intent(in) :: this

        new = FG_t(NINT(this%f), this%g)

    end function FG_nint_real


    elemental type(FG_t) function FG_floor_real(this) result(new)
        class(FG_t), intent(in) :: this

        new = FG_t(FLOOR(this%f), this%g)

    end function FG_floor_real


    elemental type(FG_t) function FG_fraction_real(this) result(new)
        class(FG_t), intent(in) :: this

        new = FG_t(FRACTION(this%f), this%g)

    end function FG_fraction_real


    elemental type(FG_t) function FG_real_real(this) result(new)
        class(FG_t), intent(in) :: this

        new = FG_t(REAL(this%f), this%g)

    end function FG_real_real


    elemental type(FG_t) function FG_Div(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new = this%FG_Mul(FG_Pow_int(that, -1))
    end function FG_Div


    elemental type(FG_t) function FG_DivObjScalar(obj, value) result(new)
        class(FG_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FG_t(obj%f / value, obj%g / value)
    end function FG_DivObjScalar


    elemental type(FG_t) function FG_DivScalarObj(value, obj) result(new)
        class(FG_t), intent(in) :: obj
        real(r8), intent(in) :: value
        type(FG_t) :: temp

        temp = FG_Pow_int(obj, -1)
        new = FG_t(value * temp%f, value * temp%g)
    end function FG_DivScalarObj


    pure subroutine FG_Assign(this, that)
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

        new = FG_t(- this%f, - this%g)
    end function FG_Neg


    elemental type(FG_t) function FG_AddObjScalar(obj, value) result(new)
        class(FG_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FG_t(obj%f + value, obj%g)
    end function FG_AddObjScalar


    elemental type(FG_t) function FG_AddScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FG_t), intent(in) :: obj

        new = FG_t(value + obj%f, obj%g)
    end function FG_AddScalarObj


    elemental type(FG_t) function FG_Add(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new = FG_t(this%f + that%f, this%g + that%g)
    end function FG_Add


    elemental type(FG_t) function FG_Pow_real(this, n) result(new)
        class(FG_t), intent(in) :: this
        real(r8), intent(in) :: n

        new = FG_t(this%f ** n, this%f ** (n-1) * n * this%g)
    end function FG_Pow_real


    elemental type(FG_t) function FG_Pow_int(this, n) result(new)
        class(FG_t), intent(in) :: this
        integer, intent(in) :: n

        new = FG_t(this%f ** n, this%f ** (n-1) * n * this%g)
    end function FG_Pow_int


    elemental type(FG_t) function FG_SubObjScalar(obj, value) result(new)
        class(FG_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FG_t(obj%f - value, obj%g)
    end function FG_SubObjScalar


    elemental type(FG_t) function FG_SubScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FG_t), intent(in) :: obj

        new = FG_t(value - obj%f, obj%g)
    end function FG_SubScalarObj


    elemental type(FG_t) function FG_Sub(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new = FG_t(this%f - that%f, this%g - that%g)
    end function FG_Sub


    elemental type(FG_t) function FG_MulObjScalar(obj, value) result(new)
        class(FG_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FG_t(obj%f * value, obj%g * value)
    end function FG_MulObjScalar


    elemental type(FG_t) function FG_MulScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FG_t), intent(in) :: obj

        new = FG_t(value * obj%f, value * obj%g)
    end function FG_MulScalarObj


    elemental type(FG_t) function FG_Mul(this, that) result(new)
        class(FG_t), intent(in) :: this, that
        integer :: a, b

        new = FG_t(this%f * that%f, [that%f * this%g, this%f * that%g])
    end function FG_Mul


    elemental type(FG_t) function FG_MulSame(this, that) result(new)
        class(FG_t), intent(in) :: this, that

        new = FG_t(this%f * that%f, this%f * that%g + that%f * this%g)
    end function FG_MulSame
end module FG_m