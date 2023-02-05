module FGL_m
    use, intrinsic :: ISO_FORTRAN_ENV, only: r8 => REAL64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_NEGATIVE_INF, IEEE_QUIET_NAN

    implicit none

    private
    public :: FGL_t, operator(+), operator(-), operator(*), r8, operator(**), operator(/), ABS, SQRT, EXP, LOG, COS, SIN
    public :: TAN, INT, NINT, FRACTION, REAL, FLOOR

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
        module procedure :: FGL_AddScalarObj, FGL_AddObjScalar
    end interface

    interface operator(-)
        module procedure :: FGL_SubScalarObj, FGL_SubObjScalar
    end interface

    interface operator(*)
        module procedure :: FGL_MulScalarObj, FGL_MulObjScalar
    end interface

    interface operator(/)
        module procedure :: FGL_DivScalarObj, FGL_DivObjScalar
    end interface

    interface operator(**)
        module procedure :: FGL_Pow_int, FGL_Pow_real
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

        new = this%FGL_MulSame(FGL_Pow_int(that, -1))
    end function FGL_DivSame


    elemental type(FGL_t) function FGL_Abs_real(this) result(new)
        class(FGL_t), intent(in) :: this

        if (this%f < 0) then
            new = -this
        else
            new = this
        end if
    end function FGL_Abs_real


    elemental type(FGL_t) function FGL_sqrt_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new = FGL_pow_real(this, 0.5_r8)

    end function FGL_sqrt_real


    elemental type(FGL_t) function FGL_exp(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new

        f_new = EXP(this%f)
        new = FGL_t(f_new, f_new * this%g, f_new * (DOT_PRODUCT(this%g, this%g) + this%l))

    end function FGL_exp


    elemental type(FGL_t) function FGL_log_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new
        real(r8)    :: dummy1(SIZE(this%g))

        if (this%f > 0) then
            f_new = LOG(this%f)
            new = FGL_t(f_new, this%g / this%f, -1 / (this%f ** 2) * DOT_PRODUCT(this%g, this%g) + 1 / this%f * this%l)
        else
            dummy1 = IEEE_VALUE(1._r8, IEEE_QUIET_NAN)
            new = FGL_t(IEEE_VALUE(1._r8, IEEE_NEGATIVE_INF), dummy1, IEEE_VALUE(1._r8, IEEE_QUIET_NAN))
        end if

    end function FGL_log_real


    elemental type(FGL_t) function FGL_cos_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new

        f_new = COS(this%f)
        new = FGL_t(f_new, - SIN(this%f) * this%g, -f_new * DOT_PRODUCT(this%g, this%g) - SIN(this%f) * this%l)

    end function FGL_cos_real


    elemental type(FGL_t) function FGL_sin_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new

        f_new = SIN(this%f)
        new = FGL_t(f_new, this%g * COS(this%f), -f_new * DOT_PRODUCT(this%g, this%g) + COS(this%f) * this%l)

    end function FGL_sin_real


    elemental type(FGL_t) function FGL_tan_real(this) result(new)
        class(FGL_t), intent(in) :: this
        real(r8)    :: f_new, sec2

        f_new = TAN(this%f)
        sec2 = 1 / cos(this%f)**2
        new = FGL_t(f_new, this%g * sec2, sec2 * (2 * f_new * DOT_PRODUCT(this%g, this%g) + this%l))

    end function FGL_tan_real


    elemental type(FGL_t) function FGL_int_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new = FGL_t(INT(this%f), this%g, this%l)

    end function FGL_int_real


    elemental type(FGL_t) function FGL_nint_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new = FGL_t(NINT(this%f), this%g, this%l)

    end function FGL_nint_real


    elemental type(FGL_t) function FGL_floor_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new = FGL_t(FLOOR(this%f), this%g, this%l)

    end function FGL_floor_real


    elemental type(FGL_t) function FGL_fraction_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new = FGL_t(FRACTION(this%f), this%g, this%l)

    end function FGL_fraction_real


    elemental type(FGL_t) function FGL_real_real(this) result(new)
        class(FGL_t), intent(in) :: this

        new = FGL_t(REAL(this%f), this%g, this%l)

    end function FGL_real_real


    elemental type(FGL_t) function FGL_Div(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new = this%FGL_Mul(FGL_Pow_int(that, -1))
    end function FGL_Div


    elemental type(FGL_t) function FGL_DivObjScalar(obj, value) result(new)
        class(FGL_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGL_t(obj%f / value, obj%g / value, obj%l / value)
    end function FGL_DivObjScalar


    elemental type(FGL_t) function FGL_DivScalarObj(value, obj) result(new)
        class(FGL_t), intent(in) :: obj
        real(r8), intent(in) :: value
        type(FGL_t) :: temp

        temp = FGL_Pow_int(obj, -1)
        new = FGL_t(value * temp%f, value * temp%g, value * temp%l)
    end function FGL_DivScalarObj


    pure subroutine FGL_Assign(this, that)
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

        new = FGL_t(- this%f, - this%g, - this%l)
    end function FGL_Neg


    elemental type(FGL_t) function FGL_AddObjScalar(obj, value) result(new)
        class(FGL_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGL_t(obj%f + value, obj%g, obj%l)
    end function FGL_AddObjScalar


    elemental type(FGL_t) function FGL_AddScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FGL_t), intent(in) :: obj

        new = FGL_t(value + obj%f, obj%g, obj%l)
    end function FGL_AddScalarObj


    elemental type(FGL_t) function FGL_Add(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new = FGL_t(this%f + that%f, this%g + that%g, this%l + that%l)
    end function FGL_Add


    elemental type(FGL_t) function FGL_Pow_real(this, n) result(new)
        class(FGL_t), intent(in) :: this
        real(r8), intent(in) :: n

        new = FGL_t(this%f ** n, this%f ** (n-1) * n * this%g, &
                n * ((n-1) * this%f ** (n-2) * DOT_PRODUCT(this%g, this%g) + this%f ** (n-1) * this%l))
    end function FGL_Pow_real


    elemental type(FGL_t) function FGL_Pow_int(this, n) result(new)
        class(FGL_t), intent(in) :: this
        integer, intent(in) :: n

        new = FGL_t(this%f ** n, this%f ** (n-1) * n * this%g, &
                n * ((n-1) * this%f ** (n-2) * DOT_PRODUCT(this%g, this%g) + this%f ** (n-1) * this%l))
    end function FGL_Pow_int


    elemental type(FGL_t) function FGL_SubObjScalar(obj, value) result(new)
        class(FGL_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGL_t(obj%f - value, obj%g, obj%l)
    end function FGL_SubObjScalar


    elemental type(FGL_t) function FGL_SubScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FGL_t), intent(in) :: obj

        new = FGL_t(value - obj%f, obj%g, obj%l)
    end function FGL_SubScalarObj


    elemental type(FGL_t) function FGL_Sub(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new = FGL_t(this%f - that%f, this%g - that%g, this%l - that%l)
    end function FGL_Sub


    elemental type(FGL_t) function FGL_MulObjScalar(obj, value) result(new)
        class(FGL_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGL_t(obj%f * value, obj%g * value, obj%l * value)
    end function FGL_MulObjScalar


    elemental type(FGL_t) function FGL_MulScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FGL_t), intent(in) :: obj

        new = FGL_t(value * obj%f, value * obj%g, value * obj%l)
    end function FGL_MulScalarObj


    elemental type(FGL_t) function FGL_Mul(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new = FGL_t(this%f * that%f, [that%f * this%g, this%f * that%g], &
        2 * DOT_PRODUCT(this%g, that%g) + this%f * that%l + that%f * this%l)
    end function FGL_Mul


    elemental type(FGL_t) function FGL_MulSame(this, that) result(new)
        class(FGL_t), intent(in) :: this, that

        new = FGL_t(this%f * that%f, this%f * that%g + that%f * this%g, &
                this%l * that%f + that%l * this%f + 2 * DOT_PRODUCT(this%g, that%g))
    end function FGL_MulSame
end module FGL_m