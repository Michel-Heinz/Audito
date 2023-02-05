module FGH_m
    use, intrinsic :: ISO_FORTRAN_ENV, only: r8 => REAL64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_NEGATIVE_INF, IEEE_QUIET_NAN

    implicit none

    private
    public :: FGH_t, operator(+), operator(-), operator(*), r8, operator(**), operator(/), ABS, SQRT, EXP, LOG, COS, SIN
    public :: TAN, INT, NINT, FRACTION, REAL, FLOOR

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
        module procedure :: FGH_AddScalarObj, FGH_AddObjScalar
    end interface

    interface operator(-)
        module procedure :: FGH_SubScalarObj, FGH_SubObjScalar
    end interface

    interface operator(*)
        module procedure :: FGH_MulScalarObj, FGH_MulObjScalar
    end interface

    interface operator(/)
        module procedure :: FGH_DivScalarObj, FGH_DivObjScalar
    end interface

    interface operator(**)
        module procedure :: FGH_Pow_int, FGH_Pow_real
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

contains

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

        new = this%FGH_MulSame(FGH_Pow_int(that, -1))
    end function FGH_DivSame


    elemental type(FGH_t) function FGH_Abs_real(this) result(new)
        class(FGH_t), intent(in) :: this

        if (this%f < 0) then
            new = -this
        else
            new = this
        end if
    end function FGH_Abs_real


    elemental type(FGH_t) function FGH_sqrt_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new = FGH_pow_real(this, 0.5_r8)

    end function FGH_sqrt_real


    elemental type(FGH_t) function FGH_exp_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g))

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        f_new = EXP(this%f)
        new = FGH_t(f_new, f_new * this%g, f_new * (out + this%h))

    end function FGH_exp_real


    elemental type(FGH_t) function FGH_log_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g))
        real(r8)    :: dummy1(SIZE(this%g)), dummy2(SIZE(this%g), SIZE(this%g))

        a = SIZE(this%g)
        if (this%f > 0) then
            out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
            f_new = LOG(this%f)
            new = FGH_t(f_new, this%g / this%f, (this%f * this%h - out) / this%f**2)
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
        f_new = COS(this%f)
        new = FGH_t(f_new, - SIN(this%f) * this%g, -f_new * out - SIN(this%f) * this%h)

    end function FGH_cos_real


    elemental type(FGH_t) function FGH_sin_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g))

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        f_new = SIN(this%f)
        new = FGH_t(f_new, this%g * COS(this%f), -f_new * out + COS(this%f) * this%h)

    end function FGH_sin_real


    elemental type(FGH_t) function FGH_tan_real(this) result(new)
        class(FGH_t), intent(in) :: this
        integer     :: a
        real(r8)    :: f_new, out(SIZE(this%g), SIZE(this%g)), sec

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)
        f_new = TAN(this%f)
        sec = 1 / cos(this%f)**2
        new = FGH_t(f_new, this%g * sec, sec * (2 * f_new * out + this%h))

    end function FGH_tan_real


    elemental type(FGH_t) function FGH_int_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new = FGH_t(INT(this%f), this%g, this%h)

    end function FGH_int_real


    elemental type(FGH_t) function FGH_nint_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new = FGH_t(NINT(this%f), this%g, this%h)

    end function FGH_nint_real


    elemental type(FGH_t) function FGH_floor_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new = FGH_t(FLOOR(this%f), this%g, this%h)

    end function FGH_floor_real


    elemental type(FGH_t) function FGH_fraction_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new = FGH_t(FRACTION(this%f), this%g, this%h)

    end function FGH_fraction_real


    elemental type(FGH_t) function FGH_real_real(this) result(new)
        class(FGH_t), intent(in) :: this

        new = FGH_t(REAL(this%f), this%g, this%h)

    end function FGH_real_real


    elemental type(FGH_t) function FGH_Div(this, that) result(new)
        class(FGH_t), intent(in) :: this, that

        new = this%FGH_Mul(FGH_Pow_int(that, -1))
    end function FGH_Div


    elemental type(FGH_t) function FGH_DivObjScalar(obj, value) result(new)
        class(FGH_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGH_t(obj%f / value, obj%g / value, obj%h / value)
    end function FGH_DivObjScalar


    elemental type(FGH_t) function FGH_DivScalarObj(value, obj) result(new)
        class(FGH_t), intent(in) :: obj
        real(r8), intent(in) :: value
        type(FGH_t) :: temp

        temp = FGH_Pow_int(obj, -1)
        new = FGH_t(value * temp%f, value * temp%g, value * temp%h)
    end function FGH_DivScalarObj


    pure subroutine FGH_Assign(this, that)
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

        new = FGH_t(- this%f, - this%g, - this%h)
    end function FGH_Neg


    elemental type(FGH_t) function FGH_AddObjScalar(obj, value) result(new)
        class(FGH_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGH_t(obj%f + value, obj%g, obj%h)
    end function FGH_AddObjScalar


    elemental type(FGH_t) function FGH_AddScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FGH_t), intent(in) :: obj

        new = FGH_t(value + obj%f, obj%g, obj%h)
    end function FGH_AddScalarObj


    elemental type(FGH_t) function FGH_Add(this, that) result(new)
        class(FGH_t), intent(in) :: this, that

        new = FGH_t(this%f + that%f, this%g + that%g, this%h + that%h)
    end function FGH_Add


    elemental type(FGH_t) function FGH_Pow_real(this, n) result(new)
        class(FGH_t), intent(in) :: this
        real(r8), intent(in) :: n
        real(r8), dimension(SIZE(this%g), SIZE(this%g)) :: out
        integer :: a

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)

        new = FGH_t(this%f ** n, this%f ** (n-1) * n * this%g, &
                n * ((n-1) * this%f ** (n-2) * out + this%f ** (n-1) * this%h))
    end function FGH_Pow_real


    elemental type(FGH_t) function FGH_Pow_int(this, n) result(new)
        class(FGH_t), intent(in) :: this
        integer, intent(in) :: n
        real(r8), dimension(SIZE(this%g), SIZE(this%g)) :: out
        integer :: a

        a = SIZE(this%g)
        out = SPREAD(this%g, 2, a) * SPREAD(this%g, 1, a)

        new = FGH_t(this%f ** n, this%f ** (n-1) * n * this%g, &
                n * ((n-1) * this%f ** (n-2) * out + this%f ** (n-1) * this%h))
    end function FGH_Pow_int


    elemental type(FGH_t) function FGH_SubObjScalar(obj, value) result(new)
        class(FGH_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGH_t(obj%f - value, obj%g, obj%h)
    end function FGH_SubObjScalar


    elemental type(FGH_t) function FGH_SubScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FGH_t), intent(in) :: obj

        new = FGH_t(value - obj%f, obj%g, obj%h)
    end function FGH_SubScalarObj


    elemental type(FGH_t) function FGH_Sub(this, that) result(new)
        class(FGH_t), intent(in) :: this, that

        new = FGH_t(this%f - that%f, this%g - that%g, this%h - that%h)
    end function FGH_Sub


    elemental type(FGH_t) function FGH_MulObjScalar(obj, value) result(new)
        class(FGH_t), intent(in) :: obj
        real(r8), intent(in) :: value

        new = FGH_t(obj%f * value, obj%g * value, obj%h * value)
    end function FGH_MulObjScalar


    elemental type(FGH_t) function FGH_MulScalarObj(value, obj) result(new)
        real(r8), intent(in) :: value
        class(FGH_t), intent(in) :: obj

        new = FGH_t(value * obj%f, value * obj%g, value * obj%h)
    end function FGH_MulScalarObj


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

        out = SPREAD(this%g, 2, SIZE(this%g)) * SPREAD(that%g, 1, SIZE(that%g))

        new = FGH_t(this%f * that%f, this%f * that%g + that%f * this%g, &
                this%h * that%f + that%h * this%f + out + TRANSPOSE(out))
    end function FGH_MulSame
end module FGH_m