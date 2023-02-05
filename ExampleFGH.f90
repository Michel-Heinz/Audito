program example
    use, intrinsic :: ISO_FORTRAN_ENV, only: r8 => REAL64
    use :: FGH_m

    implicit none
    integer :: numAtoms, seed, i, steps, maxIterations, maxLss, verbosicity
    type(FGH_t), allocatable :: coords(:)
    type(FGH_t) :: ePot
    real(r8)    :: t_0, t, stepSize
    real(r8), allocatable :: iniGrad(:), iniHess(:,:)

    !set the program variables:
    numAtoms      = 5         !number of atoms
    seed          = 1658717   !seed for the generation of initial coordinates
    maxIterations = 100000    !maximum number of iterations for the gradient descent algorithm
    maxLss        = 100       !maximum number of line search steps
    stepSize      = 1.0_r8    !initial step size for the golden section line search
    verbosicity   = 2         !0 for no printouts, 1 for minimal printouts, 2 for full printouts

    call CPU_TIME(t_0)

    allocate(iniGrad(numAtoms * 3) , coords(numAtoms * 3),&
            iniHess(numAtoms * 3, numAtoms * 3))

    !Initializing the variables
    call getIniCoords(coords, seed)
    iniGrad = 0._r8
    iniHess = 0._r8
    ePot = FGH_t(0._r8, iniGrad, iniHess)

    !Calling gradient descent
    call gradDescent(coords, ePot, stepSize, maxIterations, maxLss, verbosicity)
    call CPU_TIME(t)
    print'(a, f10.3)', 'Elapsed time (s): ', t - t_0


contains

    subroutine gradDescent(coords, ePot, stepSize, maxIter, maxLss, verbose)
        type(FGH_t), intent(inout) :: coords(:), ePot
        integer, intent(in) :: maxIter, verbose, maxLss
        real(r8), intent(in) :: stepSize

        integer :: iter, i, eval, lwork, info, evalTot
        character(len=10) :: n_str
        real(r8) :: eDiff, step, ePotOld, stepSize_, iniHess(SIZE(coords), SIZE(coords)), iniGrad(SIZE(coords))
        real(r8) :: grad(SIZE(coords)), r(INT(SIZE(coords)/3), INT(SIZE(coords)/3)), Hess(SIZE(coords), SIZE(coords))
        real(r8) :: r0(INT(SIZE(coords)/3), INT(SIZE(coords)/3)), dir(SIZE(coords)), lambda(SIZE(coords))
        real(r8) :: work2(SIZE(coords) * 3 - 1)

        !Initializing the variables
        lwork = SIZE(coords) *3 -1
        write(n_str, '(i10)') SIZE(coords)
        evalTot = 0
        iniGrad = 0._r8
        iniHess = 0._r8
        grad = 0._r8
        dir = 0._r8
        call calcPotE(coords, ePot)
        eDiff = 1.e300_r8
        call calcDist(coords, r0)

        !Gradient descent algorithm
        do iter = 1, maxIter
            ePotOld = ePot%f
            dir = ePot%g / NORM2(ePot%g)
            stepSize_ = stepSize
            call lss(coords, ePot, dir, stepSize_, maxLss, eval, iniGrad, iniHess)
            evalTot = evalTot + eval
            ePot = FGH_t(0._r8, iniGrad, iniHess)
            call calcPotE(coords, ePot)
            eDiff = ePot%f - ePotOld
            print*,''
            if (verbose >= 2) then
                print*,                   'Iteration:        ', iter
                print'(a,f13.10)',                   'Energydifference: ', eDiff
                if (verbose >= 2) print'(a,f13.10)', 'Energy:           ', ePot%f
            end if
            if (ABS(eDiff) < 1.e-10_r8) then
                if (verbose >= 1) then
                    print*,''
                    print*,'Final Function Value, Gradient, and Hessian:'
                    call ePot%print()
                    call calcDist(coords, r)
                    print*,''
                    print*,'Initial Distance Matrix:'
                    do i = 1, INT(SIZE(coords) / 3)
                        print'('//TRIM(n_str)//'f6.3)', r0(i,:)
                    end do
                    print*,''
                    print*,'Final Distance Matrix:'
                    do i = 1, INT(SIZE(coords) / 3)
                        print'('//TRIM(n_str)//'f6.3)', r(i,:)
                    end do
                    print*,''
                    print*,'Final Coordinates:'
                    do i = 1, INT(SIZE(coords) / 3)
                        print'(i2, f7.3, f7.3, f7.3)', i, coords((i-1)*3+1:(i-1)*3+3)%f
                    end do
                    Hess = ePot%h
                    call DSYEV('V', 'U', SIZE(coords), Hess, SIZE(coords), lambda, work2, lwork, info)
                    print*,''
                    print*, 'Eigenvalues of Hessian:'
                    print'('//TRIM(n_str)//'f10.3)', lambda
                    print*,''
                    print*, 'Iterations:             ', iter
                    print*, 'Line Search evaluations:', evalTot
                    print*,''
                end if
                exit
            end if
        end do
    end subroutine gradDescent

    !Simple inexact golden section line search to speed up the gradient descent algorithm
    subroutine lss(coords, ePot, dir, step, maxLss, eval, iniGrad, iniHess)
        type(FGH_t), intent(inout) :: coords(:), ePot
        real(r8), intent(inout) :: step
        real(r8), intent(in) :: iniGrad(:), iniHess(:,:), dir(:)
        integer, intent(in) :: maxLss
        integer, intent(out) :: eval

        integer  :: i
        real(r8) :: iniF

        iniF = ePot%f

        do i = 1, maxLss
            ePot = FGH_t(0, iniGrad, iniHess)
            coords = coords - step * dir
            call calcPotE(coords, ePot)
            if (ePot%f < iniF) then
                eval = i
                exit
            else
                coords = coords + step * dir
                step = step * 0.618033988749894_r8
            end if
        end do
    end subroutine

    !Initializing initial coordinates at random positions
    subroutine getIniCoords(coords, seed)
        type(FGH_t), intent(inout) :: coords(:)
        integer, intent(in) :: seed

        integer :: n, seeds(8)
        real(r8) :: coordsTemp(SIZE(coords)), iniGrad(SIZE(coords))

        !Generating random numbers
        seeds = seed
        call RANDOM_SEED(put=seeds)
        call RANDOM_NUMBER(coordsTemp)

        !Scaling the coordinates to more sensible distances
        coordsTemp = coordsTemp * 3._r8 - 1.5_r8

        !initializing coords as an FGH type
        !The gradient and the Hessian need to have the right sizes:
        !Here the gradient has the same size as the number of variables (Here the x,y ,and z coordinates of
        !each atom). The gradient needs to have a 1._r8 entry for each coordinate at the respective position,
        !because the derivative of 'x' is 1. The second derivative of 'x' is 0. Thus, the Hessian is initialized
        !as an (number of variables)x(number of variables) matrix filled with 0s.
        !Example for the x and y coordinate of atom 1 of x atoms:
        !
        !coords(1) = FGH_t(2.48_r8, [1,0, 'filled with zeros' ...], '3*numberOfAtoms x 3*numberOfAtoms Matrix
        !filled with zeros')         ^
        !                            |
        !coords(2) = FGH_t(3.94_r8, [0,1, 'filled with zeros' ...], '3*numberOfAtoms x 3*numberOfAtoms Matrix
        !filled with zeros')           ^
        !                              |
        !After initializing the variables our final function depends on the gradient and the Hessian will be
        !automatically evaluated after each mathematical operation.

        iniGrad = 0._r8
        iniHess = 0._r8
        do i = 1, SIZE(coords)
            !Here the gradient is initialized with a 1._r8 in the right place for each variable
            iniGrad(i) = 1._r8
            coords(i) = FGH_t(coordsTemp(i), iniGrad, iniHess)
            iniGrad = 0._r8
        end do
    end subroutine getIniCoords

    !Calculation of the distance matrix
    subroutine calcDist(coords, r)
        type(FGH_t), intent(in) :: coords(:)
        real(r8), intent(out) :: r(:,:)
        integer :: i, j

        do i = 1, SIZE(r(1,:))
            do j = 1, SIZE(r(1,:))
                r(j, i) = SQRT((coords(1 + 3 * (i-1))%f - coords(1 + 3 * (j-1))%f) ** 2 + &
                        (coords(2 + 3 * (i-1))%f - coords(2 + 3 * (j-1))%f) ** 2 + &
                        (coords(3 + 3 * (i-1))%f - coords(3 + 3 * (j-1))%f) ** 2)
            end do
        end do
    end subroutine

    !Calculating pair-wise Lennard-Jones potential energy, due to using automatic differentiation the gradient
    !and the Hessian is also calculated!
    pure subroutine calcPotE(coords, ePot)
        implicit none
        type(FGH_t), intent(in) :: coords(:)
        type(FGH_t), intent(inout) :: ePot
        integer :: i, j, numAtoms
        type(FGH_t) :: r

        numAtoms = INT(SIZE(coords) / 3)

        do i = 1, numAtoms
            do j = i + 1, numAtoms
                !All the operators have been overwritten to work with the FGH type, therefore formulas can be
                !implemented normally!
                r = SQRT((coords(1 + 3 * (i-1)) - coords(1 + 3 * (j-1))) ** 2 + &
                (coords(2 + 3 * (i-1)) - coords(2 + 3 * (j-1))) ** 2 + &
                (coords(3 + 3 * (i-1)) - coords(3 + 3 * (j-1))) ** 2)
                ePot = ePot + 4._r8 * (r**(-12) - r**(-6))
            end do
        end do
    end subroutine calcPotE
end program example