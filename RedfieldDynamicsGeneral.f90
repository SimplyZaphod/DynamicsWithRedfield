program RedfieldDA
    use conversion
    use Mathfunc
    use Redfield
    use UsefullIO
    use Evolutions


    implicit none
    !!!!!!!!!!!!!
    !Variables for the system
    !!!!!!!!!!!!!system
    integer dimension
    doublecomplex, allocatable :: H(:,:), adag(:,:), ahat(:,:), eigval(:), rho(:,:)
    doubleprecision, allocatable :: eigenvalues(:)
    character(len=100) inputsystem
    !!!!!!!!!!!!!
    ! Variable for time evolution
    !!!!!!!!!!!!!
    character(len=100) inputevol
    doubleprecision Deltat, T
    integer timesteps, n_dt
    character(len=100) DynType, DensityType,GammaType
    doubleprecision omega_g

    !!!!!!!!!!!!!
    !Time evolution
    !!!!!!!!!!!!!!
    doublecomplex, allocatable :: LiouvilleTensor(:,:,:,:), gammas(:,:)
    doubleprecision, allocatable :: density(:,:)
    integer nkrylov
    !!!!!!!!!!!!!!!!
    ! Save the final matrix
    !!!!!!!!!!!!!!!!
    character(len=100) name

    !!!!!!!!!!!!!!!!
    !All the necessary for cycles and weste
    !!!!!!!!!!!!!!!!
    integer i,j

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Let's read the input
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Input name:'
    read(*,'(A100)') inputsystem
    call ReadInput(inputsystem, dimension, rho, H, adag, ahat, verbose='v')
    write(*,*)

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Let's read the input for evolution
    !!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Input name:'
    read(*,'(A100)') inputevol
    call ReadTimeEvolutionInput(inputevol, deltat, timesteps, n_dt, DynType, T, densitytype, gammatype, omega_g)
    write(*,*)

    !!!!!!!!!!!!!!
    ! Fill gamma and density
    !!!!!!!!!!!!!!
    allocate(eigenvalues(dimension))
    do i=1, dimension
        eigenvalues(i) = H(i,i)
    end do

    allocate(gammas(dimension, dimension), density(dimension, dimension))
    call FilLGamma(dimension, Gammas, H, GammaType, omega_g)
    write(*,*) 'Gamma'
    call write_matrix_double(real(gammas), dimension, dimension, 'e9.2')
    call FilLDensity(dimension, density, H, DensityType, T)
    call write_matrix_double(density, dimension, dimension, 'e9.2')

    !!!!!!
    ! write the redfield tensor if it is necessary
    !!!!!!
    allocate(LiouvilleTensor(dimension, dimension, dimension,dimension))
    LiouvilleTensor = (0.d0, 0.d0)
    if((DynType=='runge-kutta-unitary').or.(DynType=='RKU'))then
        write(*,*) 'Unitary dynamics! No redfield tensor written'
    else
        write(*,*) 'Writing redfield tensor'
        LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
        call CheckDetailBalanceRedfield_Double(real(LiouvilleTensor), eigenvalues, T,dimension)
        LiouvilleTensor = Liouvillian_Complex(H/(1.d15), LiouvilleTensor, dimension)
        write(*,*) 'Liouville tensor written'
    end if
    write(*,*) "Before the evolution"
    call write_matrix_complex(rho, dimension, dimension, 'e9.2')
    write(*,*) 'Initial population'
    write(*,'(<dimension>(F6.3, X))') (real(rho(i,i)), i=1, dimension)
    write(*,*) 'Select the dynamics!'
    if((DynType=='runge-kutta-unitary').or.(DynType=='RKU'))then
        call EvolveWithRungeKuttaUnitary(rho, H, dimension, Deltat, timesteps, '*')
    elseif((DynType=='runge-kutta').or.(DynType=='RK'))then
        call EvolveWithRungeKutta(rho, LiouvilleTensor, dimension, Deltat, timesteps, '*')
    elseif ((DynType=='arnoldi').or.(DynType=='arnoldi'))then
        write(*,*) 'No krylov dimension is inserted! 10 is taken as a guess'
        call EvolveWithArnoldi(rho, LiouvilleTensor, dimension, 10, Deltat, timesteps,'*')
    elseif ((DynType=='D').or.(DynType=='diag'))then
        call EvolveWithLiouvilleDiagonalization(rho, LiouvilleTensor, dimension, Deltat, timesteps, '*')
    end if
    write(*,*) "After the evolution"
    call write_matrix_complex(rho, dimension, dimension, 'e9.2')
    write(*,*) 'Final population'
    write(*,'(<dimension>(F6.3, X))') (real(rho(i,i)), i=1, dimension)
    open(1, file='rho_final.dat', form='unformatted')
    write(1) rho
    close(1)
    write(*,*) 'End program!'
contains
    subroutine ReadInput(name, dim, rho, H, a1, a2, verbose)
        implicit none
        !> character name of the input file
        character(len=*), intent(in) :: name
        !> integer dimension of the system
        integer, intent(inout) :: dim
        !> doublecomplex(dim,dim) initial density matrix
        doublecomplex, allocatable, intent(inout) :: rho(:,:)
        !> doublecomplex (dim, dim) Hamiltonian
        doublecomplex, allocatable, intent(inout) :: H(:,:)
        !> first doublecomplex(dim, dim) operators coupled with the Bath
        doublecomplex, allocatable, intent(inout) :: a1(:,:)
        !> second doublecomplex(dim, dim) operators coupled with the Bath
        doublecomplex, allocatable, intent(inout) :: a2(:,:)
        !> character optional, if "verb", "verbose" or "v" I will print everything
        character(len=*), optional :: verbose

        character(len=50) rhoname, Hname, a1name, a2name
        logical verb

        verb=.false.
        if(present(verbose))then
            if((verbose.eq."v").or.(verbose.eq."verb").or.(verbose.eq."verbose").or.(verbose.eq."yes").or.(verbose.eq.".true."))then
                verb = .true.
            end if
        end if
        write(*,*) "Read from:", name
        open(1, file=name)
        read(1,*) dim
        read(1,*) rhoname
        read(1,*) Hname
        read(1,*) a1name
        read(1,*) a2name
        close(1)
        allocate(rho(dim,dim), H(dim,dim), a1(dim, dim), a2(dim, dim))

        open(1,file=rhoname, form='unformatted')
        read(1) rho(1:,1:)
        close(1)
        open(1,file=Hname, form='unformatted')
        read(1) H(1:,1:)
        close(1)
        open(1,file=a1name, form='unformatted')
        read(1) a1(1:,1:)
        close(1)
        open(1,file=a2name, form='unformatted')
        read(1) a2(1:,1:)
        close(1)

        if(verb)then
            write(*,*) 'rho read from: ', rhoname
            call write_matrix_complex(rho, dim, dim)
            write(*,*) 'H read from: ', Hname
            call write_matrix_complex(H, dim, dim)
            write(*,*) 'first operator read from: ',a1name
            call write_matrix_complex(a1, dim, dim)
            write(*,*) 'second operator read from: ', a2name
            call write_matrix_complex(a2, dim, dim)
        else
            write(*,*) 'rho read from: ', rhoname
            write(*,*) 'H read from: ', Hname
            write(*,*) 'first operator read from: ',a1name
            write(*,*) 'second operator read from: ', a2name
        end if
    end subroutine ReadInput

end program RedfieldDA
