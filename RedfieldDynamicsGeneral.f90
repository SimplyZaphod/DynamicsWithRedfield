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
    integer dimension, n_prop
    doublecomplex, allocatable :: H(:,:), adag(:,:), ahat(:,:), eigval(:), rho(:,:), prop(:,:,:)
    doubleprecision, allocatable :: eigenvalues(:)
    character(len=100) inputsystem
    !!!!!!!!!!!!!
    ! Variable for time evolution
    !!!!!!!!!!!!!
    character(len=100) inputevol
    doubleprecision Deltat,t_0,t_f, T
    integer timesteps, n_dt
    character(len=100) DynType, DensityType,GammaType
    doubleprecision omega_g, eta_g

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
    doubleprecision waste
    doublecomplex,allocatable :: WMComplex(:,:)
    !!!!!!!!!!!!!!!!!!!
    !Diagonalization variable
    !!!!!!!!!!!!!!!!!!!
    doublecomplex, allocatable :: work(:)
    doubleprecision, allocatable :: rwork(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Let's read the input
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Input name:'
    read(*,'(A100)') inputsystem
    call ReadInput(inputsystem, dimension, rho, H, adag, ahat,n_prop, prop)
    write(*,*)
    write(*,*) 'Hamiltonian:'
    call write_matrix_complex(H, dimension, dimension, 'e9.2')
    write(*,*) 'Diagonalize? Check the Hamiltonian...'
    do i=1, dimension
        do j=1, dimension
            if (i.ne.j) waste = waste + dreal(conjg(H(i,j))*H(i,j))
        end do
    end do
    waste = dsqrt(waste)/((dimension-1)*dimension)
    if (waste > 1.d-10) then
        write(*,*) 'Not diagonal! Diagonalize...'
        allocate(WMComplex(dimension, dimension), eigenvalues(dimension))
        allocate(work(3*dimension), rwork(3*dimension-2))
        call zheev('N', 'V', dimension, WMComplex, dimension, eigenvalues,work, 3*dimension, rwork, i)
        write(*,*) 'Diagonalized! Now I rotate everything'
        call rotate_matrix_complex(H, WMComplex, dimension)
        call rotate_matrix_complex(rho, WMComplex, dimension)
        call rotate_matrix_complex(adag, WMComplex, dimension)
        call rotate_matrix_complex(ahat, WMComplex, dimension)
        do i=1, n_prop
            call rotate_matrix_complex(prop(i,:,:), WMComplex, dimension)
        end do
        write(*,*) 'End! Deallocate'
        deallocate(work, rwork, WMComplex, eigenvalues)
        write(*,*) 'New Hamiltonian:'
        call write_matrix_complex(H, dimension, dimension, 'e9.2')
    else
        write(*,*) 'H is diagonal! Every possible error is related to the input file and not on the diagonalization...'
    end if

    write(*,*) 'Density:'
    call write_matrix_complex(rho, dimension, dimension, 'e9.2')
    write(*,*) 'Adag:'
    call write_matrix_complex(adag, dimension, dimension, 'e9.2')
    write(*,*) 'ahat:'
    call write_matrix_complex(ahat, dimension, dimension, 'e9.2')
    do i=1, n_prop
        write(*,*) 'Property', i
        call write_matrix_complex(prop(i,:,:), dimension, dimension,'e9.2')
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Let's read the input for evolution
    !!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Input name:'
    read(*,'(A100)') inputevol
    call ReadTimeEvolutionInput(inputevol, t_0,t_f, timesteps, n_dt, DynType, T, densitytype, gammatype, omega_g, eta_g)
    write(*,*)
    deltat = abs(t_f-t_0)/timesteps
    !!!!!!!!!!!!!!
    ! Fill gamma and density
    !!!!!!!!!!!!!!
    allocate(eigenvalues(dimension))
    do i=1, dimension
        eigenvalues(i) = H(i,i)
    end do

    allocate(gammas(dimension, dimension), density(dimension, dimension))
    call FilLGamma(dimension, Gammas, H, GammaType, omega_g, eta_g)
    write(*,*) 'Gamma filled!'
!    call write_matrix_double(real(gammas), dimension, dimension, 'e9.2')
    call FilLDensity(dimension, density, H, DensityType, T)
    write(*,*) 'Density filled!'
!    call write_matrix_double(density, dimension, dimension, 'e9.2')

    !!!!!!
    ! write the redfield tensor if it is necessary
    !!!!!!
    allocate(LiouvilleTensor(dimension, dimension, dimension,dimension))
    LiouvilleTensor = (0.d0, 0.d0)
    write(*,*) "Before the evolution"
    call write_matrix_complex(rho, dimension, dimension, 'e9.2')
    write(*,*) 'Initial population'
    write(*,'(<dimension>(F6.3, X))') (real(rho(i,i)), i=1, dimension)
    write(*,*) 'Select the dynamics!'
    if((DynType=='runge-kutta-unitary').or.(DynType=='RKU'))then
        call EvolveWithRungeKuttaUnitary(rho, H, dimension, Deltat, timesteps,n_prop, prop, '*')
    elseif((DynType=='runge-kutta-redfield').or.(DynType=='RKR'))then
        write(*,*) 'Writing redfield tensor'
        LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
        call EvolveWithRungeKuttaRedfield(rho, H, LiouvilleTensor, dimension, Deltat, timesteps, n_prop, prop, '*')
    elseif ((DynType=='arnoldi').or.(DynType=='A'))then
        write(*,*) 'No krylov dimension is inserted! 20 is taken as a guess'
        write(*,*) 'Writing redfield tensor'
        LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
        LiouvilleTensor = Liouvillian_Complex(H, LiouvilleTensor, dimension)
        call EvolveWithArnoldi(rho, LiouvilleTensor, dimension, 20, Deltat, timesteps, n_prop, prop, '*')
    elseif ((DynType=='D').or.(DynType=='diag'))then
        write(*,*) 'Writing redfield tensor'
        LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
        LiouvilleTensor = Liouvillian_Complex(H, LiouvilleTensor, dimension)
        call EvolveWithLiouvilleDiagonalization(rho, LiouvilleTensor, dimension, Deltat, timesteps, n_prop, prop,'*')
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
    subroutine ReadInput(name, dim, rho, H, a1, a2,n_prop, prop)
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
        !> integer number of properties to be calculated
        integer, intent(inout) :: n_prop
        !> doublecomplex(dim, dim) with properties to be calculated
        doublecomplex, allocatable, intent(inout) :: prop(:,:,:)

        character(len=50) rhoname, Hname, a1name, a2name, nameprop
        logical verb
        integer i

        write(*,*) "Read from:", name
        open(1, file=name)
        read(1,*) dim
        read(1,*) rhoname
        read(1,*) Hname
        read(1,*) a1name
        read(1,*) a2name
        read(1,*) n_prop
        if(n_prop.gt.0)then
            allocate(prop(n_prop, dim, dim))
            do i=1, n_prop
                nameprop = ''
                read(1,*) nameprop
                write(*,*) 'Property read from: ', nameprop
                open(2, file=nameprop, form='unformatted',access='stream')
                read(2) prop(i,:,:)
                close(2)
            end do
        else
            allocate(prop(0, dim, dim))
        end if

        close(1)
        allocate(rho(dim,dim), H(dim,dim), a1(dim, dim), a2(dim, dim))

        open(1,file=rhoname, form='unformatted', access='stream')
        read(1) rho(1:,1:)
        close(1)
        open(1,file=Hname, form='unformatted', access='stream')
        read(1) H(1:,1:)
        close(1)
        open(1,file=a1name, form='unformatted', access='stream')
        read(1) a1(1:,1:)
        close(1)
        open(1,file=a2name, form='unformatted', access='stream')
        read(1) a2(1:,1:)
        close(1)

        write(*,*) 'rho read from: ', rhoname
        write(*,*) 'H read from: ', Hname
        write(*,*) 'first operator read from: ',a1name
        write(*,*) 'second operator read from: ', a2name

    end subroutine ReadInput

end program RedfieldDA
