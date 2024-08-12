program RedfieldCorr
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
    doublecomplex, allocatable :: H(:,:), adag(:,:), ahat(:,:), eigval(:), rho(:,:),&
            op1(:,:), op2(:,:)
    doubleprecision, allocatable :: eigenvalues(:)
    character(len=100) inputsystem, nameout
    character(len=150) spettro !!! just in case some trash is needed
    !!!!!!!!!!!!!
    ! Variable for time evolution
    !!!!!!!!!!!!!
    character(len=100) inputevol
    doubleprecision Deltat,t_0,t_f,T
    integer timesteps, n_dt
    character(len=100) DynType, DensityType,GammaType
    integer selGF
    doubleprecision omega_g, eta_g

    !!!!!!!!!!!!!
    !Time evolution
    !!!!!!!!!!!!!!
    doublecomplex, allocatable :: LiouvilleTensor(:,:,:,:), gammas(:,:)
    doubleprecision, allocatable :: density(:,:)
    integer nkrylov
    !!!!!!!!!!!!!!!
    !Time reversal
    !!!!!!!!!!!!!!!
    logical symm
    doubleprecision, allocatable, dimension(:,:) :: vector1, vector2
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
    call ReadInput(inputsystem, dimension, rho, H, adag, ahat, op1, op2, selGF,nameout,verbose='v')
    write(*,*)
    allocate(eigenvalues(dimension))
    do i=1, dimension
        eigenvalues(i) = real(H(i,i))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Let's read the input for evolution
    !!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Input name:'
    read(*,'(A100)') inputevol
    call ReadTimeEvolutionInput(inputevol, t_0,t_f, timesteps, n_dt, DynType, T, densitytype, gammatype, omega_g, eta_g)
    write(*,*)
    write(*,*) 'ACHTUNG! For now the unitary evolution only is allowed! Take care'
    write(*,*) 'Write the generating function!'
    write(*,*) 'Selected GF:', selGF
    if(selGF.eq.1)then
        write(*,*) 'GF = -i[op2, rho]'
        call CorrelationFunction1(rho, op2, rho, dimension)
        op2 = rho
    elseif(selGF.eq.2)then
        write(*,*) 'GF = i*op2*rho'
        call CorrelationFunction2(rho, op2, rho, dimension)
        op2 = rho
    elseif(selGF.eq.3)then
        write(*,*) 'GF = i*rho*op2'
        call CorrelationFunction3(rho, op2, rho, dimension)
        op2 = rho
    end if
    !!!!!!!!!!!!!!!!!!!!!
    !In this subroutine I need to select the time such that:
    !-if t_0*t_f>0, then everything is shifted such that t_0 = 0
    !-if t_0*t_f<0, then everything is adjusted such that the 0 is in the middle of the interval
    !!!!!!!!!!!!!!!!!!!!!
    if((t_0*t_f).ge.0.d0)then
        t_f = t_f-t_0
        t_0 = 0.d0
        symm = .false.
    else
        if(mod(timesteps,2))then
            timesteps = timesteps/2 +1
            write(*,*) "Timestep is changed! A point is added"
        else
            timesteps = timesteps + 1
            timesteps = timesteps/2
        end if
        t_f = (t_f-t_0)/2.d0
        t_0 = 0.d0
        symm = .true.
    end if

    if(.not.symm)then
        write(*,*) 'Select the dynamics!'
        Deltat = t_f/timesteps
        if((DynType=='runge-kutta-unitary').or.(DynType=='RKU'))then
            call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, nameout)
        elseif((DynType=='runge-kutta-redfield').or.(DynType=='RKR'))then
            write(*,*) 'Writing redfield tensor'
            allocate(gammas(dimension, dimension), density(dimension, dimension), LiouvilleTensor(dimension, dimension, dimension, dimension))
            call FilLGamma(dimension, Gammas, H, GammaType, omega_g, eta_g)
            call FilLDensity(dimension, density, H, DensityType, T)
            LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
            call CorrelationWithRungeKuttaRedfield(op1,rho, H,LiouvilleTensor, dimension, Deltat, timesteps, nameout)
            deallocate(gammas, density, LiouvilleTensor)
        elseif ((DynType=='arnoldi').or.(DynType=='A'))then
            write(*,*) 'No krylov dimension is inserted! 20 is taken as a guess'
            write(*,*) 'Writing redfield tensor'
            allocate(gammas(dimension, dimension), density(dimension, dimension), LiouvilleTensor(dimension, dimension, dimension, dimension))
            call FilLGamma(dimension, Gammas, H, GammaType, omega_g, eta_g)
            call FilLDensity(dimension, density, H, DensityType, T)
            LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
            LiouvilleTensor = Liouvillian_Complex(H, LiouvilleTensor, dimension)
            call CorrelationWithArnoldi(op1,rho, LiouvilleTensor, dimension,20, Deltat, timesteps, nameout)
            deallocate(gammas, density, LiouvilleTensor)
        elseif ((DynType=='D').or.(DynType=='diag'))then
            write(*,*) 'Writing redfield tensor'
            write(*,*) 'This DOES NOT WORK!'
            allocate(gammas(dimension, dimension), density(dimension, dimension), LiouvilleTensor(dimension, dimension, dimension, dimension))
            call FilLGamma(dimension, Gammas, H, GammaType, omega_g, eta_g)
            call FilLDensity(dimension, density, H, DensityType, T)
            LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
            LiouvilleTensor = Liouvillian_Complex(H, LiouvilleTensor, dimension)
            call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, nameout)
            deallocate(gammas, density, LiouvilleTensor)
        end if
        write(*,*) "After the evolution"
        call write_matrix_complex(rho, dimension, dimension, 'e9.2')
        write(*,*) 'End program!'
    else
        write(*,*) 'The evolution is taken out in two steps'
        write(*,*) 'First step: from t_0=0 to t_f=t_f/2'
        spettro = ''
        write(spettro, *) nameout,'.1.tmp'
        call StripSpaces(spettro)
        Deltat = t_f/timesteps
        if((DynType=='runge-kutta-unitary').or.(DynType=='RKU'))then
            call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, spettro)
        elseif((DynType=='runge-kutta-redfield').or.(DynType=='RKR'))then
            write(*,*) 'Writing redfield tensor'
            allocate(gammas(dimension, dimension), density(dimension, dimension), LiouvilleTensor(dimension, dimension, dimension, dimension))
            call FilLGamma(dimension, Gammas, H, GammaType, omega_g, eta_g)
            call FilLDensity(dimension, density, H, DensityType, T)
            LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
            call CorrelationWithRungeKuttaRedfield(op1,rho, H,LiouvilleTensor, dimension, Deltat, timesteps, spettro)
        elseif ((DynType=='arnoldi').or.(DynType=='A'))then
            write(*,*) 'No krylov dimension is inserted! 20 is taken as a guess'
            write(*,*) 'Writing redfield tensor'
            allocate(gammas(dimension, dimension), density(dimension, dimension), LiouvilleTensor(dimension, dimension, dimension, dimension))
            call FilLGamma(dimension, Gammas, H, GammaType, omega_g, eta_g)
            call FilLDensity(dimension, density, H, DensityType, T)
            LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
            LiouvilleTensor = Liouvillian_Complex(H, LiouvilleTensor, dimension)
            call CorrelationWithArnoldi(op1,rho, LiouvilleTensor, dimension,20, Deltat, timesteps, spettro)
        elseif ((DynType=='D').or.(DynType=='diag'))then
            write(*,*) 'Writing redfield tensor'
            write(*,*) 'This DOES NOT WORK!'
            allocate(gammas(dimension, dimension), density(dimension, dimension), LiouvilleTensor(dimension, dimension, dimension, dimension))
            call FilLGamma(dimension, Gammas, H, GammaType, omega_g, eta_g)
            call FilLDensity(dimension, density, H, DensityType, T)
            LiouvilleTensor = RedfieldTensor_Complex(adag, ahat, gammas, density, eigenvalues, dimension)
            LiouvilleTensor = Liouvillian_Complex(H, LiouvilleTensor, dimension)
            call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, spettro)
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'END OF FIRST STEP'
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'Second step: from t_0=0 to t_f=-t_f/2'
        rho = op2
        spettro = ''
        write(spettro, *) nameout,'.2.tmp'
        call StripSpaces(spettro)
        Deltat = -Deltat
        if((DynType=='runge-kutta-unitary').or.(DynType=='RKU'))then
            call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, spettro)
        elseif((DynType=='runge-kutta-redfield').or.(DynType=='RKR'))then
            call CorrelationWithRungeKuttaRedfield(op1,rho, H,LiouvilleTensor, dimension, Deltat, timesteps, spettro)
            deallocate(gammas, density, LiouvilleTensor)
        elseif ((DynType=='arnoldi').or.(DynType=='A'))then
            call CorrelationWithArnoldi(op1,rho, LiouvilleTensor, dimension,20, Deltat, timesteps, spettro)
            deallocate(gammas, density, LiouvilleTensor)
        elseif ((DynType=='D').or.(DynType=='diag'))then
            call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, spettro)
            deallocate(gammas, density, LiouvilleTensor)
        end if
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'END OF FIRST STEP'
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'Now it is necessary to reorder the spectrum'
        if(GammaType.eq.'exp') write(*,*) 'Additionally I will damp with an exponential with damping 0.05 fs'
        allocate(vector2(timesteps*2-1,3))
        spettro = ''
        write(spettro, *) nameout,'.1.tmp'
        call StripSpaces(spettro)
        open(2, file=spettro)
        do i=1, timesteps
            read(2, *) vector2(i+timesteps-1,1), vector2(i+timesteps-1,2), vector2(i+timesteps-1,3)
        end do
        close(2)

        spettro = ''
        write(spettro, *) nameout,'.2.tmp'
        call StripSpaces(spettro)
        open(2, file=spettro)
        open(2, file=spettro)
        read(2,*)
        do i=1, timesteps-1
            read(2, *) vector2(timesteps-i,1), vector2(timesteps-i,2), vector2(timesteps-i,3)
        end do
        close(2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        !Now reorder
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(GammaType.eq.'exp')then
            open(1, file=nameout)
            do i=1, 2*timesteps-1
                write(1,'(F, x, 2(e18.6E5, x))') vector2(i,1), vector2(i,2)*exp(-0.05d0*abs(vector2(i,1))), vector2(i,3)*exp(-0.05d0*abs(vector2(i,1)))
            end do
            close(1)
        else
            open(1, file=nameout)
            do i=1, 2*timesteps-1
                write(1,'(F, x, 2(e18.6E5, x))') vector2(i,1), vector2(i,2), vector2(i,3)
            end do
            close(1)
        end if
    end if
contains
    subroutine ReadInput(name, dim, rho, H, a1, a2,op1,op2,GF, nameout, verbose)
        implicit none
        !> character name of the input file
        character(len=*), intent(in) :: name
        !> integer dimension of the system
        integer, intent(inout) :: dim
        !> doublecomplex(dim,dim) initial density matrix
        doublecomplex, allocatable, intent(inout) :: rho(:,:)
        !> doublecomplex (dim, dim) Hamiltonian
        doublecomplex, allocatable, intent(inout) :: H(:,:)
        !> first doublecomplex(dim, dim) operators coupled w    ith the Bath
        doublecomplex, allocatable, intent(inout) :: a1(:,:)
        !> second doublecomplex(dim, dim) operators coupled with the Bath
        doublecomplex, allocatable, intent(inout) :: a2(:,:)
        !> first doublecomplex(dim, dim) operators that does not evolve with time as a
        !> generating function
        doublecomplex, allocatable, intent(inout) :: op1(:,:)
        !> first doublecomplex(dim, dim) operators that does evolve with time as a
        !> generating function
        doublecomplex, allocatable, intent(inout) :: op2(:,:)
        !>integer type of GF to be used
        integer, intent(inout) :: GF
        !> name of the output necessary
        character(len=100) nameout
        !> character optional, if "verb", "verbose" or "v" I will print everything
        character(len=*), optional :: verbose

        character(len=50) rhoname, Hname, a1name, a2name, op1name, op2name
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
        read(1,*) op1name
        read(1,*) op2name
        read(1, *) GF
        read(1,*) nameout
        close(1)
        allocate(rho(dim,dim), H(dim,dim), a1(dim, dim), a2(dim, dim),&
        op1(dim,dim), op2(dim,dim))

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
        open(1,file=op1name, form='unformatted')
        read(1) op1(1:,1:)
        close(1)
        open(1,file=op2name, form='unformatted')
        read(1) op2(1:,1:)
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
            write(*,*) 'second operator read from: ', op1name
            call write_matrix_complex(op1, dim, dim)
            write(*,*) 'second operator read from: ', op2name
            call write_matrix_complex(op2, dim, dim)
        else
            write(*,*) 'rho read from: ', rhoname
            write(*,*) 'H read from: ', Hname
            write(*,*) 'first operator read from: ',a1name
            write(*,*) 'second operator read from: ', a2name
            write(*,*) 'first operator read from: ',op1name
            write(*,*) 'second operator read from: ', op2name
        end if
    end subroutine ReadInput
    !This function return the starting correlation function as Omega=-i[op, rho]
    subroutine CorrelationFunction1(GF, op, rho, dim)
        !dimension of the system
        integer, intent(in) :: dim
        !output correlation function
        doublecomplex, intent(inout) :: GF(dim, dim)
        !doublecomplex (dim,dim) operator
        doublecomplex, intent(in) :: op(dim, dim)
        !doublecomplex density matrix
        doublecomplex, intent(in) :: rho(dim,dim)

        doublecomplex, dimension(:,:), allocatable :: aux

        allocate(aux(dim, dim))
        call commutator_complex(aux, op, rho, dim)
        aux = -imagine*aux
        GF = aux
    end subroutine CorrelationFunction1
    !This function return the starting correlation function as Omega= op*rho
    subroutine CorrelationFunction2(GF, op, rho, dim)
        !dimension of the system
        integer, intent(in) :: dim
        !output correlation function
        doublecomplex, intent(inout) :: GF(dim, dim)
        !doublecomplex (dim,dim) operator
        doublecomplex, intent(in) :: op(dim, dim)
        !doublecomplex density matrix
        doublecomplex, intent(in) :: rho(dim,dim)

        doublecomplex, dimension(:,:), allocatable :: aux

        allocate(aux(dim, dim))
        call npmatmul_complex(aux, op, rho, dim,dim,dim)
        GF = imagine*aux
    end subroutine CorrelationFunction2
    !This function return the starting correlation function as Omega= rho*op
    subroutine CorrelationFunction3(GF, op, rho, dim)
        !dimension of the system
        integer, intent(in) :: dim
        !output correlation function
        doublecomplex, intent(inout) :: GF(dim, dim)
        !doublecomplex (dim,dim) operator
        doublecomplex, intent(in) :: op(dim, dim)
        !doublecomplex density matrix
        doublecomplex, intent(in) :: rho(dim,dim)

        doublecomplex, dimension(:,:), allocatable :: aux

        allocate(aux(dim, dim))
        call npmatmul_complex(aux, rho, op, dim,dim,dim)
        GF = -imagine*aux
    end subroutine CorrelationFunction3
end program RedfieldCorr
