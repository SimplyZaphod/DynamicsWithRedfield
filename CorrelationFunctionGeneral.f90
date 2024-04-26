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
    call ReadInput(inputsystem, dimension, rho, H, adag, ahat, op1, op2, nameout,verbose='v')
    write(*,*)

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Let's read the input for evolution
    !!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Input name:'
    read(*,'(A100)') inputevol
    call ReadTimeEvolutionInput(inputevol, deltat, timesteps, n_dt, DynType, T, densitytype, gammatype, omega_g)
    write(*,*)
    write(*,*) 'ACHTUNG! For now the unitary evolution only is allowed! Take care'
    write(*,*) 'Write the generating function'
    call commutator_complex(rho, op2, rho, dimension)
    rho = -imagine*rho
!    call npmatmul_complex(rho, op2, rho, dimension, dimension, dimension)

    write(*,*) 'Select the dynamics!'
    call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, nameout)

    !if((DynType=='runge-kutta-unitary').or.(DynType=='RKU'))then
!        call CorrelationWithRungeKuttaUnitary(op1,rho, H, dimension, Deltat, timesteps, nameout)
!    elseif ((DynType=='arnoldi').or.(DynType=='arnoldi'))then
!        write(*,*) 'No krylov dimension is inserted! 10 is taken as a guess'
!        call EvolveWithArnoldi(rho, LiouvilleTensor, dimension, 10, Deltat, timesteps,'*')
!    elseif ((DynType=='D').or.(DynType=='diag'))then
!        call EvolveWithLiouvilleDiagonalization(rho, LiouvilleTensor, dimension, Deltat, timesteps, '*')
!    end if
    write(*,*) "After the evolution"
    call write_matrix_complex(rho, dimension, dimension, 'e9.2')
    write(*,*) 'End program!'
contains
    subroutine ReadInput(name, dim, rho, H, a1, a2,op1,op2,nameout, verbose)
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
        !> first doublecomplex(dim, dim) operators that does not evolve with time as a
        !> generating function
        doublecomplex, allocatable, intent(inout) :: op1(:,:)
        !> first doublecomplex(dim, dim) operators that does evolve with time as a
        !> generating function
        doublecomplex, allocatable, intent(inout) :: op2(:,:)
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

end program RedfieldCorr
