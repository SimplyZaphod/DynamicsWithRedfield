module Evolutions
    use conversion
    use Mathfunc
    use Redfield
    use UsefullIO
    contains

    !> @todo check if I imposed the optional constant in a correct way
    subroutine FilLGamma(dim, Gammas, H, GammaType, const)
        !> integer dimension of Gamma matrix
        integer, intent(in) :: dim
        !> doublecomplex (dim, dim) Gamma matrix
        doublecomplex, intent(out) :: Gammas(dim,dim)
        !> doublecomplex (dim, dim) hamiltonian
        doublecomplex, intent(in) :: H(dim,dim)
        !> character type of gamma
        character(len=*), intent(in) :: GammaType
        !> doubleprecision optional constant in case of constant density
        doubleprecision,optional, intent(in) :: const

        integer i
        doubleprecision, allocatable :: values(:), gammatoken(:,:)
        write(*,*) 'Write the Gamma Matrix'
        allocate(values(dim), gammatoken(dim,dim))
        do i=1, dim
            values(i) = real(H(i,i))
        end do
        Gammas = (0.d0, 0.d0)
        if((GammaType=='constant').or.(GammaType=='const'))then
            write(*,*) 'Constant spectral density'
            call GammaConstant(gammatoken, const, values, dim)
        elseif((GammaType=='emission').or.(GammaType=='radiation'))then
            write(*,*) 'Radiative spectral density'
            call GammaRadiationField(gammatoken, values, dim)
        elseif((GammaType=='null').or.(GammaType=='zero'))then
            write(*,*) 'Null spectral density'
            call GammaConstant(gammatoken, 0.d0, values, dim)
        endif
        Gammas = gammatoken*reality
    end subroutine FilLGamma

    subroutine FilLDensity(dim, density, H, DensityType, Temperature)
        !> integer dimension of Gamma matrix
        integer, intent(in) :: dim
        !> doubleprecision (dim) density matrix
        doubleprecision, intent(out) :: density(dim,dim)
        !> doublecomplex (dim) hamiltonian matrix
        doublecomplex, intent(in) :: H(dim,dim)
        !> character requested density
        character(len=*), intent(in) :: DensityType
        !> doubleprecision temperature
        doubleprecision, intent(in) :: Temperature

        integer i, j
        doubleprecision omega
        write(*,*) 'Fill the bath density matrix!'
        write(*,*) 'Density type: ', DensityType

        density = 0.d0
        do i=1, dim
            do j=1, dim
                omega = real(H(i,i)-H(j,j))
                density(i,j) = choose_density(DensityType, omega, Temperature)
            end do
        end do

    end subroutine FilLDensity

    subroutine EvolveWithRungeKuttaUnitary(rho, H, dim, Deltat, timesteps, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim) Hamiltonian
        doublecomplex, intent(in) :: H(dim,dim)
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name
        write(*,*) 'Start evolution with Runge Kutta'
        if(present(name))then
            if(name.ne."*")then
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call Runge_Kutta_density(rho, H, dim, Deltat)
                end do
                write(1) rho(1:,1:)
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F,A,F)") 'Time: ', Deltat*i, 'trace=', real(trace_complex(rho, dim))
                    call write_matrix_complex(rho, dim, dim, 'e11.3')
                    call Runge_Kutta_density(rho, H, dim, Deltat)
                end do
                write(*,"(A, F)") 'Time: ', Deltat*(i+1)
                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call Runge_Kutta_density(rho, H, dim, Deltat)
            end do
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
    end subroutine EvolveWithRungeKuttaUnitary

    subroutine EvolveWithRungeKutta(rho, liouv, dim, Deltat, timesteps, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim) Hamiltonian
        doublecomplex, intent(in) :: liouv(dim,dim, dim, dim)
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name
        write(*,*) 'Start evolution with Runge Kutta'
        if(present(name))then
            if(name.ne."*")then
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call Runge_Kutta_Liouville(rho, liouv, dim, Deltat)
                end do
                write(1) rho(1:,1:)
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F,A,F)") 'Time: ', Deltat*i, 'trace=', real(trace_complex(rho, dim))
                    call write_matrix_complex(rho, dim, dim, 'e11.3')
                    call Runge_Kutta_Liouville(rho, liouv, dim, Deltat)
                end do
                write(*,"(A, F)") 'Time: ', Deltat*(i+1)
                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call Runge_Kutta_Liouville(rho, liouv, dim, Deltat)
            end do
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
    end subroutine EvolveWithRungeKutta

    subroutine EvolveWithLiouvilleDiagonalization(rho, Liouv, dim, Deltat, timesteps, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim, dim, dim) Liouvillian
        doublecomplex, intent(in) :: Liouv(dim,dim)
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), intent(in), optional :: name

        integer dimsq, i, j, k,l
        doublecomplex, allocatable :: Lliouv(:,:), rhoLiouv(:), RVec(:,:), LVec(:,:), Values(:),&
                rho_t(:), proj(:)
        dimsq = dim**2
        allocate(Lliouv(dimsq, dimsq), rhoLiouv(dimsq), rho_t(dimsq))
        rhoLiouv = ReturnOperatorsInLiouvillian(rho, dim)
        Lliouv = ReturnSuperOperatorsInLiouvillian(Liouv, dim)
        !!!!!!!!!!!!!!!
        !Diagonalizing the liouvillian
        !!!!!!!!!!!!!!!
        call LiovullianDiagonalization(LVec, RVec, Values, Lliouv, dimsq)
        !!!!!!!!!!!!!!!
        ! Projecting the rho into the left eigenstate
        !!!!!!!!!!!!!!!
        allocate(proj(dimsq))
        do i=1, dimsq
            proj(i) = (0.d0, 0.d0)
            do k=1, dimsq
                proj(i) = proj(i)+conjg(Lvec(k,i))*rhoLiouv(k)
            end do
        end do
        allocate(rho_t(dimsq))
        write(*,*) 'Start evolution by Diagonalizing the Liouvillian'
        if(present(name))then
            if(name.ne."*")then
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    do j=1, dimsq
                        rho_t(j) = (0.d0, 0.d0)
                        do k=1, dimsq
                            rho_t(j) = rho_t(j) + RVec(j,k)*proj(k)*e_nepero**(Values(k)*i*Deltat)
                        end do
                    end do
                    write(1) rho(1:,1:)
                end do
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F)") 'Time: ', Deltat*i
                    do j=1, dimsq
                        rho_t(j) = (0.d0, 0.d0)
                        do k=1, dimsq
                            rho_t(j) = rho_t(j) + RVec(j,k)*proj(k)*e_nepero**(Values(k)*i*Deltat)
                        end do

                    end do
                    do j=1, dim
                        write(*,'(<dim>(2(e9.2, x), x))') (rho_t((j-1)*dim+k), k=1, dim)
                    end do
                end do
            end if
        else
            do j=1, dimsq
                rho_t(j) = (0.d0, 0.d0)
                do k=1, dimsq
                    rho_t(j) = rho_t(j) + RVec(j,k)*proj(k)*e_nepero**(Values(k)*timesteps*Deltat)
                end do
            end do
        end if
        rho = ReturnOperatorsInRealSpace(rho_t, dim)
        write(*,*) 'Evolved with the Diagonalized Liouvillian!'
        deallocate(rhoLiouv, lLiouv, proj, RVec, LVec, Values, rho_t)
    end subroutine EvolveWithLiouvilleDiagonalization

    subroutine EvolveWithArnoldi(rho, red, dim, nkrylov,Deltat, timesteps, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> optional doublecomplex (dim,dim,dim,dim) redfield tensor
        doublecomplex, optional, intent(in) :: red(dim,dim,dim,dim)
        !> optional integer dimension of the krylov space
        integer, optional, intent(in) :: nkrylov
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name
        integer dimsq
        doublecomplex, allocatable :: Liouvillian(:,:), rhovec(:)
        dimsq = dim**2
        allocate(Liouvillian(dimsq, dimsq), rhovec(dimsq))
        Liouvillian = ReturnSuperOperatorsInLiouvillian(red, dim)
        rhovec = ReturnOperatorsInLiouvillian(rho, dim)
        write(*,*) 'Start evolution with Arnoldi'
        if(present(name))then
            if(name.ne."*")then
                write(*,*) 'I write to ', name
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call arnoldi(rhovec, Liouvillian, n_krylov,dim, Deltat)
                    rho = ReturnOperatorsInRealSpace(rhovec, dim)
                end do
                write(1) rho(1:,1:)
                close(1)
            else
                write(*,*) 'Print in output'
                do i=0, timesteps-1
                    write(*,"(A, F)") 'Time: ', Deltat*i
                    call write_matrix_complex(rho, dim, dim, 'e11.3')
                    call arnoldi(rhovec, Liouvillian, n_krylov,dim, Deltat)
                    rho = ReturnOperatorsInRealSpace(rhovec, dim)
                end do
                write(*,"(A, F)") 'Time: ', Deltat*(i+1)
                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call arnoldi(rhovec, Liouvillian, n_krylov,dim, Deltat)
            end do
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
        rho = ReturnOperatorsInRealSpace(rhovec, dim)
        deallocate(Liouvillian, rhovec)
    end subroutine EvolveWithArnoldi

    subroutine CorrelationWithRungeKuttaUnitary(op1, rho, H, dim, Deltat, timesteps, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) operator for the correlation function
        doublecomplex, intent(inout) :: op1(dim,dim)
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim) Hamiltonian
        doublecomplex, intent(in) :: H(dim,dim)
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name

        integer i
        doublecomplex correlation
        doublecomplex, allocatable :: correlation_mat(:,:)
        allocate(correlation_mat(dim,dim))

        write(*,*) 'Start evolution with Runge Kutta'
        if(present(name))then
            if(name.ne."*")then
                open(1, file=name)
                do i=0, timesteps-1
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(1, '(F, x, 2(e15.3E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call Runge_Kutta_density(rho, H, dim, Deltat)
                end do
                call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                correlation = trace_complex(correlation_mat, dim)
                write(1, '(F, x, 2(e15.3E5, x))') deltat*i,real(correlation), aimag(correlation)
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F,A,F)") 'Time: ', Deltat*i, 'trace=', real(trace_complex(rho, dim))
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(*, '(F, x, 2(e15.3E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call Runge_Kutta_density(rho, H, dim, Deltat)
                end do
                write(*,"(A, F)") 'Time: ', Deltat*(i+1)
                call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                correlation = trace_complex(correlation_mat, dim)
                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            write(*,*) 'ERROR! I need a valid file in which to write'
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
        deallocate(correlation_mat)
    end subroutine CorrelationWithRungeKuttaUnitary

end module Evolutions