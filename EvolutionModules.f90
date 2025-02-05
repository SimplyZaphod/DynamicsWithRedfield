module Evolutions
    use conversion
    use Mathfunc
    use Redfield
    use UsefullIO
    contains



    !This routine fill a Gamma matrix given as input with the desired spectral density
    !> @todo check if I imposed the optional constant in a correct way
    subroutine FilLGamma(dim, Gammas, H, GammaType, const, eta)
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
        !> doubleprecision optional constant in case of debye density
        doubleprecision,optional, intent(in) :: eta

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
        elseif((GammaType=='Debye').or.(GammaType=='Db'))then
            write(*,*) 'Debye spectral density'
            call GammaDebye(gammatoken, values, dim, const, eta)
        elseif((GammaType=='null').or.(GammaType=='zero'))then
            write(*,*) 'Null spectral density'
            call GammaConstant(gammatoken, 0.d0, values, dim)
        else
            write(*,*) 'Unrecognized keyword! Null is guessed'
            call GammaConstant(gammatoken, 0.d0, values, dim)
        endif
        Gammas = gammatoken*reality
    end subroutine FilLGamma

    !This routine fill a density matrix given as input with the desired population distribution
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

    subroutine EvolveWithRungeKuttaUnitary(rho, H, dim, Deltat, timesteps,n_prop, prop, name)
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
        !>integer additional numbero f properties i want to calculate
        integer, intent(in) :: n_prop
        !> doublecomplex (n_prop, dim, dim) matrix of properties
        doublecomplex, intent(in) :: prop(n_prop, dim, dim)
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name

        doublecomplex, allocatable :: propts(:,:)
        character(len=50) :: spettro
        doublecomplex property

        do i=1, n_prop
            spettro= ''
            write(spettro,*) 'prop', i,'.dat'
            call StripSpaces(spettro)
            open(i, file=spettro)
        end do
        write(*,*) 'Start evolution with Runge Kutta'
        if(present(name))then
            write(*,*) 'Save at: ', name
            if(name.ne."*")then
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call Runge_Kutta_Unitary(rho, H, dim, Deltat,1)
                end do
                write(1) rho(1:,1:)
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F12.4,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    write(*,*) ' Population:'
                    write(*,'(<dim>(F7.3, X))') (real(rho(j,j)), j=1, dim)
                    if(n_prop.gt.0)then
                        allocate(propts(dim, dim))
                        do j=1, n_prop
                            propts = prop(j,:,:)
                            call npmatmul_complex(propts, propts, rho, dim, dim, dim)
                            property = trace_complex(propts, dim)
                            write(*,'(A, I2)') 'Property number ', j
                            write(*,'(A, e18.6E5, 2X, e18.6E5)') "Value ", real(property), aimag(property)
                            write(j,'(F12.4, x, e18.6E5, x, e18.6E5)') Deltat*i, real(property), aimag(property)
                        end do
                        deallocate(propts)
                    end if
                    call Runge_Kutta_Unitary(rho, H, dim, Deltat,1)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call Runge_Kutta_Unitary(rho, H, dim, Deltat,1)
            end do
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
        do  i=1, n_prop
            close(i)
        end do
    end subroutine EvolveWithRungeKuttaUnitary

    subroutine EvolveWithUnitary(rho, H, dim, Deltat, timesteps,n_prop, prop, name)
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
        !>integer additional numbero f properties i want to calculate
        integer, intent(in) :: n_prop
        !> doublecomplex (n_prop, dim, dim) matrix of properties
        doublecomplex, intent(in) :: prop(n_prop, dim, dim)
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name

        doublecomplex, allocatable :: propts(:,:)
        character(len=50) :: spettro
        doublecomplex property

        do i=1, n_prop
            spettro= ''
            write(spettro,*) 'prop', i,'.dat'
            call StripSpaces(spettro)
            open(i, file=spettro)
        end do
        write(*,*) 'Start unitary evolution'
        if(present(name))then
            write(*,*) 'Save at: ', name
            if(name.ne."*")then
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call Unitary_Evolution(rho, H, dim, Deltat,1)
                end do
                write(1) rho(1:,1:)
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F12.4,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    write(*,*) ' Population:'
                    write(*,'(<dim>(F7.3, X))') (real(rho(j,j)), j=1, dim)
                    if(n_prop.gt.0)then
                        allocate(propts(dim, dim))
                        do j=1, n_prop
                            propts = prop(j,:,:)
                            call npmatmul_complex(propts, propts, rho, dim, dim, dim)
                            property = trace_complex(propts, dim)
                            write(*,'(A, I2)') 'Property number ', j
                            write(*,'(A, e18.6E5, 2X, e18.6E5)') "Value ", real(property), aimag(property)
                            write(j,'(F12.4, x, e18.6E5, x, e18.6E5)') Deltat*i, real(property), aimag(property)
                        end do
                        deallocate(propts)
                    end if
                    call Unitary_Evolution(rho, H, dim, Deltat,1)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call Unitary_Evolution(rho, H, dim, Deltat,1)
            end do
        end if
        write(*,*) 'Evolved with Unitary Evolution!'
        do  i=1, n_prop
            close(i)
        end do
    end subroutine EvolveWithUnitary



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
            write(*,*) 'Save at: ', name
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
                    write(*,"(A, F12.4,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    write(*,*) ' Population:'
                    write(*,'(<dim>(F7.3, X))') (real(rho(j,j)), j=1, dim)
!                    call write_matrix_complex(rho, dim, dim, 'e11.3')
                    call Runge_Kutta_Liouville(rho, liouv, dim, Deltat)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
!                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call Runge_Kutta_Liouville(rho, liouv, dim, Deltat)
            end do
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
    end subroutine EvolveWithRungeKutta

        subroutine EvolveWithRungeKuttaRedfield(rho, H, red, dim, Deltat, timesteps,n_prop, prop, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim) Hamiltonian
        doublecomplex, intent(in) :: H(dim, dim)
        !> doublecomplex(dim,dim, dim, dim) redfield tensor
        doublecomplex, intent(in) :: red(dim,dim, dim, dim)
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !>integer additional numbero f properties i want to calculate
        integer, intent(in) :: n_prop
        !> doublecomplex (n_prop, dim, dim) matrix of properties
        doublecomplex, intent(in) :: prop(n_prop, dim, dim)
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name

        doublecomplex, allocatable :: propts(:,:)
        character(len=50) :: spettro
        doublecomplex property

        do i=1, n_prop
            spettro= ''
            write(spettro,*) 'prop', i,'.dat'
            call StripSpaces(spettro)
            open(i, file=spettro)
        end do
        write(*,*) 'Start evolution with Runge Kutta'
        if(present(name))then
            write(*,*) 'Save at: ', name
            if(name.ne."*")then
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call Runge_Kutta_Redfield(rho,H, red, dim, Deltat,1)
                end do
                write(1) rho(1:,1:)
                close(1)
            else
                do i=0, timesteps-1
                    write(*,"(A, F12.4,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    write(*,*) ' Population:'
                    write(*,'(<dim>(F7.3, X))') (real(rho(j,j)), j=1, dim)
!                    call write_matrix_complex(rho, dim, dim, 'e11.3')
                    if(n_prop.gt.0)then
                        allocate(propts(dim, dim))
                        do j=1, n_prop
                            propts = prop(j,:,:)
                            call npmatmul_complex(propts, propts, rho, dim, dim, dim)
                            property = trace_complex(propts, dim)
                            write(*,'(A, I2)') 'Property number ', j
                            write(*,'(A, e18.6E5, 2X, e18.6E5)') "Value ", real(property), aimag(property)
                            write(j,'(F12.4, x, e18.6E5, x, e18.6E5)') Deltat*i, real(property), aimag(property)
                        end do
                        deallocate(propts)
                    end if
                    call Runge_Kutta_Redfield(rho,H, red, dim, Deltat,1)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
!                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call Runge_Kutta_Redfield(rho,H, red, dim, Deltat,1)
            end do
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
        do  i=1, n_prop
            close(i)
        end do
    end subroutine EvolveWithRungeKuttaRedfield

    subroutine EvolveWithLiouvilleDiagonalization(rho, Liouv, dim, Deltat, timesteps, n_prop, prop,name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim, dim, dim) Liouvillian
        doublecomplex, intent(in) :: Liouv(dim,dim, dim, dim)
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !>integer additional numbero f properties i want to calculate
        integer, intent(in) :: n_prop
        !> doublecomplex (n_prop, dim, dim) matrix of properties
        doublecomplex, intent(in) :: prop(n_prop, dim, dim)
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), intent(in), optional :: name

        integer dimsq, i, j, k,l
        doublecomplex tracee
        doublecomplex, allocatable :: Lliouv(:,:), rhoLiouv(:), RVec(:,:), LVec(:,:), Values(:),&
                rho_t(:), proj(:)
        doublecomplex, allocatable :: propts(:,:)
        character(len=50) :: spettro
        doublecomplex property


        dimsq = dim**2
        allocate(Lliouv(dimsq, dimsq), rhoLiouv(dimsq), rho_t(dimsq))
        rhoLiouv = ReturnOperatorsInLiouvillian(rho, dim)
        Lliouv = ReturnSuperOperatorsInLiouvillian(Liouv, dim)
        !!!!!!!!!!!!!!!
        !Diagonalizing the liouvillian to obtain left eigenvectors {l_i}, right eigenvectors {r_i} and eigenvalues {\lambda_i}
        ! so that I can make use of |rho(t)> = \sum_i |l_i><r_i|rho(0)> exp^(t*\lambda_i)
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

        do i=1, n_prop
            spettro= ''
            write(spettro,*) 'prop', i,'.dat'
            call StripSpaces(spettro)
            open(i, file=spettro)
        end do

        write(*,*) 'Start evolution by Diagonalizing the Liouvillian'
        if(present(name))then
            write(*,*) 'Save at: ', name
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
                    tracee = (0.d0, 0.d0)
                    do j=1, dim
                        tracee =tracee+ rho_t((j-1)*dim + j)
                    end do
                    write(*,"(A, F12.4, A, F12.4)") 'Time: ', Deltat*i, ',trace= ', real(tracee)
                    write(*,*) ' Population:'
                    write(*,'(<dim>(F7.3, X))') (real(rho_t((j-1)*dim + j)), j=1, dim)
                    do j=1, dimsq
                        rho_t(j) = (0.d0, 0.d0)
                        do k=1, dimsq
                            rho_t(j) = rho_t(j) + RVec(j,k)*proj(k)*e_nepero**(Values(k)*i*Deltat)
                        end do
                    end do
!                    do j=1, dim
!                        write(*,'(<dim>(2(e11.4, x), x))') (rho_t((j-1)*dim+k), k=1, dim)
!                    end do
                    rho = ReturnOperatorsInRealSpace(rho_t, dim)
                    if(n_prop.gt.0)then
                        allocate(propts(dim, dim))
                        do j=1, n_prop
                            propts = prop(j,:,:)
                            call npmatmul_complex(propts, propts, rho, dim, dim, dim)
                            property = trace_complex(propts, dim)
                            write(*,'(A, I2)') 'Property number ', j
                            write(*,'(A, e18.6E5, 2X, e18.6E5)') "Value ", real(property), aimag(property)
                            write(j,'(F12.4, x, e18.6E5, x, e18.6E5)') Deltat*i, real(property), aimag(property)
                        end do
                        deallocate(propts)
                    end if
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

    subroutine EvolveWithArnoldi(rho, red, dim, nkrylov,Deltat, timesteps, n_prop, prop, name)
        implicit none
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex (dim,dim,dim,dim) redfield tensor
        doublecomplex, intent(in) :: red(dim,dim,dim,dim)
        !> optional integer dimension of the krylov space
        integer, intent(in) :: nkrylov
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !>integer additional numbero f properties i want to calculate
        integer, intent(in) :: n_prop
        !> doublecomplex (n_prop, dim, dim) matrix of properties
        doublecomplex, intent(in) :: prop(n_prop, dim, dim)
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name

        integer dimsq, i, j, k, l, kl
        doublecomplex, allocatable :: Liouvillian(:,:), rhovec(:)
        doublecomplex property
        character(len=50) spettro

        dimsq = dim**2
        allocate(Liouvillian(dimsq, dimsq), rhovec(dimsq))
        Liouvillian = ReturnSuperOperatorsInLiouvillian(red, dim)
        rhovec = ReturnOperatorsInLiouvillian(rho, dim)

        do i=1, n_prop
            spettro= ''
            write(spettro,*) 'prop', i,'.dat'
            call StripSpaces(spettro)
            open(i, file=spettro)
        end do

        write(*,*) 'Start evolution with Arnoldi'
        if(present(name))then
            write(*,*) 'Save at: ', name
            if(name.ne."*")then
                write(*,*) 'I write to ', name
                open(1, file=name, form='unformatted')
                do i=0, timesteps-1
                    write(1) rho(1:,1:)
                    call arnoldi(rhovec, Liouvillian,dimsq, nkrylov, Deltat)
                    rho = ReturnOperatorsInRealSpace(rhovec, dim)
                end do
                write(1) rho
                close(1)
            else
                write(*,*) 'Print in output'
                do i=0, timesteps-1
                    property = (0.d0, 0.d0)

!                   ! $omp parallel do default(none), &
!                  !  $omp shared(dim, rhovec),&
!                   ! $omp private(j),&
!                   ! $omp reduction(+:property)
                    do j=1, dim
                        property = property + rhovec((j-1)*dim+j)
                    end do
                    !!$omp end parallel do
                    write(*,"(A, F12.4,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(property)
                    write(*,*) ' Population:'
                    write(*,'(<dim>(F7.3, X))') (real(rhovec((j-1)*dim+j)), j=1, dim)

                    if(n_prop.gt.0)then
                        do j=1, n_prop
                            property = (0.d0, 0.d0)
                            do k=1, dim
                                do l=1, dim
                                    kl = (k-1)*dim+l
                                    property = property + prop(j,l,k)*rhovec(kl)
                                end do
                            end do
                            write(*,'(A, I2)') 'Property number ', j
                            write(*,'(A, e18.6E5, 2X, e18.6E5)') "Value ", real(property), aimag(property)
                            write(j,'(F12.4, x, e18.6E5, x, e18.6E5)') Deltat*i, real(property), aimag(property)
                        end do
                    end if
                    call arnoldi(rhovec, Liouvillian, dimsq,nkrylov, Deltat)
                end do
!                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
!                call write_matrix_complex(rho, dim, dim, 'e11.3')
                write(*,*)
            end if
        else
            do i=0, timesteps-1
                call arnoldi(rhovec, Liouvillian, dimsq, nkrylov, Deltat)
            end do
        end if
        write(*,*) 'Evolved with Runge-Kutta!'
        rho = ReturnOperatorsInRealSpace(rhovec, dim)
        deallocate(Liouvillian, rhovec)
    end subroutine EvolveWithArnoldi

    subroutine CorrelationWithArnoldi(op1, rho, Liouv, dim,nkrylov, Deltat, timesteps, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) operator for the correlation function
        doublecomplex, intent(inout) :: op1(dim,dim)
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim,dim,dim) Liouvillian
        doublecomplex, intent(in) :: Liouv(dim,dim,dim,dim)
        !> dimension of the krylov space
        integer, intent(in) :: nkrylov
        !> doubleprecision timesteps in femtoseconds
        doubleprecision, intent(in) :: Deltat
        !> integer number of timesteps
        integer, intent(in) :: timesteps
        !> do you want to write? if "*" write on output, if name, write on filename
        character(len=*), optional :: name

        integer i, dimsq
        doublecomplex correlation
        doublecomplex, allocatable :: correlation_mat(:,:), lliouv(:,:), lrho(:)
        allocate(correlation_mat(dim,dim))

        dimsq = dim**2
        allocate(lliouv(dimsq, dimsq), lrho(dimsq))

        lliouv = ReturnSuperOperatorsInLiouvillian(Liouv, dim)
        lrho = ReturnOperatorsInLiouvillian(rho, dim)
        write(*,*) 'Start evolution with Arnoldi'
        if(present(name))then
            write(*,*) 'Save at: ', name
            if(name.ne."*")then
                open(1, file=name)
                do i=0, timesteps-1, 1
                    rho = ReturnOperatorsInRealSpace(lrho, dim)
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(1, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call arnoldi(lrho, lliouv, dimsq, nkrylov, Deltat)
                end do
                call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                correlation = trace_complex(correlation_mat, dim)
                write(1, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                close(1)
            else
                do i=0, timesteps-1, 1
                    write(*,"(A, F12.4 ,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(*, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call arnoldi(lrho, lliouv, dimsq, nkrylov, Deltat)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
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
    end subroutine CorrelationWithArnoldi

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
                do i=0, timesteps-1, 1
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(1, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call Runge_Kutta_Unitary(rho, H, dim, Deltat, 1)
                end do
                call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                correlation = trace_complex(correlation_mat, dim)
                write(1, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                close(1)
            else
                do i=0, timesteps-1, 1
                    write(*,"(A, F12.4 ,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(*, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call Runge_Kutta_Unitary(rho, H, dim, Deltat, 1)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
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

    subroutine CorrelationWithRungeKuttaRedfield(op1, rho, H, red, dim, Deltat, timesteps, name)
        !> integer dimension of all the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) operator for the correlation function
        doublecomplex, intent(inout) :: op1(dim,dim)
        !> doublecomplex(dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doublecomplex(dim,dim) Hamiltonian
        doublecomplex, intent(in) :: H(dim,dim)
        !> doublecomplex(dim,dim, dim, dim) Redfield
        doublecomplex, intent(in) :: red(dim,dim, dim, dim)
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
                do i=0, timesteps-1, 1
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(1, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call Runge_Kutta_Redfield(rho, H, red, dim, Deltat, 1)
                end do
                call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                correlation = trace_complex(correlation_mat, dim)
                write(1, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                close(1)
            else
                do i=0, timesteps-1, 1
                    write(*,"(A, F12.4 ,A,F12.4)") 'Time: ', Deltat*i, ',trace= ', real(trace_complex(rho, dim))
                    call npmatmul_complex(correlation_mat, op1, rho, dim, dim, dim)
                    correlation = trace_complex(correlation_mat, dim)
                    write(*, '(F, x, 2(e18.6E5, x))') deltat*i,real(correlation), aimag(correlation)
                    call Runge_Kutta_Redfield(rho, H, red, dim, Deltat, 1)
                end do
                write(*,"(A, F12.4)") 'Time: ', Deltat*(i+1)
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
    end subroutine CorrelationWithRungeKuttaRedfield

end module Evolutions