module Redfield
    !> @brief This module contains the routine to execute the time evolution of a matrix/operator in contact with a Bath. it also contains some usefull function like
    !> the writing of Redfield tensor, some spectral density, the reading of some input, some integration.
    !> It needs the module conversion.o, InputOutput.o and MathFunc.o
    use conversion
    use usefullIO
    use mathfunc
    implicit none
contains
    !> @brief This subroutine simply read the info to execute the dynamics from input file
    subroutine load_time_evolution_input(deltat, timesteps, n_dt, n_krylov, nameinput)
        !> @param doublepreciswion timestep
        doubleprecision, intent(out) :: deltat
        !> @param integer Number of desidered timesteps for the dynamic
        integer, intent(out) :: timesteps
        !> @param integer number of steps before save
        integer, intent(out) :: n_dt
        !> @param integer dimension of the krylov space for arnoldi
        integer, intent(out) :: n_krylov
        !> @param nameinput Optional name of the input. If not inserted, "generaltime.inp" will be selected
        character(len=*), intent(in), optional :: nameinput

        character(len=150) name
        if(present(nameinput))then
            name = nameinput
        else
            write(name,*) 'generaltime.inp'
        end if
        open(1, file=name)
        read(1, *) deltat, timesteps, n_dt
        read(1,*) n_krylov
        write(*,*) 'Debu2'
        write(*,'(A, A30)') 'Time evolution parameters read from', name
        write(*,*) 'Deltat:', deltat,'fs, n° of timesteps:', timesteps
        write(*,'(A,F,A)') 'from t=0 to t=', deltat*(timesteps-1),' fs'
        write(*,*) 'I will save every ', n_dt,' snapshots'
        write(*,*) 'Dimension of the Krylov space:', n_krylov
        write(*,*) '========================================='
        close(1)
    end subroutine load_time_evolution_input

    subroutine ReadTimeEvolutionInput(nameinput,deltat, timesteps, n_dt, DynType,T, densitytype, gammatype, omega_g)
        !> From where I need to read the information
        character(len=*), intent(in) :: nameinput
        !> doubleprecision Timestep
        doubleprecision, intent(out) :: deltat
        !> integer number of timesteps
        integer, intent(out) :: timesteps
        !> Integer number of timesteps before printing
        integer, intent(out) :: n_dt
        !> character type of the dynamics. Allowed value: 'arnoldi', 'runge-kutta', 'diag', 'unitary'
        character(len=100), intent(out) :: DynType
        !> doubleprecision temperature
        doubleprecision, intent(out) :: T
        !> character density needed. Allowed value: 'BE' or 'BoseEinstein', 'dark' or 'darkness'
        character(len=100), intent(out) :: DensityType
        !> character type of spectral density. Allowed value: 'const' or 'constant' or 'G', 'emissive' or 'radiation'
        character(len=100), intent(out) :: GammaType
        !> doubleprecision constant in case of costant or Debye spectral density. Write whatever you want if emissive or radiation is requested
        doubleprecision, intent(out) :: omega_g

        open(1, file=nameinput)
        read(1, *) deltat
        read(1,*) timesteps
        read(1,*) n_dt
        read(1,*) DynType
        read(1,*) T
        read(1,*) DensityType
        read(1,*) GammaType
        read(1,*) omega_g

        write(*,'(A, A30)') 'Time evolution parameters read from: ', nameinput
        write(*,*) 'Deltat:', deltat,'fs, n° of timesteps:', timesteps
        write(*,'(A,F,A)') 'from t=0 to t=', deltat*(timesteps-1),' fs'
        write(*,*) 'I will save every ', n_dt,' snapshots'
        write(*,*) 'Integration method:', DynType
        write(*,'(A, F7.2,A)') 'Temperature: ', T, ' K'
        write(*,*) 'Type of density:', DensityType
        write(*,*) 'GammaType: ', GammaType
        write(*,*) 'Eventual constant: ', omega_g
        write(*,*) '========================================='
        close(1)
    end subroutine ReadTimeEvolutionInput


    !> Calculation of the weight for the Simpson integration
    !> @todo I need to put this routine in the MathFunc.o module
    double precision function weight(i,start,end)
        implicit none
        !> Integer step that I am integrating
        integer, intent(in) :: i
        !> Integer first number of the vectors
        integer, intent(in) :: start
        !> Integer last number of the vector
        integer, intent(in) :: end

        if (i.eq.start.or.i.eq.end)then
            weight=1.d0/3.d0
        else
            if (mod(i,2).eq.0) then
                weight=4.d0/3.d0
            else
                weight=2.d0/3.d0
            endif
        endif
        return
    end function weight
    !prova
    !> Function to create the Redfield tensor in dobleprecision
    !> Write the Redfield Tensor as a double N*N*N*N matrix having given the generic operators to write all the elements
    !> @note val is in femtosecond!
    function RedfieldTensor_Double(adag, ahat, gammas, density, val, dim)
        !> integer dimension of the Redfield tensor and all the operator
        integer, intent(in) :: dim
        !> Doubleprecision (dim,dim) matrix, it is the firs operator that need to be coupled to the bath
        doubleprecision, intent(in) :: adag(dim, dim)
        !> Doubleprecision (dim,dim) matrix, it is the second operator that need to be coupled to the bath
        doubleprecision, intent(in) :: ahat(dim, dim)
        !> Spectral density as a doubleprecision (dim,dim) matrix
        doubleprecision, intent(in) :: gammas(dim, dim)
        !> Density of the field as a doubleprecision (dim, dim) matrix
        doubleprecision, intent(in) :: density(dim, dim)
        !> Eigenvalues doubleprecision (dim) matrix
        doubleprecision, intent(in) :: val(dim)
        !> Redfield Tensor doubleprecision (dim, dim, dim,dim) matrix as output
        doubleprecision :: RedfieldTensor_Double(dim,dim,dim,dim)

        integer a,b,c,d,m

        RedfieldTensor_Double = 0.d0
        do a=1, dim
            do b=1, dim
                do c=1, dim
                    do d=1, dim
                        if(b==d) then
                            do m=1, dim
                                RedfieldTensor_Double(a,b,c,d) = RedfieldTensor_Double(a,b,c,d) -&
                                        GammaPlus_double(adag, ahat, gammas, val,density, a,m,m,c, dim)
                            end do
                        end if
                        if(a==c)then
                            do m=1, dim
                                RedfieldTensor_Double(a,b,c,d) = RedfieldTensor_Double(a,b,c,d) -&
                                        GammaMinus_double(adag, ahat, gammas, val,density, d,m,m,b, dim)
                            end do
                        end if
                        RedfieldTensor_Double(a,b,c,d) = RedfieldTensor_Double(a,b,c,d) +&
                                GammaPlus_double(adag, ahat, gammas, val,density, d, b,a,c, dim)+&
                                GammaMinus_double(adag, ahat, gammas, val, density,d,b,a,c, dim)

                    end do
                end do
            end do
        end do
    end function RedfieldTensor_Double

    !> write the GammaPlus as double matrix element for the writing of Redfield Tensor
    !> @note take care of the order of the element!
    !> @note val in femtosecond
    function GammaPlus_double(op1, op2, gammas, val, density, a,b,c,d, dim)
        integer, intent(in) :: dim,a,b,c,d
        doubleprecision, intent(in) :: op1(dim,dim), op2(dim, dim), gammas(dim, dim), val(dim), density(dim,dim)
        doubleprecision :: GammaPlus_double
        doubleprecision omega
        GammaPlus_double = 0.d0
        omega = val(c)-val(d)
        if(omega.lt.0.d0)then
            GammaPlus_double = op1(a,b)*op2(c,d)*gammas(d,c)*(density(d,c)+1)
        elseif(omega.gt.0.d0) then
            GammaPlus_double = op2(a,b)*op1(c,d)*gammas(c,d)*(density(c,d))
        end if
    end function GammaPlus_double

    !> write the GammaMinus as double matrix element for the writing of Redfield Tensor
    !> @note take care of the order of the element!
    !> @note val in femtosecond
    function GammaMinus_double(op1, op2, gammas, val, density, a,b,c,d, dim)
        integer, intent(in) :: dim, a,b,c,d
        doubleprecision, intent(in) :: op1(dim,dim), op2(dim, dim), gammas(dim, dim), val(dim),density(dim,dim)
        doubleprecision :: GammaMinus_double
        doubleprecision omega

        GammaMinus_double = 0.d0
         omega = val(a)-val(b)
        if(omega.gt.0.d0)then
            GammaMinus_double = op1(a,b)*op2(c,d)*gammas(a,b)*(density(a,b)+1)
        else if(omega.lt.0.d0)then
        GammaMinus_double = GammaMinus_double+ op2(a,b)*op1(c,d)*gammas(b,a)*(density(b,a))
        end if
    end function GammaMinus_double

    !> This subroutine check the detail balance condition R_iijj/R_jjii = exp(-h omega_ij/kT)
    !> and then print the output to compare
    !> @note eigenvalues in eV
    subroutine CheckDetailBalanceRedfield_Double(Redfield, val, T, dim)
        !> dimension of all the vectors
        integer, intent(in) :: dim
        !> Redfield tensor as a (dim,dim,dim,dim) double matrix
        doubleprecision, intent(in) :: Redfield(dim,dim,dim,dim)
        !> Eigenvalues doubleprecision (dim) array
        doubleprecision, intent(in) :: val(dim)
        !> Temperature doubleprecision DEPRECATED
        doubleprecision, intent(in) :: T

        integer i,j,k,l
        doubleprecision wastedouble1, wastedouble2
        write(*,*) 'Check of the Redfield tensor using the Detail balance'
        write(*,*) 'If the distribution is the Bose-Einstein, I must have that'
        write(*,*) 'R_iijj/R_jjii *1/exp(-homega/kbT) = 1'
        write(*,*) '  i  |  j  |  Ratio  '
        do i=1, dim
            do j=i, dim
                wastedouble1 = e_nepero**(-(val(i)-val(j))/(KB_eV*T))
                wastedouble2 = 0.d0
                if (Redfield(j,j,i,i).ne.0.d0) wastedouble2 = Redfield(i,i,j,j)/Redfield(j,j,i,i)
                if (Redfield(i,i,j,j).eq.0.d0.and.Redfield(j,j,i,i).eq.0.d0) wastedouble2 = wastedouble1
                write(*,'(2X,i3,A3,i3,A3,F7.4)') i,' | ',j,' | ', wastedouble2/wastedouble1
            end do
        end do
    end subroutine CheckDetailBalanceRedfield_Double

    !> Function to create the Redfield tensor in doblecomplex
    !> Write the Redfield Tensor as a double N*N*N*N matrix having given the generic operators to write all the elements
    !> @note val is in femtosecond!
    function RedfieldTensor_Complex(adag, ahat, gammas, density, val, dim)
        !> integer dimension of the Redfield tensor and all the operator
        integer, intent(in) :: dim
        !> doublecomplex (dim,dim) matrix, it is the firs operator that need to be coupled to the bath
        doublecomplex, intent(in) :: adag(dim, dim)
        !> doublecomplex (dim,dim) matrix, it is the second operator that need to be coupled to the bath
        doublecomplex, intent(in) :: ahat(dim, dim)
        !> Spectral density as a doublecomplex (dim,dim) matrix
        doublecomplex, intent(in) :: gammas(dim, dim)
        !> Density of the field as a doubleprecision (dim, dim) matrix
        doubleprecision, intent(in) :: density(dim, dim)
        !> Eigenvalues doubleprecision (dim) matrix
        doubleprecision, intent(in) :: val(dim)
        !> Redfield Tensor doublecomplex (dim, dim, dim,dim) matrix as output
        doublecomplex, allocatable :: RedfieldTensor_complex(:,:,:,:)
        integer a,b,c,d,m

        allocate(RedfieldTensor_Complex(dim, dim, dim, dim))
        RedfieldTensor_Complex = (0.d0, 0.d0)
        do a=1, dim
            do b=1, dim
                do c=1, dim
                    do d=1, dim
                        if(b==d) then
                            do m=1, dim
                                RedfieldTensor_Complex(a,b,c,d) = RedfieldTensor_Complex(a,b,c,d) -&
                                        GammaPlus_complex(adag, ahat, gammas, val,density, a,m,m,c, dim)
                            end do
                        end if
                        if(a==c)then
                            do m=1, dim
                                RedfieldTensor_Complex(a,b,c,d) = RedfieldTensor_Complex(a,b,c,d) -&
                                        GammaMinus_complex(adag, ahat, gammas, val,density, d,m,m,b, dim)
                            end do
                        end if
                        RedfieldTensor_Complex(a,b,c,d) = RedfieldTensor_Complex(a,b,c,d) +&
                                GammaPlus_complex(adag, ahat, gammas, val,density, d, b,a,c, dim)+&
                                GammaMinus_complex(adag, ahat, gammas, val, density,d,b,a,c, dim)

                    end do
                end do
            end do
        end do
    end function RedfieldTensor_Complex

    !> write the GammaPlus as double complex element for the writing of Redfield Tensor
    !> @note take care of the order of the element!
    !> @note val in femtosecond
    function GammaPlus_Complex(op1, op2, gammas, val, density, a,b,c,d, dim)
        integer dim, a,b,c,d
        doubleprecision  density(dim, dim), omega, val(dim)
        doublecomplex op1(dim,dim), op2(dim, dim), gammas(dim, dim),GammaPlus_Complex

        GammaPlus_Complex = (0.d0, 0.d0)
        omega = val(c)-val(d)
        if(omega.lt.0.d0)then
            GammaPlus_Complex = op1(a,b)*op2(c,d)*gammas(d,c)*(density(d,c)+1)
        elseif(omega.gt.0.d0) then
            GammaPlus_Complex = op2(a,b)*op1(c,d)*gammas(c,d)*(density(c,d))
        end if
    end function GammaPlus_Complex

    !> write the GammaMinux as doublecomplex matrix element for the writing of Redfield Tensor
    !> @note take care of the order of the element!
    !> @note val in femtosecond
    function GammaMinus_Complex(op1, op2, gammas, val, density, a,b,c,d, dim)
        integer dim, a,b,c,d
        doubleprecision density(dim, dim), omega, val(dim)
        doublecomplex op1(dim,dim), op2(dim, dim), gammas(dim, dim),GammaMinus_Complex

        GammaMinus_Complex = (0.d0, 0.d0)
        omega = val(a)-val(b)
        if(omega.gt.0.d0)then
            GammaMinus_Complex = op1(a,b)*op2(c,d)*gammas(a,b)*(density(a,b)+1)
        else if(omega.lt.0.d0)then
            GammaMinus_Complex = GammaMinus_Complex+ op2(a,b)*op1(c,d)*gammas(b,a)*(density(b,a))
        end if
    end function GammaMinus_Complex

    !> This algotirhm execute the rotate wave approximation on a redfield vector
    subroutine RotatingWaveApproximation_complex(redfield, dim)
        !> integer dimension of the array
        integer, intent(in) :: dim
        !> doublecomplex (dim,dim,dim,dim) redfield tensor
        doublecomplex, intent(inout) :: redfield(dim, dim, dim, dim)
        integer i,j,k,l
        do i=1, dim
            do j=1, dim
                if(i.ne.j) then
                    do k=1, dim
                        do l=1, dim
                            if (k.ne.i.and.l.ne.j) redfield(i,j,k,l) = (0.d0, 0.d0)
                        end do
                    end do
                end if
            end do
        end do
    end subroutine RotatingWaveApproximation_complex

    !> Create a Gamma spectral density with real constant amplitude and null imaginary amplitude
    subroutine GammaConstant(gammas,gamma, eigenvalues, dim)
        !> integer dimension of all the arrays
        integer, intent(in) :: dim
        !> doubleprecision (dim, dim) spectral density as output
        doubleprecision, intent(inout) :: gammas(dim, dim)
        !> double constant to be put in the spectra density
        doubleprecision, intent(in) :: gamma
        !> Doubleprecision (dim) array of eigenvectors
        doubleprecision, intent(in) :: eigenvalues(dim)

        integer i,j
        doubleprecision wastedouble1
        gammas = 0.d0
        do i=1, dim
            do j=1, dim
                wastedouble1 = eigenvalues(i) - eigenvalues(j)
                if((abs(wastedouble1).gt.1.d-9)) then
                    if(eigenvalues(i).gt.eigenvalues(j)) gammas(i,j) = gamma
                    if(eigenvalues(j).gt.eigenvalues(i)) gammas(i,j) = -gamma
                end if
            end do
        end do
    end subroutine GammaConstant

    !> Create a Gamma spectral density with real amplitude with intensity proportional to the radiational field
    !> @note the value are in femtosecond
    subroutine GammaRadiationField(gammas, eigenvalues, dim)
        integer, intent(in) :: dim
        doubleprecision, intent(in) :: eigenvalues(dim)
        doubleprecision, intent(inout) :: gammas(dim, dim)
        integer i,j
        doubleprecision hbar_evfs, wastedouble1
        hbar_evfs = hbar_ev*1.d15
        gammas = 0.d0
        do i=1, dim
            do j=1, dim
                wastedouble1 = eigenvalues(i) - eigenvalues(j)
                if(abs(wastedouble1).gt.1.d-14) then
                    gammas(i,j) = 14.4*(wastedouble1/hbar_ev)**3/(3*c_SI**3*1.d30*hbar_evfs)
                end if
            end do
        end do
    end subroutine GammaRadiationField

    !> Create a Gamma spectral density with Ohmic density
    subroutine GammaOhmic(gammas, eigenvalues, eta, omega_c, dim)
        !> dimension of the gamma matrix
        integer, intent(in) :: dim
        !> doubleprecision gamma matrix (dim,dim)
        doubleprecision, intent(inout) :: gammas(dim, dim)
        !> eigenvalues as a doublepreicsion(dim) array
        doubleprecision, intent(in) :: eigenvalues(dim)
        !> doubleprecision eta
        doubleprecision, intent(in) :: eta
        !> doubleprecision omega_c
        doubleprecision, intent(in) :: omega_c

        integer i,j
        doubleprecision hbar_evfs, wastedouble1

        gammas = 0.d0
        do i=1, dim
            do j=1, dim
                wastedouble1 = eigenvalues(i) - eigenvalues(j)
                gammas(i,j) = eta*wastedouble1*e_nepero**(-abs(wastedouble1)/omega_c)
            end do
        end do
    end subroutine GammaOhmic

    !> Having given the difference in energu in eV and the desired distribution, this routine gives the population value
    function choose_density(what, omega, T)
        !> What distribution do you need? BE or BoseEinstein for BoseEinsteinn; dark or darkness for null
        character(len=*), intent(in) :: what
        !> The energy difference as doubleprecision number
        doubleprecision omega
        !> The temperature as a doubleprecision value
        doubleprecision, intent(in) ::T
        !> the output of the function
        doubleprecision choose_density

        if(what.eq.'dark'.or.what.eq.'darkness')then
            choose_density=darkness(omega)
        elseif(what.eq.'BoseEinstein'.or.what.eq.'BE')then
            choose_density = Bose_Einstein(omega, T)
        else
            choose_density = Bose_Einstein(omega, T)
        end if
    end function choose_density
    !> Having given the difference in energy (in eV) and the temperature this subroutine gives the population ratio with the Bose Einstein
    doubleprecision function Bose_Einstein(omega, T)
        doubleprecision omega, T
        Bose_Einstein = 0.d0
        if(omega.ne.0.d0) Bose_Einstein = 1.d0/(dexp(omega/(kB_ev*T))-1)
    end function Bose_Einstein

    !> Having given the difference in energy (in eV) and the temperature this subroutine gives 0
    doubleprecision function darkness(omega)
        doubleprecision omega
        darkness = 0.d0
    end function darkness

    !> @todo deprecated bose einstein!
    doubleprecision function n_distr(a1, a2)
        doubleprecision a1, a2
        if(a1.lt.a2)then
           n_distr=0.d0
        elseif (a1.eq.a2)then
           n_distr=0.d0
        else
           n_distr=e_nepero**(eVtoJ(a1-a2)/(kB_SI*298.15))-1
           n_distr=1/n_distr
        endif
    end function n_distr


    !> Generate the vector in the Liouvillian Space
    !> @note R_ij becomes the element V_{j+(i-1)*dim}
    function ReturnKrylovVector_Complex(vector, dim)
        !> integer dimension of the matrix
        integer, intent(in) :: dim
        !> doublecomplex(dim, dim) matrix
        doublecomplex, intent(in) :: vector(dim, dim)
        !> doublecomplex(dim**2) array in Liouvillian space
        doublecomplex ReturnKrylovVector_Complex(dim**2)
        integer i,j, ij
        do i=1, dim
            do j=1, dim
                ij = j +dim*(i-1)
                ReturnKrylovVector_Complex(ij) = vector(i,j)
            end do
        end do
    end function ReturnKrylovVector_Complex

    !> Generate the vector in the Real space from the Liouvillian space
    !> @note V_{j+(i-1)*dim} becomes the element R_ij
    function ReturnKrylovtoRealVector_Complex(vector, dim)
        !> integer dimension of the vector IN REAL SPACE!
        integer, intent(in) :: dim
        !> doublecomplex Liovullian vector of dimension (dim**2)
        doublecomplex, intent(in) :: vector(dim**2)
        !> doublecomplex vector in real space
        doublecomplex ReturnKrylovtoRealVector_Complex(dim, dim)
        integer i,j, ij
        do i=1, dim
            do j=1, dim
                ij = j +dim*(i-1)
                ReturnKrylovtoRealVector_Complex(i,j) = vector(ij)
            end do
        end do
    end function ReturnKrylovtoRealVector_Complex

    !> This subroutine returns the Liouvillian in the Liouvillian space
    !> The routine needs the eigenvalues and the redfield tensor, so that the rho is written
    !> in the eigenstate of H
    !> @note the eigenvalues need to be inserted in frequency an preferable in femtosecond
    function ReturnKrylovLiouvillian_Complex(vector, val, dim)
        !> integer dimension of the vector in real space
        integer, intent(in) :: dim
        !> Doublecomplex Redfield tensor of dimension (dim, dim, dim, dim)
        doublecomplex, intent(in) :: vector(dim, dim, dim, dim)
        !> doubleprecision eigenvalues frequency of dimension (dim)
        doubleprecision, intent(in) :: val(dim)
        !> doublecomplex Liovullian in Liouvillian space of dimension (dim**2, dim**2)
        doublecomplex ReturnKrylovLiouvillian_Complex(dim**2, dim**2)
        integer i,j,k,l, ij, kl
        ReturnKrylovLiouvillian_Complex = (0.d0, 0.d0)
        do i=1, dim
            do j=1, dim
                ij = j+(i-1)*dim
                ReturnKrylovLiouvillian_Complex(ij,ij) = ReturnKrylovLiouvillian_Complex(ij,ij) -imagine*(val(i)-val(j))
                do k=1, dim
                    do l=1, dim
                        kl = l +dim*(k-1)
                        ReturnKrylovLiouvillian_Complex(ij, kl) =  ReturnKrylovLiouvillian_Complex(ij, kl) + vector(i,j,k,l)
                    end do
                end do
            end do
        end do
    end function ReturnKrylovLiouvillian_Complex

    !> This subroutine returns the Liouvillian in the Liouvillian space
    !> The routine needs the Hamiltonian and the redfield tensor, so that the rho it is not necessary writte in
    !> in the eigenstate of H
    !> @note the commutator need to be inserted in frequency an preferable in femtosecond
    !> @todo IT IS NOT WARKING! SImple reason: the commutato IS NOT independent of time
    function ReturnKrylovLiouvillianComm_Complex(vector, comm, dim)
        !> integer dimension of the matrices
        integer, intent(in) :: dim
        !> doublecomplex redfield tensor of dimension (dim, dim,dim,dim)
        doublecomplex, intent(in) :: vector(dim, dim, dim, dim)
        !> doublecomplex commutator of dimension (dim, dim)
        doublecomplex, intent(in) :: comm(dim, dim)
        !> doublecomplex Liouvillian in Liouville space of dimension (dim**2, dim**2)
        doublecomplex ReturnKrylovLiouvillianComm_Complex(dim**2, dim**2)
        integer i,j,k,l, ij, kl
        ReturnKrylovLiouvillianComm_Complex = (0.d0, 0.d0)
        do i=1, dim
            do j=1, dim
                ij = j+(i-1)*dim
                ReturnKrylovLiouvillianComm_Complex(ij,ij) = ReturnKrylovLiouvillianComm_Complex(ij,ij) -imagine*(comm(i,j))
                do k=1, dim
                    do l=1, dim
                        kl = l +dim*(k-1)
                        ReturnKrylovLiouvillianComm_Complex(ij, kl) =  ReturnKrylovLiouvillianComm_Complex(ij, kl) + vector(i,j,k,l)
                    end do
                end do
            end do
        end do
    end function ReturnKrylovLiouvillianComm_Complex

    !> this subroutine returns the Liouvillian operators defined as:
    !> L = -i/hbar [H, ...] + R:...
    !> to operate in a good way, I insert the commutator inside Redfield employing:
    !> <i|L|j> = sum_kl (R_ijkl -i/hbar \delta_jl H_ik + i/hbar \delta_ki H_lj) <k|...|l>
    !> @note return value in seconds!
    function Liouvillian_Complex(H,red,dim)
        !> integer dimension of the operators
        integer, intent(in) :: dim
        !> doublecomplex (dim,dim) hamiltonian
        doublecomplex, intent(in) :: H(dim, dim)
        !> doublecomplex (dim,dim, dim, dim) redfield operator
        doublecomplex, intent(in) :: red(dim, dim, dim, dim)
        !>Liouvillian as a doublecomplex(dim, dim, dim, dim) matrix
        doublecomplex, allocatable :: Liouvillian_Complex(:,:,:,:)

        integer i,j,k,l
        allocate(Liouvillian_Complex(dim, dim, dim, dim))
        do i=1, dim
            do j=1, dim
                do k=1, dim
                   do l=1, dim
                       Liouvillian_Complex(i,j,k,l) = red(i,j,k,l)
                       if(j.eq.l) Liouvillian_Complex(i,j,k,l) = Liouvillian_Complex(i,j,k,l) - &
                               imagine/hbar_ev * H(i,k)
                       if(i.eq.k) Liouvillian_Complex(i,j,k,l) = Liouvillian_Complex(i,j,k,l) + &
                               imagine/hbar_ev * H(l,j)
                    end do
                end do
            end do
        end do
    end function Liouvillian_Complex

    !> This function return a superoperator into the liouvillian space
    !> @note R_ijkl = V((i-1)*dim+j, (k-1)*dim+l)
    function ReturnSuperOperatorsInLiouvillian(Op, dim)
        !> dimension of the superoperator
        integer, intent(in) :: dim
        !> superoperators as a doublecomplex(dim,dim,dim,dim) matrix
        doublecomplex, intent(in) :: Op(dim,dim,dim,dim)
        !> doublecomplex (dim**2, dim**2) matrix as output
        doublecomplex, allocatable :: ReturnSuperOperatorsInLiouvillian(:,:)
        integer i,j,ij,k,l,kl, dimsq
        dimsq = dim**2
        allocate(ReturnSuperOperatorsInLiouvillian(dimsq, dimsq))
        do i=1, dim
            do j=1, dim
                ij = (i-1)*dim+j
                do k=1, dim
                    do l=1, dim
                        kl = (k-1)*dim + l
                        ReturnSuperOperatorsInLiouvillian(ij,kl) = Op(i,j,k,l)
                    end do
                end do
            end do
        end do
    end function ReturnSuperOperatorsInLiouvillian

    !> This function return an operator into the liouvillian space
    !> @note Op_ij = V((i-1)*dim+j)
    function ReturnOperatorsInLiouvillian(Op, dim)
        !> dimension of the operator
        integer, intent(in) :: dim
        !> superoperators as a doublecomplex(dim,dim) matrix
        doublecomplex, intent(in) :: Op(dim,dim)
        !> doublecomplex (dim**2) vector as output
        doublecomplex, allocatable :: ReturnOperatorsInLiouvillian(:)
        integer i,j,ij,dimsq
        dimsq = dim**2
        allocate(ReturnOperatorsInLiouvillian(dimsq))
        do i=1, dim
            do j=1, dim
                ij = (i-1)*dim+j
                ReturnOperatorsInLiouvillian(ij) = Op(i,j)
            end do
        end do
    end function ReturnOperatorsInLiouvillian

        !> This function return an operator into the liouvillian space
    !> @note Op_ij = V((i-1)*dim+j)
    function ReturnOperatorsInRealSpace(Op, dim)
        !> dimension of the operator
        integer, intent(in) :: dim
        !> superoperators as a doublecomplex(dim,dim) matrix
        doublecomplex, intent(in) :: Op(dim**2)
        !> doublecomplex (dim**2) vector as output
        doublecomplex, allocatable :: ReturnOperatorsInRealSpace(:,:)
        integer i,j,ij,dimsq
        allocate(ReturnOperatorsInRealSpace(dim,dim))
        do i=1, dim
            do j=1, dim
                ij = (i-1)*dim+j
                ReturnOperatorsInRealSpace(i,j) = Op(ij)
            end do
        end do
    end function ReturnOperatorsInRealSpace
    !!!!!!!!!!!!!!!!!!!!!!!!
    !NOW TIME EVOLUTION
    !!!!!!!!!!!!!!!!!!!!!!!
    !> Time evolution with Arnoldi, written by A.D.Phan
    subroutine arnoldi(denvec,redfield,dimsq,nkrylov,dt)
        !> integer dimension of the Liouville space
        integer, intent(in) :: dimsq
        !> doublecomplex density matrix of dimension (dimsq) in the liouvillian space.
        !> In output it will contain the evolved density matrix
        doublecomplex, intent(inout), dimension(dimsq) :: denvec
        !> doublecomplex Liouvillian in the Liouvillian space of dimension (dimsq, dimsq)
        doublecomplex, intent(in), dimension(dimsq,dimsq) :: redfield
        !> integer dimension of the krylov space in which the dynamic is executed
        integer, intent(in) :: nkrylov
        !> doubleprecision timestep
        double precision, intent(in) :: dt
        doublecomplex overlapc
        integer a,b,c,i,j,k,ikr,jkr,kkr,info
        integer, allocatable :: ipiv(:,:)
        double precision, allocatable :: rwork(:)
        doublecomplex,allocatable :: kspace(:,:), hessenberg(:,:),rhoone(:)
        doublecomplex, allocatable :: work(:),wmat(:,:),winv(:,:),lambda(:),coeffs(:),&
                temp1(:),temp2(:,:),vl(:,:)
        allocate(kspace(dimsq,nkrylov), hessenberg(nkrylov,nkrylov),rhoone(dimsq))
        !!! Frobenius norm of the density matrix
        overlapc = (0.d0, 0.d0)
        do a = 1,dimsq
            overlapc = overlapc + dconjg(denvec(a))*denvec(a)
        end do
        denvec = denvec/zsqrt(overlapc)
        !!! First column of the krylov vectors
        kspace = (0.d0, 0.d0)
        do a = 1,dimsq
            kspace(a,1) = denvec(a)
        end do
        !!! BEGIN CONSTRUCTION OF KRYLOV SPACE
        hessenberg = (0.d0, 0.d0)

        do jkr = 1,nkrylov
            !building rhoone
            rhoone = (0.d0,0.d0)
            !$omp parallel do default(none), &
            !$omp private(a,b), &
            !$omp shared(dimsq,redfield,kspace,jkr), &
            !$omp reduction(+:rhoone)
            do a = 1,dimsq
                do b = 1,dimsq
                    rhoone(a) = rhoone(a) + redfield(a,b)*kspace(b,jkr)
                end do
            end do
            !$omp end parallel do

            !writing the updiagonal elements of the heseenberg matrix
            !$omp parallel do default(none), &
            !$omp private(ikr,a), &
            !$omp shared(jkr,dimsq,kspace,rhoone), &
            !$omp reduction(+:hessenberg)
            do ikr = 1,jkr
                do a = 1,dimsq
                    hessenberg(ikr,jkr) = hessenberg(ikr,jkr) + dconjg(kspace(a,ikr))*rhoone(a)
                end do
            end do
            !$omp end parallel do


            !building rhotwo
            !$omp parallel do default(none), &
            !$omp private(ikr,a), &
            !$omp shared(jkr,dimsq,kspace,hessenberg), &
            !$omp reduction(+:rhoone)
            do ikr = 1,jkr
                do a = 1,dimsq
                    rhoone(a) = rhoone(a) - hessenberg(ikr,jkr)*kspace(a,ikr)
                end do
            end do

            !$omp end parallel do

            if (jkr==nkrylov) go to 999

            !writing the lower diagonakspace(a,1)l element of the hessenberg matrix
            do a = 1,dimsq
                hessenberg(jkr+1,jkr) = hessenberg(jkr+1,jkr) + dconjg(rhoone(a))*rhoone(a)
            end do
            hessenberg(jkr+1,jkr) = zsqrt(hessenberg(jkr+1,jkr))

            !writing the jkr+1 krylov space
            do a = 1,dimsq
                if (hessenberg(jkr+1,jkr)==0.d0) then
                    write(*,*) "WARNING! Zero value in Hessenberg matrix!"
                else
                    kspace(a,jkr+1) = rhoone(a) / hessenberg(jkr+1,jkr)
                end if
            end do
            999 continue
        end do
        !!! END CONSTRUCTION OF KRYLOV SPACE
    !    do ikr = 1,nkrylov
    !            write(9999,*) (dreal(hessenberg(ikr,jkr)),jkr=1,nkrylov)
    !    end do
    !    do ikr = 1,nkrylov
    !            write(8888,*) (dimag(hessenberg(ikr,jkr)),jkr=1,nkrylov)
    !    end do
        !!! diagonalization of hessenberg matrix

        allocate(lambda(nkrylov),wmat(nkrylov,nkrylov),winv(nkrylov,nkrylov),vl(nkrylov,nkrylov))
        allocate(work(2*nkrylov),rwork(2*nkrylov))
        call zgeev('N','V',nkrylov,hessenberg,nkrylov,lambda,vl,1,wmat,nkrylov,work,&
                2*nkrylov,rwork,info)
        !write(*,*) 'zgeev info ', info
        winv = wmat !saving wmat. winv, after zgetrf and zgetri, it will contain the inverse of wmat
        deallocate(work,rwork,vl)
        !!! end diagonalization of hessenberg matrix

        !!! matrix inversion of W (here vr)
        allocate(ipiv(nkrylov,nkrylov))
        call zgetrf(nkrylov,nkrylov,winv,nkrylov,ipiv,info)
        !write(*,*) 'zgetrf info', info
        allocate(work(nkrylov))
        call zgetri(nkrylov,winv,nkrylov,ipiv,work,nkrylov,info)
        !write(*,*) 'zgetri info', info
        deallocate(ipiv,work)
        !!! end matrix inversion

        !!! getting <rho_i|rho_0>
        allocate(temp1(nkrylov))
        temp1 = (0.d0, 0.d0)
        !$omp parallel do default(none), &
        !$omp private(ikr,a), &
        !$omp shared(nkrylov,dimsq,kspace), &
        !$omp reduction(+:temp1)
        do ikr = 1,nkrylov
            do a = 1,dimsq
                temp1(ikr) = temp1(ikr) + dconjg(kspace(a,ikr))*kspace(a,1)
            end do
        end do
        !$omp end parallel do
        !!!

        !!! getting W*exp(lambda*dt)*W-1
        allocate(temp2(nkrylov,nkrylov))
        temp2 = (0.d0, 0.d0)
        !$omp parallel do default(none), &
        !$omp private(ikr,jkr,kkr), &
        !$omp shared(nkrylov,wmat,lambda,dt,winv), &
        !$omp reduction(+:temp2)
        do ikr = 1,nkrylov
            do jkr = 1,nkrylov
                do kkr = 1,nkrylov
                    temp2(ikr,jkr) = temp2(ikr,jkr) + wmat(ikr,kkr)*zexp(lambda(kkr)*dt)&
                                                     *winv(kkr,jkr)
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(wmat,lambda,winv)
        !!!

        !!! getting the c_i(t) coefficients
        allocate(coeffs(nkrylov))
        coeffs = (0.d0, 0.d0)
        !$omp parallel do default(none), &
        !$omp private(ikr,jkr), &
        !$omp shared(nkrylov,temp2,temp1), &
        !$omp reduction(+:coeffs)
        do ikr = 1,nkrylov
            do jkr = 1,nkrylov
                coeffs(ikr) = coeffs(ikr) + temp2(ikr,jkr)*temp1(jkr)
            end do
        end do
        !$omp end parallel do
        coeffs = zsqrt(overlapc) * coeffs
        deallocate(temp1,temp2)
        !!!

        !!! writing the density matrix at time t+dt
        denvec = (0.d0, 0.d0)
        !$omp parallel do default(none), &
        !$omp shared(dimsq,nkrylov,coeffs,kspace), &
        !$omp private(a,ikr), &
        !$omp reduction(+:denvec)
        do a = 1,dimsq
            do ikr = 1,nkrylov
                denvec(a) = denvec(a) + coeffs(ikr)*kspace(a,ikr)
            end do
        end do
        !$omp end parallel do
        deallocate(coeffs)
        !!!

    end subroutine arnoldi

    !> This subroutine Diagonalize the Liouvillian
    !> Firstly, I need to diagonalize the liouvillian. However, the obtained vectors are not
    !> bioorthogonal. I then use https://joshuagoings.com/2015/04/03/biorthogonalizing-left-and-right-eigenvectors-the-easy-lazy-way/
    !> and rotate the left and right eigenvectors to obtain bioorthogonal vectors
    !> @note The output consist of the eigenvalues, the right eige(nvectors, and the complex conjugate of the left eigenvectors such that
    !> L r(:,i) = lambda_i r(:,i) and conjg(l(:,i)) L = conjg(l(:,i)) lambda_i. This is because fortran
    subroutine LiovullianDiagonalization(VLeft,VRight, Eigenvalues, Liouvillian, dimsq)
        !> integer dimension of all the matricex and arrays
        integer, intent(in) :: dimsq
        !> doublecomplex right eigenvectors NOT ALLOCATED
        doublecomplex, allocatable, intent(out) :: VRight(:,:)
        !> doublecomplex left eigenvectors NOT ALLOCATED
        doublecomplex, allocatable, intent(out) :: VLeft(:,:)
        !> doublecomplex eigenvalues NOT ALLOCATED
        doublecomplex, allocatable, intent(out) :: Eigenvalues(:)
        !> doublecomplex Liouvillian in Liouville space of dimension (dimsq, dimsq)
        doublecomplex, intent(in) :: Liouvillian(dimsq, dimsq)
        doublecomplex, allocatable :: wasteee(:,:)
        !!!!!!!!!!!!!
        !variable for zgeev
        !!!!!!!!!!!!!
        integer info, i, j
        doublecomplex, allocatable :: work(:), RWORK(:)
        !!!!!!!!!!!!!!!!
        !Variable for decomposition
        !!!!!!!!!!!!!!!!
        doublecomplex, allocatable :: M_tot(:,:), M_L(:,:), M_U(:,:)
        integer, allocatable :: ipiv(:)

        !!!!!!!!!!!!!!!!!!!!
        !Firstly I diagonalize the liouvillian
        !!!!!!!!!!!!!!!!!!!!
        allocate(VRight(dimsq, dimsq),VLeft(dimsq,dimsq), Eigenvalues(dimsq),wasteee(dimsq, dimsq))
        allocate(work(3*dimsq), Rwork(2*dimsq))
        call zgeev('V', 'V', dimsq, Liouvillian, dimsq, &
                Eigenvalues, VLeft, dimsq, VRight, dimsq,&
                work, 3*dimsq, RWORK, info)
        write(*,*) 'Liouvillian diagonalization:', info
        deallocate(work, rwork)
        !!!!!!!!!!!!!!!!!
        !Then I want to decompose the Product of L and R
        !!!!!!!!!!!!!!!!!
        allocate(M_tot(dimsq, dimsq), M_L(dimsq, dimsq), M_U(dimsq, dimsq), ipiv(dimsq))
        call conjgtranspose(VLeft, VLeft, dimsq)
        call npmatmul_complex(M_tot, VLeft, VRight, dimsq, dimsq, dimsq)
        write(*,*) 'M_TOT PRE'
        call write_matrix_complex(M_tot, dimsq, dimsq, 'e7.1')
        call ZGETRF(dimsq, dimsq, M_tot, dimsq, ipiv, info)
        write(*,*) 'inversion info:', info
        write(*,*) 'ipiv:'
        write(*,*) ipiv
        deallocate(ipiv)
        do i=1, dimsq
            do j =1, i-1
                M_L(i,j) = M_tot(i,j)
            end do
            M_L(i,i) = (1.d0, 0.d0)
            do j=i, dimsq
                M_U(i,j) = M_tot(i,j)
            end do
        end do
        write(*,*) 'M_TOT'
        call write_matrix_complex(M_tot, dimsq, dimsq, 'e7.1')
        write(*,*) "M_L"
        call write_matrix_complex(M_L, dimsq, dimsq, 'e7.1')
        write(*,*) "M_U"
        call write_matrix_complex(M_U, dimsq, dimsq, 'e7.1')
        deallocate(M_tot)
        !!!!!!!!!!!!!!!
        ! now I want to invers M_L and M_U
        !!!!!!!!!!!!!!!
        call ZTRTRI('U', 'N', dimsq, M_U, dimsq, info)
        call ZTRTRI('L', 'N', dimsq, M_L, dimsq, info)
        write(*,*) "M_L"
        call write_matrix_complex(M_L, dimsq, dimsq, 'e7.1')
        write(*,*) "M_U"
        call write_matrix_complex(M_U, dimsq, dimsq, 'e7.1')
        !!!!!!!!!!!!!!!!!!!
        !Finally I multiplicate the eigenvectors
        !I want L' = conjg(M_L-1 L).T and R' = RM_U-1
        !!!!!!!!!!!!!!!!!!!
        call npmatmul_complex(VRight, VRight, M_U, dimsq, dimsq, dimsq)
        call npmatmul_complex(VLeft, M_L, VLeft, dimsq, dimsq, dimsq)
        call conjgtranspose(VLeft, VLeft, dimsq)
        deallocate(M_L, M_U)
    end subroutine LiovullianDiagonalization

    !> This routine evolved the density matrix by diagonalizing the Liouvillian and projecting the density matrix
    !> on the right eigenvectors.
    !> @todo apparently this does not work...
    subroutine DiagonalizedLiouvilleEvolution(rho_t, rho_0,LVec, RVec, Value, dimsq,t)
        !> integer dimension of the matrices
        integer, intent(in) :: dimsq
        !> doublecomplex evolved density matrix of dimension (dimsq, dimsq)
        doublecomplex, intent(out) :: rho_t(dimsq)
        !> doublecomplex density matrix before the evolution of dimension (dimsq, dimsq)
        doublecomplex, intent(in) ::  rho_0(dimsq)
        !> doublecomplex left eigenvectors of dimension (dimsq, dimsq)
        doublecomplex, intent(in) :: LVec(dimsq, dimsq)
        !> doublecomplex right eigenvectors of dimension (dimsq, dimsq)
        doublecomplex, intent(in) :: RVec(dimsq, dimsq)
        !> doublecomplex eigenvalues of dimension (dimsq)
        doublecomplex, intent(in) :: Value(dimsq)
        !> doubleprecision final time t
        doubleprecision, intent(in) :: t

        integer i, j, k
        doublecomplex waste
        doublecomplex, allocatable :: waste2(:)
        allocate(waste2(dimsq))

        do k=1, dimsq
         waste2(k) = (0.d0, 0.d0)
         do j=1, dimsq
             waste2(k) = waste2(k) + conjg(Lvec(j,k))*rho_0(j)
         end do
        end do
        do i=1, dimsq
         rho_t(i) = (0.d0, 0.d0)
         waste = (0.d0 ,0.d0)
         do k=1, dimsq
            rho_t(i) = rho_t(i) + Rvec(i,k)*waste2(k)*e_nepero**(Value(k)*t)
         end do
        end do
        deallocate(waste2)
    end subroutine DiagonalizedLiouvilleEvolution

    !> This subroutine save the autocorrelation function of an operator by evolving it with arnoldi algorithm
    !> @note DEPRECATED? I don't know. For now I will not document this subroutine
    subroutine CorrelationFunctionWithGeneratingFunction(op, rho, L_krylov,  dim, dt, timesteps, n_krylov, name)
        integer dim, dimsq, n_krylov, check, timesteps,i
        doubleprecision dt, val(dim)
        doublecomplex, dimension(dim, dim) :: rho, op, GF, op0
        doublecomplex R(dim,dim,dim,dim), L_krylov(dim**2, dim**2), GF_krylov(dim**2),&
                trace
        character(len=*) name

        write(*,*) 'Initialize the calculation of the correlation function:'
        write(*,'(A,F10.3)') 'Deltat (fs) = ', dt
        write(*,'(A,I12)') 'timesteps = ', timesteps
        dimsq = dim**2
        op0 = op
        call commutator_complex(GF, op, rho, dim)
        GF_krylov= ReturnKrylovVector_Complex(GF, dim)

        check = timesteps/20
        write(*,*) 'I will print a control every', check, 'iteration'
        open(9999, file=name)
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !TIME INTEGRATION
        !!!!!!!!!!!!!!!!!!!!!!!!!
        do i=0, timesteps
            GF = ReturnKrylovtoRealVector_Complex(GF_krylov, dim)
            !!!!!!!!!!!!!!!!!!
            !DEBUG
            !!!!!!!!!!!!!!!!!!
!            if (i.lt.20)then
!                write(*,*) '------------------------'
!                write(*,*) 'STEP ', i
!                write(*,*) '------------------------'
!                write(*,*)  'Generating function', i
!                call write_matrix_complex(GF, dim,dim,'e9.2')
!            end if
            call npmatmul_complex(GF, op0, GF, dim, dim, dim)
            !!!!!!!!!!!!!!!!!!
            !DEBUG
            !!!!!!!!!!!!!!!!!!
!            if (i.lt.20)then
!                write(*,*)  'op(0)Omega(t)', i
!                call write_matrix_complex(GF, dim,dim,'e9.2')
!            end if
!            if (i.lt.20)then
!                write(*,*) 'Diag mu0omegat', i
!                write(*,'(<dim_system>(e8.1, 2X))') (real(GeneratingFunction(j,j)),j=1, dim_system )
!                write(*,'(<dim_system>(e8.1, 2X))') (aimag(GeneratingFunction(j,j)),j=1, dim_system )
!            end if
            trace = trace_complex(GF, dim)
            write(9999,*) i*dt, real(trace), aimag(trace)
            if(mod(i,check)==0) then
                write(*,'(A, F11.2, A, e10.3, 2X, e10.3)') 'Trace step',i*dt, ' fs:',real(trace), aimag(trace)
            end if
            call arnoldi(GF_krylov, L_krylov, dimsq, n_krylov, dt)
        end do
        close(9999)
    end subroutine CorrelationFunctionWithGeneratingFunction


    function drho_dt(rho, val, R_ten, dim)
        integer, intent(in) :: dim
        doubleprecision, intent(in) :: R_ten(dim, dim, dim, dim), val(dim)
        doublecomplex, intent(in) :: rho(dim, dim)
        doublecomplex drho_dt(dim, dim)
        integer i,j,k,l,m
        do i=1, dim
            do j=1, dim
                drho_dt(i,j) = -imagine*(val(j)-val(i))*rho(i,j)
                do k=1, dim
                    do l=1, dim
                        drho_dt(i,j) = drho_dt(i,j) + R_ten(i,j,k,l)*rho(k,l)
                    end do
                end do
            end do
        end do
    end function drho_dt

    !> Unitary time evolution with the Runge-Kutta algorithm
    !> @note timestep and val must be inserted in femtosecond
    subroutine Runge_Kutta_Unitary(rho,val, dim,timestep, n_timesteps)
        !> integer dimension of the operators
        integer, intent(in) :: dim
        !> complex operator or complex matrix of dimension(dim, dim) to be evolvedù
        doublecomplex, intent(inout) :: rho(dim,dim)
        !> doubleprecision timestep
        doubleprecision, intent(in) :: timestep
        !> doubleprecision array of dimension (dim) with the eigenvalues
        doubleprecision, intent(in) :: val(dim)
        !> integer number of timesteps before the output
        integer, intent(in) :: n_timesteps
        integer i,j, i_time
        doublecomplex drho(dim,dim), rho1(dim, dim)
        do i_time=1, n_timesteps
            !k1
            drho = drho_dt_unitary(rho, val, dim)
            rho1 = rho + timestep*drho/6
            !k2
            drho = drho_dt_unitary(rho+timestep*drho/2, val, dim)
            rho1 = rho1 + 2*timestep*drho/6
            !k3
            drho = drho_dt_unitary(rho+timestep*drho/2, val, dim)
            rho1 = rho1 + 2*timestep*drho/6
            !k4
            drho = drho_dt_unitary(rho+timestep*drho, val, dim)
            rho1 = rho1 + timestep*drho/6
            !I completed the evolution!
            rho = rho1
        end do
    end subroutine Runge_Kutta_Unitary

    !> Calculate the unitary first derivative of the density matrix using Liouville equation
    !> @note val must be inserted in femtosecond
    function drho_dt_unitary(rho, val, dim)
        !> integer dimension of rho and val
        integer, intent(in) :: dim
        !> doublecomplex operator or complex matrix of dimension(dim, dim) to be evolved
        doublecomplex, intent(in) :: rho(dim, dim)
        !> doubleprecision array of dimension (dim) with the frequency of the eigenvalues
        doubleprecision, intent(in) :: val(dim)
        integer i,j,k,l,m
        doublecomplex drho_dt_unitary(dim, dim)
        do i=1, dim
            do j=1, dim
                drho_dt_unitary(i,j) = -imagine*(val(j)-val(i))*rho(i,j)
            end do
        end do
    end function drho_dt_unitary


    subroutine Runge_Kutta_Redfield_Comm(rho,comm, R_tens, dim,timestep, n_timesteps)
        integer i,j, dim, i_time, n_timesteps
        doubleprecision timestep,R_tens(dim, dim, dim, dim)
        doublecomplex drho(dim,dim), rho(dim,dim), rho1(dim, dim), comm(dim, dim)
        do i_time=1, n_timesteps
            !k1
            drho = drho_dt_comm(rho, comm, R_tens, dim)
            rho1 = rho + timestep*drho/6
            !k2
            drho = drho_dt_comm(rho+timestep*drho/2, comm, R_tens, dim)
            rho1 = rho1 + 2*timestep*drho/6
            !k3
            drho = drho_dt_comm(rho+timestep*drho/2, comm, R_tens, dim)
            rho1 = rho1 + 2*timestep*drho/6
            !k4
            drho = drho_dt_comm(rho+timestep*drho, comm, R_tens, dim)
            rho1 = rho1 + timestep*drho/6
            !I completed the evolution!
            rho = rho1
        end do
    end subroutine Runge_Kutta_Redfield_Comm

    function drho_dt_comm(rho, comm, R_ten, dim)
        integer dim, i,j,k,l,m
        doubleprecision R_ten(dim, dim, dim, dim)
        doublecomplex drho_dt_comm(dim, dim), rho(dim, dim), comm(dim, dim)
        do i=1, dim
            do j=1, dim
                drho_dt_comm(i,j) = -imagine*comm(i,j)
                do k=1, dim
                    do l=1, dim
                        drho_dt_comm(i,j) = drho_dt_comm(i,j) + R_ten(i,j,k,l)*rho(k,l)
                    end do
                end do
            end do
        end do
    end function drho_dt_comm

    subroutine Runge_Kutta_Operator(Op, H, dim,timestep, n_timesteps)
        integer a,b,c,jkr,i,j, dim, i_time, n_timesteps
        doubleprecision timestep, hbar_ev,R_tens(dim, dim, dim, dim), val(dim)
        doublecomplex dop(dim,dim), op(dim,dim), op1(dim, dim), H(dim, dim)
        do i_time=1, n_timesteps
            !k1
            call commutator_complex(dop, Op, H, dim)
            dop = -dop*imagine/hbar_ev
            op1 = op + timestep*dop/6
            !k2
            call commutator_complex(dop, Op+timestep*dop/2, H, dim)
            dop = -dop*imagine/hbar_ev
            op1 = op1 + 2*timestep*dop/6
            !k3
            call commutator_complex(dop, Op+timestep*dop/2, H, dim)
            dop = -dop*imagine/hbar_ev
            op1 = op1 + 2*timestep*dop/6
            !k4
            call commutator_complex(dop, Op+timestep*dop, H, dim)
            dop = -dop*imagine/hbar_ev
            op1 = op1 + timestep*dop/6
            !I completed the evolution!
            op = op1
        end do
    end subroutine Runge_Kutta_Operator

    !>This is the very first runge kutta subroutine. This is the only one that I trust
    !> However, this ruoitine INCORRECTLY do the unitary evolution only. I keep the name
    !> just becouse I do not want to mess up all codes
    !> @note time steps in femtosceond, hamiltonian in eV
    subroutine Runge_Kutta_density(rho,H,dim,timestep)
        !> dimension of the operators
        integer, intent(in) :: dim
        !> doublecomplex (dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim, dim)
        !> doublecomplex (dim, dim) hamitlonian for the unitary evolution
        doublecomplex, intent(in) ::  H(dim, dim)
        !> doubleprecision timestpes in femtosecond!
        doubleprecision, intent(in) :: timestep
        doubleprecision hbar_evfs
        complex*16  rho_n(dim, dim), k1(dim, dim), k2(dim, dim), k3(dim, dim), k4(dim, dim)
        integer i,j
        hbar_evfs = hbar_ev*1.d15
        call commutator_complex(k1, H, rho, dim)
        do i=1, dim
            do j=1, dim
                k1(i,j) = -k1(i,j)*imagine/hbar_evfs
                rho_n(i,j)= rho(i,j) + timestep*k1(i,j)/2
            end do
        end do
        call commutator_complex(k2, H, rho_n, dim)
        do i=1, dim
            do j=1, dim
                k2(i,j) = -k2(i,j)*imagine/hbar_evfs
                rho_n(i,j)= rho(i,j) + timestep*k2(i,j)/2
            end do
        end do
        call commutator_complex(k3, H, rho_n, dim)
        do i=1, dim
            do j=1, dim
                k3(i,j) = -k3(i,j)*imagine/hbar_evfs
                rho_n(i,j)= rho(i,j) + timestep*k3(i,j)
            end do
        end do
        call commutator_complex(k4, H, rho_n, dim)
        do i=1, dim
            do j=1, dim
                k4(i,j) = -k4(i,j)*imagine/hbar_evfs
            end do
        end do
        do i=1, dim
            do j=1, dim
                rho(i,j) = rho(i,j) +(timestep/6.d0)*(k1(i,j)+2*k2(i,j)+2*k3(i,j)+k4(i,j))
            end do
        end do
    end subroutine Runge_Kutta_density

        !>This is the very first runge kutta subroutine. This is the only one that I trust
    !> However, this ruoitine INCORRECTLY do the unitary evolution only. I keep the name
    !> just becouse I do not want to mess up all codes
    !> @note time steps in femtosceond, hamiltonian in eV
    subroutine Runge_Kutta_Liouville(rho,Liouv,dim,timestep)
        !> dimension of the operators
        integer, intent(in) :: dim
        !> doublecomplex (dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim, dim)
        !> doublecomplex (dim, dim, dim, dim) liouvillian
        doublecomplex, intent(in) ::  Liouv(dim, dim, dim, dim)
        !> doubleprecision timestpes in femtosecond!
        doubleprecision, intent(in) :: timestep
        doubleprecision hbar_evfs
        integer i,j, k, l
        doublecomplex, allocatable :: rho_n(:,:), k1(:,:), k2(:,:), k3(:,:), k4(:,:)
        allocate(rho_n(dim, dim), k1(dim, dim), k2(dim, dim), k3(dim, dim), k4(dim, dim))


        do i=1, dim
            do j=1, dim
                k1(i,j) = (0.d0, 0.d0)
                do k=1, dim
                    do l=1, dim
                        k1(i,j) = k1(i,j) + Liouv(i,j,k,l)*rho(k,l)
                    end do
                end do
                rho_n(i,j)= rho(i,j) + timestep*k1(i,j)/2
            end do
        end do

        do i=1, dim
            do j=1, dim
                k2(i,j) = (0.d0, 0.d0)
                do k=1, dim
                    do l=1, dim
                        k2(i,j) = k2(i,j) + Liouv(i,j,k,l)*rho_n(k,l)
                    end do
                end do
                rho_n(i,j)= rho(i,j) + timestep*k2(i,j)/2
            end do
        end do

        do i=1, dim
            do j=1, dim
                k3(i,j) = (0.d0, 0.d0)
                do k=1, dim
                    do l=1, dim
                        k3(i,j) = k3(i,j) + Liouv(i,j,k,l)*rho_n(k,l)
                    end do
                end do
                rho_n(i,j)= rho(i,j) + timestep*k3(i,j)
            end do
        end do

        do i=1, dim
            do j=1, dim
                k4(i,j) = (0.d0, 0.d0)
                do k=1, dim
                    do l=1, dim
                        k4(i,j) = k4(i,j) + Liouv(i,j,k,l)*rho_n(k,l)
                    end do
                end do
            end do
        end do
        do i=1, dim
            do j=1, dim
                rho(i,j) = rho(i,j) +(timestep/6.d0)*(k1(i,j)+2*k2(i,j)+2*k3(i,j)+k4(i,j))
            end do
        end do
        deallocate(rho_n, k1, k2, k3, k4)
    end subroutine Runge_Kutta_Liouville


    !>This is the very first runge kutta subroutine. This is the only one that I trust
    !> @note time steps in femtosceond, hamiltonian in eV
    subroutine Runge_Kutta_redfield(rho,H,red,dim,timestep)
        !> dimension of the operators
        integer, intent(in) :: dim
        !> doublecomplex (dim, dim) density matrix to be evolved
        doublecomplex, intent(inout) :: rho(dim, dim)
        !> doublecomplex (dim, dim) hamitlonian for the unitary evolution
        doublecomplex, intent(in) ::  H(dim, dim)
        !> doublecomplex (dim, dim, dim, dim) redfield tensor
        doublecomplex, intent(in) ::  red(dim, dim, dim, dim)
        !> doubleprecision timestpes in femtosecond!
        doubleprecision, intent(in) :: timestep
        doubleprecision hbar_evfs
        complex*16  rho_n(dim, dim), k1(dim, dim), k2(dim, dim), k3(dim, dim), k4(dim, dim)
        integer i,j,k,l
        hbar_evfs = hbar_ev*1.d15
        call commutator_complex(k1, H, rho, dim)
        do i=1, dim
            do j=1, dim
                k1(i,j) = -k1(i,j)*imagine/hbar_evfs
                do k=1, dim
                    do l=1, dim
                        k1(i,j) = k1(i,j) + red(i,j,k,l)*rho(k,l)
                    end do
                end do
                k1(i,j) = -k1(i,j)*imagine/hbar_evfs
                rho_n(i,j)= rho(i,j) + timestep*k1(i,j)/2
            end do
        end do
        call commutator_complex(k2, H, rho_n, dim)
        do i=1, dim
            do j=1, dim
                k2(i,j) = -k2(i,j)*imagine/hbar_evfs
                do k=1, dim
                    do l=1, dim
                        k2(i,j) = k2(i,j) + red(i,j,k,l)*rho_n(k,l)
                    end do
                end do
                rho_n(i,j)= rho(i,j) + timestep*k2(i,j)/2
            end do
        end do
        call commutator_complex(k3, H, rho_n, dim)
        do i=1, dim
            do j=1, dim
                k3(i,j) = -k3(i,j)*imagine/hbar_evfs
                do k=1, dim
                    do l=1, dim
                        k3(i,j) = k3(i,j) + red(i,j,k,l)*rho_n(k,l)
                    end do
                end do
                rho_n(i,j)= rho(i,j) + timestep*k3(i,j)
            end do
        end do
        call commutator_complex(k4, H, rho_n, dim)
        do i=1, dim
            do j=1, dim
                k4(i,j) = -k4(i,j)*imagine/hbar_evfs
                do k=1, dim
                    do l=1, dim
                        k4(i,j) = k4(i,j) + red(i,j,k,l)*rho_n(k,l)
                    end do
                end do
            end do
        end do
        do i=1, dim
            do j=1, dim
                rho(i,j) = rho(i,j) +(timestep/6.d0)*(k1(i,j)+2*k2(i,j)+2*k3(i,j)+k4(i,j))
            end do
        end do
    end subroutine Runge_Kutta_redfield


subroutine Save_Time_Correlation_unitary(operator, rho_0, val, dim,Deltat, t0, timefinal, n_dt, name)
    integer dim, i, j, timesteps, n_dt
        doubleprecision Deltat, tnow, timefinal, t0,val(dim)
        logical save
        doublecomplex operator(dim, dim), rho_0(dim, dim), operator_0(dim,dim), opp
        character(len=*) name
        open(42, file=name)
        call npmatmul_complex(operator, operator, rho_0, dim, dim, dim)
        operator_0 = operator
        timesteps = int((timefinal-t0)/Deltat)
        do i=0, timesteps, n_dt
            opp = expval_Liouville_complex(operator_0, operator, dim)
            write(42,*) i*Deltat, real(opp), aimag(opp)
            call Runge_Kutta_Unitary(operator, val, dim,Deltat, n_dt)
        end do
        close(42)
    end subroutine Save_Time_Correlation_unitary

    subroutine Save_Time_Correlation(operator, rho_0, H, dim,Deltat, t0, timefinal, n_dt, name)
        integer dim, i, j, timesteps, n_dt
        doubleprecision Deltat, tnow, timefinal, t0
        doublecomplex H(dim, dim)
        logical save
        doublecomplex operator(dim, dim), rho_0(dim, dim), operator_0(dim,dim), opp
        character(len=*) name
        open(42, file=name)
        call npmatmul_complex(operator, operator, rho_0, dim, dim, dim)
        operator_0 = operator
        timesteps = int((timefinal-t0)/Deltat)
        do i=0, timesteps, n_dt
            opp = expval_Liouville_complex(operator_0, operator, dim)
            write(42,*) i*Deltat, real(opp), aimag(opp)
            call Runge_Kutta_Operator(operator, H, dim,Deltat, n_dt)
        end do
        close(42)
    end subroutine Save_Time_Correlation

    !> this subroutine calculate the expectation value of an operator by multiplying it
    !> for the density matrix and tracing over the basis
    doublecomplex function expval_Liouville_complex(A, rho, dim)
        implicit none
        !> integer dimension of the operators
        integer, intent(in) :: dim
        !> doublecomplex operator of dimension (dim,dim) to be evaluated
        doublecomplex, intent(in) :: A(dim,dim)
        !> doublecomplex density matrix of dimension (dim,dim)
        doublecomplex, intent(in) :: rho(dim,dim)
        doublecomplex, allocatable :: token1(:,:)
        allocate(token1(dim, dim))
        call npmatmul_complex(token1, A, rho, dim, dim, dim)
        expval_Liouville_complex = trace_complex(token1, dim)
        deallocate(token1)
    end function expval_Liouville_complex


end module Redfield
