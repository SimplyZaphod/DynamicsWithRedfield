program DA
    use mathfunc
    use UsefullIO

    implicit none

    integer i,j, dim_s, dim_v, N_el
    doublecomplex, allocatable :: H(:,:), H_el(:,:), H_v(:,:), &
            Q(:,:), rho(:,:), mu(:,:), &
            waste(:,:)
    doubleprecision z, tau, g, omega
    !!!!!!!!!!
    !DIAG
    !!!!!!!!!!!
    integer info
    doublecomplex, allocatable :: vec(:,:), work(:)
    doubleprecision, allocatable :: values(:), rwork(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PARAMETERS!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tau = 0.95
    z = 0.88
    omega = 0.14
    g = 0.80

    N_el = 2
    write(*,*) 'Number of modes'
    read(*,*) dim_v
    dim_s = n_el * dim_v
    allocate(H(dim_s, dim_s), H_el(N_el, N_el), H_v(dim_v, dim_v),&
    Q(dim_v, dim_v), rho(n_el, n_el))

    H_el = (0.d0, 0.d0)
    H_el(2,2) = 2*z
    H_el(1,2) = -tau
    H_el(2,1) = -tau
    write(*,*) 'H_el'
    call write_matrix_double(real(H_el), N_el, N_el)
    H_v = (0.d0, 0.d0)
    do i=1, dim_v
        H_v(i,i) = (i-1)*omega
     end do
     write(*,*) 'homega n'
    call write_matrix_double(real(H_v), dim_v, dim_v)

    Q = (0.d0, 0.d0)
    do i=1, dim_v -1
        Q(i,i+1) = sqrt((i+1)*1.d0)
        Q(i+1, i) = sqrt((i+1)*1.d0)
     end do
     write(*,*) 'g(a^dag+a)'
    call write_matrix_double(real(Q), dim_v, dim_v)

    rho = (0.d0, 0.d0)
    rho(2,2) = (1.d0, 0.d0)

    H = (0.d0, 0.d0)
    H = direct_product_complex(H_el, npones(dim_v)*reality, N_el, dim_v) &
            - g * direct_product_complex(npones(N_el)*reality, Q, N_el, dim_v) &
            + direct_product_complex(npones(N_el)*reality, H_v, N_el, dim_v)
    write(*,*) 'H_tot'
    call write_matrix_double(real(H), dim_s, dim_s)


    allocate(work(2*dim_s-1), values(dim_s), vec(dim_s, dim_s), rwork(3*dim_s-2))
    vec = H
    call zheev('V', 'U', dim_s, vec,  dim_s, values, work,2*dim_s-1, rwork, info)
    H = (0.d0, 0.d0)
    do i=1, dim_s
        H(i,i) = values(i)
    end do
    open(1, file='Hamiltonian_rot.dat', form='unformatted')
    write(1) H(1:,1:)
    close(1)

    allocate(waste(dim_s, dim_s))
    waste = direct_product_complex(npones(N_el)*reality,Q, N_el, dim_v)
   
    call rotate_matrix_complex(waste, vec, dim_s)
    open(1, file='Q_rot.dat', form='unformatted')
    write(1) waste(1:,1:)
    close(1)

    allocate(mu(dim_s, dim_s))
    mu = direct_product_complex(rho, npones(dim_v)*reality, n_el, dim_v)
    call rotate_matrix_complex(mu, vec, dim_s)
    write(*,*) 'mu_rot'
    call write_matrix_double(real(mu), dim_s, dim_s)
    open(1, file='mu_rot.dat', form='unformatted')
    write(1) mu(1:,1:)
    close(1)
    deallocate(work)
    allocate(work(dim_s))
    work = (0.d0, 0.d0)
    work(1) = reality
    call npmatmul_complex(work, mu, work, dim_s, dim_s, 1)
    write(*,*) 'mu*|G>'
    call write_matrix_double(real(work), dim_s, 1)
    work = work/(npnorm_complex(work, dim_s))
    write(*,*) 'mu*|G>/norm'
    call write_matrix_double(real(work), dim_s, 1)
    waste = (0.d0, 0.d0)
    do i=1, dim_s
        do j=1, dim_s
            waste(i,j) = conjg(work(j))*work(i)
        end do
     end do
     write(*,*) 'rho'
    call write_matrix_double(real(waste), dim_s, dim_s)
    open(1, file='rho_0.dat', form='unformatted')
    write(1) waste(1:,1:)
    close(1)
    waste = (0.d0, 0.d0)
    waste(1,1) = reality
    open(1, file='rho_G.dat', form='unformatted')
    write(1) waste(1:,1:)
    close(1)
    waste = (0.d0, 0.d0)
    waste(dim_s,dim_s) = reality
    open(1, file='rho_E.dat', form='unformatted')
    write(1) waste(1:,1:)
    close(1)  
    
    waste = mu
    do i=1, dim_s
       do j=1, i
          waste(i,j) = (0.d0, 0.d0)
       end do
    end do
    open(1, file='muU.dat', form='unformatted')
    write(1) waste(1:,1:)
    close(1)
                                   
    waste = mu
    do i=1, dim_s
       do j=i, dim_s
          waste(i,j) = (0.d0, 0.d0)
       end do
    end do
    open(1, file='muL.dat', form='unformatted')
    write(1) waste(1:,1:)
    close(1)

end program DA
