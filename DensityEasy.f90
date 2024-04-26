program DensityEasy
    implicit none
    integer i, state, dim
    doublecomplex, dimension(:,:), allocatable :: rho

    write(*,*) 'Insert the required dimension'
    read(*,*) dim
    write(*,*) 'Which state do you want?'
    read(*,*) state
    allocate(rho(dim, dim))
    rho = (0.d0, 0.d0)
    rho(state, state) = (1.d0 ,0.d0)
    open(1, file='rho.dat', form='unformatted')
    write(1) rho
    close(1)
end program DensityEasy
