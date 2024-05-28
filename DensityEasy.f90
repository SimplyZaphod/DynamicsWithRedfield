program DensityEasy
    implicit none
    integer i, state, dim
    doublecomplex, dimension(:,:), allocatable :: rho
    doubleprecision pop, pop_tot
    character(len=1) con
    
    write(*,*) 'Insert the required dimension'
    read(*,*) dim
    allocate(rho(dim, dim))
    rho = (0.d0, 0.d0)
    pop_tot = 0.d0
10    continue 
    write(*,*) 'Which state do you want?'
    read(*,*) state
    write(*,*) 'Population'
    read(*,*) pop
    pop_tot = pop_tot+pop
    rho(state, state) = (1.d0 ,0.d0)*pop
    write(*,*) 'Do you want to continue [y/n]?'
    read(*,*) con
    if(con.eq.'y'.or.con.eq.'Y') go to 10
    write(*,*) 'Ended!'
    open(1, file='rho.dat', form='unformatted')
    write(1) rho
    close(1)
end program DensityEasy
