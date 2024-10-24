module mathfunc
  implicit none
  doubleprecision, parameter :: pi=4.D0*DATAN(1.D0), e_nepero=2.7182818284590d0
  complex*16, parameter :: imagine=(0.d0, 1.d0), reality=(1.d0, 0.d0)
contains
    !> this function return a random integer number that is between a1 and b1
    integer function random_integer(a1, b1)
        !> integer first number of the interval
        integer, intent(in) :: a1
        !> integer second number of the interval
        integer, intent(in) :: b1

        integer i
        doubleprecision u

        call random_seed()
        call random_number(u)
        random_integer = a1 + FLOOR((b1+1-a1)*u)

    end function random_integer

    !> this function return a random doubleprecision number having given an optional dispersion and center
    doubleprecision function random_stdnormal_double(sigma, mu)
        doubleprecision u1, u2, sigma1, mu1
        !> doubleprecision optional dispersion of the gaussian
        doubleprecision, optional, intent(in) :: sigma
        !> doubleprecision optional centre of the gaussian
        doubleprecision, optional, intent(in) :: mu
        if(present(sigma))then
        sigma1 = sigma
        else
            sigma1 = 1.d0
        end if
        if(present(mu))then
            mu1 = mu
        else
            mu1 = 0.d0
        end if
        call random_seed()
        call random_number(u1)
        call random_number(u2)

        random_stdnormal_double = sigma1*dsqrt(-2*dlog(u1))*dcos(2*pi*u2) + mu1

    end function random_stdnormal_double

    !> this function compute the doubleprecision distance between 2 3X1 array
    function distance3D(a, b)
        doubleprecision distance3D, a(3), b(3)
        distance3D = sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)
        endfunction distance3D

        doubleprecision  function npdot_double(a1, a2, dim)
        integer i, dim
        doubleprecision a1(dim), a2(dim)
        npdot_double = 0.d0
        do i=1, dim
            npdot_double = npdot_double +a1(i)*a2(i)
        end do

    end function npdot_double

    !>this function compute the dot product between two doublecomplex (dim) array
    doubleprecision  function npdot_complex(a1, a2, dim)
        integer i, dim
        doublecomplex a1(dim), a2(dim)
        npdot_complex = (0.d0, 0.d0)
        do i=1, 3
            npdot_complex = npdot_complex +a1(i)*dconjg(a2(i))
        end do
    end function npdot_complex

    function Extract_Vector_double(vector, pos, len, dim)
        integer pos, i, len, dim
        doubleprecision Extract_Vector_double(dim), vector(len,dim)

        do i = 1, dim
            Extract_vector_double(i) = vector(pos, i)
        enddo
    end function Extract_Vector_double

    !>this function compute the factorial of an input integer number
    integer function factorial(num)
        integer num, i

        factorial = 1
        if(num.eq.0)then
            return
        end if

        do i=1, num
            factorial = factorial*i
        end do
    end function factorial

    !>this function compute the binomial (n m) between two integer n and m numbers
    integer function binomial(n, m)
        integer n,m, i,count1, count2, count3
        count1=1
        do i=n-m+1, n
            count1 = count1 * i
        end do
        count2 = factorial(m)
        binomial = count1/count2
    end function binomial
  
    function mod_positive(a, b)
        integer a, b, mod_positive
        mod_positive = a
        do while (mod_positive.lt.0)
           mod_positive = mod_positive + b
        enddo
        do while (mod_positive.ge.b)
           mod_positive = mod_positive - b
        enddo
    endfunction mod_positive
  
    subroutine assign_Vector_Double(vec_to_return,vector, pos, len, dim)
        integer pos, i, len, dim
        doubleprecision vec_to_return(dim), vector(len,dim)

        do i = 1, dim
            vec_to_return(i) = vector(pos, i)
        enddo
    endsubroutine assign_Vector_Double

  
    function npnorm(vec, dim)
        integer i, dim
        doubleprecision npnorm, vec(dim)
        npnorm=0.d0
        do i=1, dim
            npnorm = npnorm + vec(i)**2
        enddo
        npnorm = sqrt(npnorm)
    endfunction npnorm

    function npnorm_complex(vec, dim)
        integer i, dim
        doublecomplex npnorm_complex, vec(dim)
        npnorm_complex=(0.d0, 0.d0)
        do i=1, dim
            npnorm_complex = npnorm_complex + conjg(vec(i))*vec(i)
        enddo
        npnorm_complex = sqrt(npnorm_complex)
    endfunction npnorm_complex

    !> this matrix rotate a double matrix M using a double unitary rotation U
    !> @note parallelized with openmp!
    subroutine rotate_matrix_double(mat, U, dim_base)
        !> integer dimension of the matrix to be rotated
        integer, intent(in) :: dim_base
        !> doubleprecision (dim,dim) matrix of rotation
        doubleprecision, intent(in) :: U(dim_base, dim_base)
        !> doubleprecision (dim, dim) matrix that need to be rotated
        doubleprecision, intent(inout) :: mat(dim_base, dim_base)

        doubleprecision, allocatable :: mat_rot(:,:)
        integer i, j, k, l

        allocate(mat_rot(dim_base, dim_base))
        !$omp parallel do default(none), &
        !$omp private(i,j), &
        !$omp shared(mat, mat_rot, dim_base)
        do i=1, dim_base
            do j=1, dim_base
                mat_rot(i,j) = mat(i,j)
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do default(none), &
        !$omp private(i,j,k,l), &
        !$omp shared(mat_rot, U, mat, dim_base)
        do i=1, dim_base
            do j=1, dim_base
                do k=1, dim_base
                    do l=1, dim_base
                        mat(i,j) = U(k,i)*U(l,j)*mat_rot(k,l)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do

        deallocate(mat_rot)
    end subroutine rotate_matrix_double

    !> this matrix rotate a doublecomplex matrix M using a doublecomplex unitary rotation U
    !> @note parallelized with openmp!
    subroutine rotate_matrix_complex(mat, U, dim_base)
        !> integer dimension of the matrix to be rotated
        integer, intent(in) :: dim_base
        !> doubleprecision (dim,dim) matrix of rotation
        doublecomplex, intent(in) :: U(dim_base, dim_base)
        !> doubleprecision (dim, dim) matrix that need to be rotated
        doublecomplex, intent(inout) :: mat(dim_base, dim_base)

        doublecomplex, allocatable :: mat_rot(:,:)
        integer i, j, k, l

        allocate(mat_rot(dim_base, dim_base))
        call npmatmul_complex(mat, mat, U, dim_base, dim_base, dim_base)
        call conjgtranspose(mat_rot, U, dim_base)
        call npmatmul_complex(mat, mat_rot, mat, dim_base, dim_base, dim_base)

    end subroutine rotate_matrix_complex

    !>this subroutine compute the matrix product between 3 doublecomplex matrix
    !> @note parallelized with openmp
    subroutine npmatmul_complex(prod,h1, h2, a, b, c)
        integer, intent(in) :: a, b, c
        integer i, j, k
        doublecomplex :: h1(a,b), h2(b,c)
        doublecomplex, intent(inout) ::prod(a,c)
        doublecomplex, allocatable :: prod1(:,:)
        allocate(prod1(a,c))
        prod1 = (0.d0, 0.d0)
        !$OMP parallel do default(none), &
        !$OMP private(i,j,k), &
        !$OMP shared(h1,h2,a,b,c),&
        !$OMP reduction(+:prod1)
        do i=1, a
            do j=1, c
                do k=1, b
                    prod1(i,j) = prod1(i,j) + h1(i, k)*h2(k,j)
                enddo
            enddo
        enddo
        !$OMP end parallel do

        !$OMP parallel do default(none), &
        !$OMP private(i,j), &
        !$OMP shared(prod1, prod,a,c)
        do i=1, a
            do j=1, c
                prod(i,j) = prod1(i,j)
            end do
        end do
        !$OMP end parallel do
    endsubroutine npmatmul_complex

    function direct_product_complex(A,B,dimA,dimB)
        !> integer dimension of matrix A and B
        integer, intent(in) :: dimA, dimB
        !> doublecomplex matrix A and B
        doublecomplex, intent(in) :: A(dimA,dimA), B(dimB,dimB)
        !> direct product output AXB
        doublecomplex direct_product_complex(dimA*dimB, dimA*dimB)
        integer dim, iA, jA, iB, jB, iC, jC
        if(dimA.eq.1)then
            direct_product_complex = B
        elseif(dimB.eq.1) then
            direct_product_complex = A
        else
            do iA = 1, dimA
                do jA=1, dimA
                    do iB=1, dimB
                        iC = iB+(iA-1)*dimB
                        do jb=1, dimB
                            jC = jB+(jA-1)*dimB
                            direct_product_complex(iC, jC) = A(iA, jA)*B(iB, jB)
                        end do
                    end do
                end do
            end do
        end if
    end function direct_product_complex

    function npones(dim)
        !> dimension of the output
        integer, intent(in) :: dim
        !> integer (dim) matrix of ones on the diagonal
        integer, allocatable :: npones(:,:)
        integer i
        allocate(npones(dim, dim))
        npones = 0
        do i=1, dim
            npones(i,i) = 1
        end do
    end function npones

    subroutine conjgtranspose(res, mat, dim)
        integer, intent(in) :: dim
        doublecomplex, intent(inout) :: res(dim, dim)
        doublecomplex, intent(in) :: mat(dim,dim)
        integer i, j
        doublecomplex, dimension(:,:), allocatable :: aus(:,:)
        allocate(aus(dim, dim))

        do i=1, dim
            do j=1, dim
                aus(j,i) = dconjg(mat(i,j))
            end do
        end do

        do i=1, dim
            do j=1, dim
                res(i,j) = aus(i,j)
            end do
        end do
    end subroutine conjgtranspose
  
    function exp_value_complex(vec1, mat, vec2, dim)
        integer, intent(in) :: dim
        doublecomplex, intent(in) :: vec1(dim), mat(dim, dim), vec2(dim)

        doublecomplex, dimension(:), allocatable :: vec3(:)
        doublecomplex exp_value_complex

        integer  i, j

        allocate(vec3(dim))
        vec3=(0.d0, 0.d0)
        exp_value_complex = (0.d0, 0.d0)
        do i=1, dim
            do j=1, dim
                vec3(i)=vec3(i)+mat(i,j)*vec2(j)
            enddo
        end do

        do j=1, dim
            exp_value_complex=exp_value_complex+dconjg(vec1(j))*vec3(j)
        enddo
    endfunction exp_value_complex

    subroutine npcopy_matrix_complex(mat1, mat2, n1, n2)
    integer n1, n2, i ,j
    complex*16 mat1(n1, n2), mat2(n1,n2)
    do i=1, n1
        do j=1, n2
            mat1(i,j) = mat2(i,j)
        enddo
    enddo
    end subroutine npcopy_matrix_complex

    subroutine npcopy_matrix_double(mat1, mat2, n1, n2)
        integer n1, n2, i ,j
        doubleprecision mat1(n1, n2), mat2(n1,n2)
        do i=1, n1
           do j=1, n2
              mat1(i,j) = mat2(i,j)
           enddo
        enddo
    end subroutine npcopy_matrix_double

    subroutine commutator_complex(comm,mat1, mat2, dim)
        integer, intent(in) ::  dim
        integer i, j
        doublecomplex, dimension(dim, dim):: mat1, mat2
        doublecomplex, dimension(dim, dim), intent(inout) :: comm

        doublecomplex, dimension(:,:), allocatable :: prod1, prod2
        allocate(prod1(dim, dim), prod2(dim, dim))
        call npmatmul_complex(prod1, mat1, mat2, dim, dim, dim)
        call npmatmul_complex(prod2, mat2, mat1, dim, dim, dim)
        comm = prod1 - prod2
        deallocate(prod1, prod2)
    end subroutine commutator_complex

    subroutine nplinspace(vector,start, end, steps)
        integer steps
        double precision start, end, vector(steps), delta, i
        delta = (end-start)/(max(steps-1, 1))
        do i = 0, steps-1
            vector(i+1) = start + delta*i
        enddo
    end subroutine nplinspace

    !> this function computes the trace of a complex operator
    complex*16 function trace_complex(mat, dim)
        integer i
        !> integer dimension of the operator
        integer, intent(in) :: dim
        !> doublecomplex(dim,dim) operator
        complex*16, intent(in) :: mat(dim, dim)

        trace_complex = (0.d0, 0.d0)
        !$OMP parallel do default(none),&
        !$OMP private(i),&
        !$OMP shared(dim, mat), &
        !$omp reduction(+:trace_complex)
        do i=1, dim
            trace_complex = trace_complex + mat(i,i)
        enddo
        !$omp end parallel do
    end function trace_complex

    doubleprecision function trace_real(mat, dim)
        integer i, dim
        doubleprecision mat(dim, dim)
        trace_real = 0.d0
        !$OMP parallel do default(none),&
        !$OMP private(i),&
        !$OMP shared(dim, mat), &
        !$omp reduction(+:trace_real)
        do i=1, dim
            trace_real = trace_real + mat(i,i)
        enddo
        !$omp end parallel do
    end function trace_real
endmodule mathfunc
