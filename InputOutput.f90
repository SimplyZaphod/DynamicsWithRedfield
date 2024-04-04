module usefullIO
  implicit none

contains
  subroutine write_matrix_double(mat, a1, a2, fin)
    implicit none
    !> dimension of th matrix
    integer a1,a2
    !> matrix to be printed
    doubleprecision mat(a1, a2)
    !> optional format (default f8.3)
    character(len=*), optional :: fin
    character*50 f
    integer i1,i2
    if(present(fin)) then
       write(f,*) '('//fin//',A2)'
    else
       f = '(f8.3, A2)'
    end if
    write(*,*)
    write(*,*)
    do i1=1, a1
       do i2=1, a2
          write(*,f, advance='no') mat(i1, i2)
       enddo
       write(*,*) ''
    enddo
    write(*,*)
    write(*,*)
  end subroutine write_matrix_double
  subroutine write_int_mat(mat, a1, a2)
    implicit none
    integer i1,i2,a1,a2
    integer mat(a1, a2)
    
    do i1=1, a1
       do i2=1, a2
          write(*,'(I3, A1)', advance='no') mat(i1, i2)
       enddo
       write(*,*) ''
    enddo
  end subroutine write_int_mat
    
    subroutine write_matrix_complex(mat, a1, a2, fin)
      integer i,j,a1,a2
      complex*16 mat(a1, a2)
      character(len=50) f
      character(len=*), optional :: fin
      if(present(fin))then
         write(f,*) '(A1,'//fin//',A1,'//fin//',A1)'
      else
         f='(A1,f5.2, A1, f5.2, A1)'
      end if
      do i=1, a1
         do j=1, a2
           write(*,f, advance='no') '(',real(mat(i,j)),',',aimag(mat(i,j)),')'
         enddo
         write(*,*) ''
      enddo
      write(*,*)
    end subroutine write_matrix_complex
    subroutine save_matrix_complex(filename,mat, a1, a2, fin)
      character(len=*) filename
      integer i,j,a1,a2
      complex*16 mat(a1, a2)
      character(len=50) f
      character(len=*), optional :: fin
      if(present(fin))then
         write(f,*) '(A1,'//fin//',A1,'//fin//',A1)'
      else
         f='(A1,f5.2, A1, f5.2, A1)'
      end if
      open(9213, file=filename)
      write(*,*) 'Saved at: ', filename
      do i=1, a1
         do j=1, a2
           write(9213,f, advance='no') '(',real(mat(i,j)),',',aimag(mat(i,j)),')'
         enddo
         write(*,*) ''
      enddo
      close(9213)
      
    end subroutine save_matrix_complex
    
    
    subroutine write_vector_double(vec, a1, fin)
      integer i,a1
      doubleprecision vec(a1)
      character(len=50) f
      character(len=*), optional :: fin
      if(present(fin))then
         write(f,*) '('//fin//')'
      else
         f='(F7.3)'
      end if
      do i=1, a1
         write(*,f) vec(i)
      enddo
      write(*,*) ''
      
    end subroutine write_vector_double
    
    
        subroutine read_xyz(xyz_file_name, atom_list_complete, xyz_atom_coordinates, xyz_coordinates_dimension)
      integer xyz_coordinates_dimension, i
      character(len=1), allocatable ::  atom_list_complete(:)
      doubleprecision, allocatable :: xyz_atom_coordinates(:,:)
      character(len=*) xyz_file_name
      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '     MODULE: usefullIO, subroutine read_xyz            '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*) "Let's read some files!"

      open(6589,file=xyz_file_name)
      write(*,*) 'xyz file: '//xyz_file_name
      read(6589,*) xyz_coordinates_dimension
      write(*,define_a_format(xyz_coordinates_dimension)) xyz_coordinates_dimension
      write(*,*)
      read(6589,*)
      allocate(atom_list_complete(xyz_coordinates_dimension))
      allocate(xyz_atom_coordinates(xyz_coordinates_dimension, 3))
      do i = 1, xyz_coordinates_dimension
         read(6589, *) atom_list_complete(i),  xyz_atom_coordinates(i, 1), xyz_atom_coordinates(i, 2), xyz_atom_coordinates(i, 3)
         write(*,*)  atom_list_complete(i),  xyz_atom_coordinates(i, 1), xyz_atom_coordinates(i, 2), xyz_atom_coordinates(i, 3)
      enddo
      close(6589)
      write(*,*)
      write(*,*) '##########   End of file!    ##########'
    endsubroutine read_xyz


    subroutine write_xyz(xyz_file_name, atom_list_complete, xyz_atom_coordinates, xyz_coordinates_dimension)
      integer xyz_coordinates_dimension, i
      character(len=1)  atom_list_complete(xyz_coordinates_dimension)
      doubleprecision xyz_atom_coordinates(xyz_coordinates_dimension,3)
      character*100 xyz_file_name
      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '      MODULE: usefullIO, subroutine write_xyz       '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*)


 8194 format (A1,A3, F10.5, A3, F10.5, A3, F10.5)
      open(6589,file=xyz_file_name)
      write(6589,define_a_format(xyz_coordinates_dimension))  xyz_coordinates_dimension
      write(6589,'(A37)') 'Generated within the module usefullIO'
      do i = 1, xyz_coordinates_dimension
         write(6589,8194)  atom_list_complete(i),'   ',  xyz_atom_coordinates(i, 1),'   ', xyz_atom_coordinates(i, 2),'   ', xyz_atom_coordinates(i, 3)
      enddo
      close(6589)
      write(*,*) 'xyz file: WRITTEN!'
    endsubroutine write_xyz

    subroutine print_xyz(xyz_file_name, atom_list_complete, xyz_atom_coordinates, xyz_coordinates_dimension)
      integer xyz_coordinates_dimension, i
      character(len=1)  atom_list_complete(xyz_coordinates_dimension)
      doubleprecision xyz_atom_coordinates(xyz_coordinates_dimension,3)
      character*100 xyz_file_name
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '                  ', xyz_file_name,'                   '      
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'


 8194 format (A1,A3, F10.5, A3, F10.5, A3, F10.5)
      write(*,define_a_format(xyz_coordinates_dimension))  xyz_coordinates_dimension
      write(*,'(A37)') ''
      do i = 1, xyz_coordinates_dimension
         write(*,8194)  atom_list_complete(i),'   ',  xyz_atom_coordinates(i, 1),'   ', xyz_atom_coordinates(i, 2),'   ', xyz_atom_coordinates(i, 3)
      enddo
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      
    endsubroutine print_xyz


    subroutine save_vector_integer(filename, vector, vector_dimension)
      character(len=*) filename
      integer vector_dimension, i
      integer vector(vector_dimension)
      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) ' MODULE: usefullIO, subroutine save_vector_integer  '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*)
      open(9213, file=filename)
      write(*,*) 'Saved at: ', filename
      do i=1, vector_dimension
         write(9213,define_a_format(vector(i))) vector(i)
      enddo
      close(9213)
      
    endsubroutine save_vector_integer

    
    subroutine save_matrix_double(filename, mat, mat_dimension)
      character(len=*) filename

      integer mat_dimension, i, j
      doubleprecision mat(mat_dimension, mat_dimension)

      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '   MODULE: usefullIO, subroutine save_mat_double    '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*)
      open(9213, file=filename)
      write(*,*) 'Saved at: ', filename
      do i=1, mat_dimension
         do j=1, mat_dimension
            write(9213, '(F10.3, A1)', advance='no') mat(i, j)
         enddo
         write(9213, *) ''
      enddo
      close(9213)
      
    endsubroutine save_matrix_double

subroutine save_matrix_integer(filename, mat, mat_dimension1, mat_dimension2)
      character(len=*) filename

      integer mat_dimension1,mat_dimension2, i, j
      integer mat(mat_dimension1, mat_dimension2)

      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '    MODULE: usefullIO, subroutine save_mat_integer     '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*)
      open(9213, file=filename)
      write(*,*) 'Saved at: ', filename
      do i=1, mat_dimension1
         do j=1, mat_dimension2
            write(9213, '(I3, A1)', advance='no') mat(i, j)
         enddo
         write(9213, *) ''
      enddo
      close(9213)
      
    endsubroutine save_matrix_integer

    subroutine save_vector_double(filename, vector, vector_dimension)
      character(len=*) filename
      integer vector_dimension, i
      doubleprecision vector(vector_dimension)
      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) ' MODULE: usefullIO, subroutine save_vector_double   '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*)
      open(9213, file=filename)
      do i=1, vector_dimension
         write(9213, '(F15.5)') vector(i)
      enddo
      close(9213)
      
    endsubroutine save_vector_double
    
    subroutine save_vector_binary(filename, vector, vector_dimension, formatt)
      character(len=*) filename
      character*100 formatt
      integer vector_dimension, i
      integer vector(vector_dimension)
      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) ' MODULE: usefullIO, subroutine save_vector_integer  '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*)
      open(9213, file=filename)
      write(*,*) 'Saved at: ', filename
      do i=1, vector_dimension
         write(9213,formatt) vector(i)
      enddo
      close(9213)
      
    endsubroutine save_vector_binary
    
     subroutine write_vector_binary(vector, vector_dimension, formatt)
      character*100 formatt
      integer vector_dimension, i
      integer vector(vector_dimension)

      write(*,*)
      write(*,*)
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '   MODULE: usefullIO, subroutine write_vector_binary   '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      write(*,*)
      do i=1, vector_dimension
         write(*,formatt) vector(i)
      enddo
      
      
    endsubroutine write_vector_binary
    function define_a_format(num)
      
      integer num, i
      character*4 define_a_format
      character*1 ichar
      
      i = 0
      
      do while(.true.)
         if(num .lt. 10**i) exit
         i = i + 1
      end do
      write(ichar,'(I1)') i
      define_a_format = '(I'//ichar//')'

    end function define_a_format



    subroutine write_vector_complex(vec, a1)
      integer i,a1
      complex*16 vec(a1)
      
      do i=1, a1
         write(*,'(A1,2f5.2,A2)', advance='no') '(',vec(i),')'
      enddo
      write(*,*) ''
      
    end subroutine write_vector_complex



    
    subroutine StripSpaces(string)				!thank you stackoverflow!
      character(len=*) :: string
      integer :: stringLen
      integer :: last, actual
      
      stringLen = len (string)
      last = 1
      actual = 1
      
      do while (actual < stringLen)
         if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
         else
            last = last + 1
            if (actual < last) &
                 actual = last
         endif
      end do
      
    end subroutine StripSpaces
    
 
end module usefullIO
