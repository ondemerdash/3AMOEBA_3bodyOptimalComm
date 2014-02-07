      call mpi_pack_size(800,mpi_integer,mpi_comm_world,
     &     s1,ierr)

      call mpi_pack_size(1,mpi_integer,mpi_comm_world,
     &     s2,ierr)

      totsize1=numtasks*offset*(s1+s2+2*mpi_bsend_overhead)
      allocate (buf1(totsize1))

      call mpi_buffer_attach(buf1,totsize1,ierr)

      if (taskid.le.remainder-1) then
       call mpi_pack_size(800,mpi_integer,mpi_comm_world,
     &     s3,ierr)
       call mpi_pack_size(1,mpi_integer,mpi_comm_world,
     &     s4,ierr)
       totsize2=remainder*(s3+s4+2*mpi_bsend_overhead)
       allocate (buf2(totsize2))
       call mpi_buffer_attach(buf2,totsize2,ierr)
      end if

      if(taskid.eq.master) then
           do i=1,numtasks
              do j=1,offset
                call mpi_bsend(mollst3mod(1:800,(i-1)*offset+j:
     &           (i-1)*offset+j),800,mpi_integer,i-1,
     &           2*offset*(i-1)+2*(j-1),mpi_comm_world,ierr)
                call mpi_bsend(nmollst3mod( (i-1)*offset+j:
     &           (i-1)*offset+j),1,mpi_integer,i-1,
     &          2*offset*(i-1)+2*(j-1)+1,mpi_comm_world,ierr)
              end do
           end do
           do i=1,remainder
             call mpi_bsend(mollst3mod(1:800,(numtasks*offset+i):
     &       (numtasks*offset+i)),800,mpi_integer,i-1,
     &       2*numtasks*offset+2*(i-1),mpi_comm_world,ierr)
             call mpi_bsend(nmollst3mod((numtasks*offset+i):
     &        (numtasks*offset+i)),1,mpi_integer,i-1,
     &         2*numtasks*offset+2*(i-1)+1,mpi_comm_world,ierr)
           end do
      end if 
