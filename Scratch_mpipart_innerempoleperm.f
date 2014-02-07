             offsetem=int((npole-1)/numtasks)
             remainderem=mod((npole-1),numtasks)

             if(taskid.le.remainderem-1) then
               start=taskid*offsetem+1
               do atomind =start,start+offsetem-1
                 call Innerloop_ereal1c_3b_Perm(atomind,emt,viremt,demt)
                 do i=1,3
                   do j=1,3
                 viremt_tot(i,j)=viremt_tot(i,j)+viremt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 demt_tot(j,i)=demt_tot(j,i)+demt(j,i)
                    end do
                 end do

                 emt_tot=emt_tot+emt
               end do

                atomind=numtasks*offsetem+taskid+1
                call Innerloop_ereal1c_3b_Perm(atomind,emt,viremt,demt)

                 do i=1,3
                   do j=1,3
                 viremt_tot(i,j)=viremt_tot(i,j)+viremt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 demt_tot(j,i)=demt_tot(j,i)+demt(j,i)
                    end do
                 end do

                 emt_tot=emt_tot+emt

                  call mpi_reduce(emt_tot,emreal,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(demt_tot,demreal,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(viremt_tot,viremreal,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

             else
               start=taskid*offsetem+1
               do atomind =start,start+offsetem-1
                 call Innerloop_ereal1c_3b_Perm(atomind,emt,viremt,demt)

                 do i=1,3
                   do j=1,3
                    viremt_tot(i,j)=viremt_tot(i,j)+viremt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                    demt_tot(j,i)=demt_tot(j,i)+demt(j,i)
                    end do
                 end do

                emt_tot=emt_tot+emt
               end do
                  call mpi_reduce(emt_tot,emreal,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(demt_tot,demreal,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(viremt_tot,viremreal,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if

