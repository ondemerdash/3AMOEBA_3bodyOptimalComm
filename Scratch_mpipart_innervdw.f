             offsetvdw=int((nvdw-1)/numtasks)
             remaindervdw=mod((nvdw-1),numtasks)

             if(taskid.le.remaindervdw-1) then
               start=taskid*offsetvdw+1
               do atomind =start,start+offsetvdw-1
                  call Innerloop_ehal1a(atomind,evt,virevt,devt)
                 do i=1,3
                   do j=1,3
                 virevt_tot(i,j)=virevt_tot(i,j)+virevt(i,j)
                   end do
                 end do

                 do i=1,nvdw
                    do j=1,3
                 devt_tot(j,i)=devt_tot(j,i)+devt(j,i)
                    end do
                 end do

                 evt_tot=evt_tot+evt
               end do

                atomind=numtasks*offsetvdw+taskid+1
                  call Innerloop_ehal1a(atomind,evt,virevt,devt)
                 do i=1,3
                   do j=1,3
                   virevt_tot(i,j)=virevt_tot(i,j)+virevt(i,j)
                   end do
                 end do
                 do i=1,nvdw
                    do j=1,3
                    devt_tot(j,i)=devt_tot(j,i)+devt(j,i)
                    end do
                 end do
                evt_tot=evt_tot+evt

                  call mpi_reduce(evt_tot,ev,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(devt_tot,dev,nvdw*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virevt_tot,virev,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

             else
               start=taskid*offsetvdw+1
               do atomind =start,start+offsetvdw-1
                  call Innerloop_ehal1a(atomind,evt,virevt,devt)

                 do i=1,3
                   do j=1,3
                    virevt_tot(i,j)=virevt_tot(i,j)+virevt(i,j)
                   end do
                 end do

                 do i=1,nvdw
                    do j=1,3
                    devt_tot(j,i)=devt_tot(j,i)+devt(j,i)
                    end do
                 end do

                evt_tot=evt_tot+evt
               end do
                  call mpi_reduce(evt_tot,ev,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(devt_tot,dev,nvdw*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virevt_tot,virev,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if

