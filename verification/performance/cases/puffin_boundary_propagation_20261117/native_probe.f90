! Invoke the frozen production diffraction/boundary routines on source-free fields.
! The production archive and its Fortran sources are never rebuilt or changed.
program native_probe
  use mpi
  use puffin_kinds, only: wp, ip
  use puffin_constants, only: pi
  use Globals
  use puffin_mpiInfo, only: tProcInfo_G
  use puffin_fftwInfo, only: tTransInfo_G
  use transforms, only: getTransformPlans4FEL,GetKValues,clearTransformPlans
  use PDiff, only: DiffractionStep
  use GlobalTypes, only: tSimulationContext
  implicit none
  integer :: ierr,rank,nproc,ix,iy,iz,izg,ind,n,steps,stride,nloc,nt,uh,uf,ul,um,ur,frames
  integer :: ix0,iy0,ic,jc,mx,my
  integer(ip),allocatable :: recvs(:),displs(:)
  real(wp),allocatable :: ar(:),ai(:),oldr(:),oldi(:),line(:),xs(:),ys(:)
  complex(wp),allocatable :: crop(:,:)
  type(tSimulationContext) :: ctx
  real(wp) :: beta,h,lg,lc,lt,dx,dy,dz2,x,y,ph,amp,local_u,global_u,u0,roi,cellw
  real(wp) :: theta,kx_phys,kabs,lambda_phys,delta,local_err,global_err,local_max,global_max
  complex(wp) :: factor,actual,expect
  logical :: ok,plane
  character(1024) :: arg,out,mode,hostname
  integer :: hostlen

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  tProcInfo_G%comm=MPI_COMM_WORLD;tProcInfo_G%rank=rank
  tProcInfo_G%size=nproc;tProcInfo_G%qroot=rank==0
  call get_command_argument(1,out)
  call get_command_argument(2,arg);read(arg,*) NX_G
  call get_command_argument(3,arg);read(arg,*) NY_G
  call get_command_argument(4,arg);read(arg,*) beta
  call get_command_argument(5,arg);read(arg,*) steps
  call get_command_argument(6,arg);read(arg,*) stride
  call get_command_argument(7,mode);plane=trim(mode)=='plane'
  call get_command_argument(8,arg);read(arg,*) lg
  call get_command_argument(9,arg);read(arg,*) lc
  call get_command_argument(10,arg);read(arg,*) lt
  if (NX_G<81.or.NY_G<81.or.mod(NX_G,2)/=1.or.mod(NY_G,2)/=1) call MPI_Abort(MPI_COMM_WORLD,11,ierr)
  NZ2_G=21;NBX_G=16;NBY_G=16;NBZ2_G=37
  fieldMesh=iPeriodic;sperwaves_G=1.0_wp
  qDiffraction_G=.true.;qFilter=.false.;sfilt=0.0_wp;sBeta_G=beta;ioutInfo_G=0
  dx=40e-6_wp;dy=dx
  sLengthOfElmX_G=dx/lt;sLengthOfElmY_G=dy/lt
  dz2=4.0_wp*pi*.005_wp/20.0_wp;sLengthOfElmZ2_G=dz2
  h=.25_wp*4.0_wp*pi*.005_wp
  ffact=real(NX_G,wp)*real(NY_G,wp)*20.0_wp
  ctx%integration%time_start=0.0_wp
  call getTransformPlans4FEL([NX_G,NY_G,NZ2_G],.false.,ok)
  if (.not.ok) call MPI_Abort(MPI_COMM_WORLD,12,ierr)
  nloc=tTransInfo_G%loc_nz2;nt=NX_G*NY_G
  allocate(kx_G(0:NX_G-1),ky_G(0:NY_G-1),kz2_loc_G(0:nloc-1))
  allocate(recvs(nproc),displs(nproc));call GetKValues(recvs,displs,ok)
  if (.not.ok.or.nloc<=0) call MPI_Abort(MPI_COMM_WORLD,13,ierr)
  allocate(ar(nt*nloc),ai(nt*nloc),oldr(nt*nloc),oldi(nt*nloc))
  allocate(crop(81,81),line(81),xs(NX_G),ys(NY_G))
  do ix=1,NX_G;xs(ix)=real(ix-(NX_G+1)/2,wp)*dx;end do
  do iy=1,NY_G;ys(iy)=real(iy-(NY_G+1)/2,wp)*dy;end do
  ix0=(NX_G-81)/2+1;iy0=(NY_G-81)/2+1
  kabs=(2.0_wp*acos(-1.0_wp)/(20.0_wp*dz2))/lc
  kx_phys=30000.0_wp;theta=kx_phys/kabs;lambda_phys=2*acos(-1.0_wp)/kabs
  mx=8;my=3
  do iz=1,nloc
    izg=tTransInfo_G%loc_z2_start+iz-1
    do iy=1,NY_G
      y=ys(iy)
      do ix=1,NX_G
        x=xs(ix);ind=ix+NX_G*(iy-1)+nt*(iz-1)
        if (plane) then
          amp=1.0_wp
          ph=2*acos(-1.0_wp)*(real(mx*(ix-1),wp)/NX_G+real(my*(iy-1),wp)/NY_G-real(izg,wp)/20)
        else
          amp=exp(-((x+300e-6_wp)**2+y*y)/(4*(200e-6_wp)**2))
          ph=kx_phys*x-2*acos(-1.0_wp)*real(izg,wp)/20
        end if
        ar(ind)=amp*cos(ph);ai(ind)=amp*sin(ph)
      end do
    end do
  end do
  oldr=ar;oldi=ai
  call MPI_Get_processor_name(hostname,hostlen,ierr)
  print '(a,i0,a,a,a,i0,a,i0)', 'BIND rank=',rank,' host=',trim(hostname),' slab_start=',tTransInfo_G%loc_z2_start,' slab_count=',nloc
  if (rank==0) then
    open(newunit=uh,file=trim(out)//'/history.csv',status='new')
    write(uh,'(a)') 'step,z_m,total_norm2_m2,roi_norm2_m2'
    open(newunit=uf,file=trim(out)//'/field_xy.bin',access='stream',form='unformatted',status='new')
    open(newunit=ul,file=trim(out)//'/line_x.bin',access='stream',form='unformatted',status='new')
    open(newunit=um,file=trim(out)//'/META.json',status='new')
    write(um,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a,es25.16e3,a)') &
      '{"nx":',NX_G,',"ny":',NY_G,',"steps":',steps,',"snapshot_stride":',stride,',"mpi_size":',nproc, &
      ',"dx_m":',dx,',"h_bar":',h,',"step_m":',h*lg,',"kx_per_m":',kx_phys,',"theta_x_rad":',theta,'}'
    close(um)
    print *, 'STARTUP native DiffractionStep, plane=',plane,' beta=',beta,' theta=',theta,' lambda=',lambda_phys
    flush(6)
  end if
  frames=0
  do n=0,steps
    if (n>0) then
      ctx%lattice%cumulative_steps=n
      call DiffractionStep(h,ar,ai,ctx,ok)
      if (.not.ok.or.any(ar/=ar).or.any(ai/=ai)) call MPI_Abort(MPI_COMM_WORLD,14,ierr)
    end if
    local_u=sum(ar*ar+ai*ai)*dx*dy/20.0_wp
    call MPI_Allreduce(local_u,global_u,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(n==0)u0=global_u
    if(rank==0)then
      if(tTransInfo_G%loc_z2_start/=0)call MPI_Abort(MPI_COMM_WORLD,15,ierr)
      do jc=1,81
        iy=iy0+jc-1
        do ic=1,81
          ix=ix0+ic-1;ind=ix+NX_G*(iy-1)
          crop(ic,jc)=cmplx(ar(ind),ai(ind),wp)
        end do
      end do
      roi=0.0_wp
      do jc=28,54
        do ic=28,54
          cellw=1.0_wp
          if(ic==28.or.ic==54)cellw=cellw*.5_wp
          if(jc==28.or.jc==54)cellw=cellw*.5_wp
          roi=roi+cellw*abs(crop(ic,jc))**2*dx*dy
        end do
      end do
      write(uh,'(i0,3(a,es26.17e3))') n,',',real(n,wp)*h*lg,',',global_u,',',roi
      line=abs(crop(:,41))**2;write(ul)line
      if(mod(n,stride)==0.or.n==steps)then
        write(uf)crop;frames=frames+1
      end if
      if(mod(n,50)==0.or.n==steps)then
        print '(a,i0,a,f10.5,a,es14.5)', 'PROGRESS step=',n,' z=',real(n,wp)*h*lg,' norm2/U0=',global_u/u0
        flush(6);flush(uh)
      end if
    end if
  end do
  if(plane)then
    delta=h*(kx_G(mx)**2+ky_G(my)**2)/(2*(-2*acos(-1.0_wp)/(20*dz2)))
    factor=exp(cmplx(0.0_wp,steps*delta,wp));local_err=0.0_wp;local_max=0.0_wp
    do ind=1,size(ar)
      actual=cmplx(ar(ind),ai(ind),wp);expect=factor*cmplx(oldr(ind),oldi(ind),wp)
      local_err=max(local_err,abs(actual-expect));local_max=max(local_max,abs(expect))
    end do
    call MPI_Allreduce(local_err,global_err,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_err/global_max>2e-12_wp)call MPI_Abort(MPI_COMM_WORLD,16,ierr)
    if(rank==0)then
      open(newunit=ur,file=trim(out)//'/PLANE_CHECK.json',status='new')
      write(ur,'(a,es26.17e3,a)') '{"complete":true,"native_vs_exact_relative_max":',global_err/global_max,'}'
      close(ur)
    end if
  end if
  if(rank==0)then
    close(uh);close(uf);close(ul)
    open(newunit=ur,file=trim(out)//'/COMPLETE.json',status='new')
    write(ur,'(a,i0,a,i0,a,es26.17e3,a)') '{"complete":true,"steps":',steps,',"frames":',frames,',"initial_norm2_m2":',u0,'}'
    close(ur)
  end if
  call clearTransformPlans(ok)
  call MPI_Finalize(ierr)
end program native_probe
