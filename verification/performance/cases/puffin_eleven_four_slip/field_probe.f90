program field_probe11
  use mpi
  use puffin_kinds
  use puffin_constants, only: pi
  use puffin_mpiInfo, only: tProcInfo_G
  use globals, only: qUndEnds_G,qscaled_G,zUndType_G
  use lattice, only: setupMods,initUndulator
  use setupcalcs, only: calcScaling
  use AdapterGlobals, only: PopulateUndulatorFromGlobals
  use GlobalTypes, only: tSimulationContext
  use bfields, only: getBFields
  implicit none
  type(tSimulationContext) :: ctx
  character(len=1024) :: latt,points
  character(len=32) :: typ
  real(kind=wp) :: ux,uy,taper,dz,kx,ky,sz,szl,x(1),y(1),bx(1),by(1),bz(1),z,ex,ey,ez,err,trans
  integer :: nsteps,ierr,ios,n
  call MPI_INIT(ierr)
  tProcInfo_G%comm=MPI_COMM_WORLD
  call MPI_COMM_RANK(MPI_COMM_WORLD,tProcInfo_G%rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,tProcInfo_G%size,ierr)
  tProcInfo_G%qroot=tProcInfo_G%rank==0
  call get_command_argument(1,latt)
  call get_command_argument(2,points)
  typ='aprl-sparc11'
  qscaled_G=.true.
  qUndEnds_G=.true.
  zUndType_G=typ
  ux=0.0_wp;uy=1.0_wp;taper=0.0_wp;kx=0.0_wp;ky=0.0_wp
  call calcScaling(.005_wp,2.060074987322757_wp,298.26088438022981_wp,.028_wp,typ,ux,uy,ctx)
  nsteps=990*80;dz=4*pi*.005_wp/80
  call setupMods(latt,taper,.005_wp,nsteps,dz,ux,uy,kx,ky,ctx%frame)
  sz=0.0_wp
  call initUndulator(1,sz,szl,ctx%frame,ctx%und)
  call PopulateUndulatorFromGlobals(ctx%und)
  if (ctx%und%undulator_type/=typ) call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
  if (.not.ctx%und%model_undulator_ends) call MPI_ABORT(MPI_COMM_WORLD,3,ierr)
  trans=sqrt(ctx%frame%gain_length*ctx%frame%cooperation_length)
  open(21,file=points,status='old')
  n=0;err=0.0_wp
  do
    read(21,*,iostat=ios) z,x(1),y(1),ex,ey,ez
    if (ios/=0) exit
    x=x/trans;y=y/trans
    call getBFields(x,y,z/ctx%frame%gain_length,bx,by,bz,ctx%und,ctx%frame)
    err=max(err,abs(bx(1)-ex),abs(by(1)-ey),abs(bz(1)-ez))
    n=n+1
  end do
  close(21)
  print*, 'FIELD_ROUTE_TYPE ',trim(ctx%und%undulator_type)
  print*, 'FIELD_ROUTE_MODEL_ENDS ',ctx%und%model_undulator_ends
  print*, 'FIELD_POINTS_ERROR ',n,err
  if (n<10000 .or. err>2e-10_wp) call MPI_ABORT(MPI_COMM_WORLD,4,ierr)
  print*, 'FIELD_ROUTE_PASS'
  call MPI_FINALIZE(ierr)
end program field_probe11
