program route_probe
  use mpi
  use puffin_kinds
  use puffin_constants, only: pi
  use puffin_mpiInfo, only: tProcInfo_G
  use globals, only: qUndEnds_G, qscaled_G, zUndType_G, fx_G, fy_G
  use lattice, only: setupMods, initUndulator
  use setupcalcs, only: calcScaling
  use AdapterGlobals, only: PopulateUndulatorFromGlobals
  use GlobalTypes, only: tSimulationContext
  use bfields, only: getBFields
  implicit none
  type(tSimulationContext) :: ctx
  character(len=1024) :: filename
  character(len=32) :: undtype
  real(kind=wp) :: ux,uy,taper,dz,kbx,kby,sz,szl,x(1),y(1),bx(1),by(1),bz(1)
  integer :: nsteps,ierr
  call MPI_INIT(ierr)
  tProcInfo_G%comm=MPI_COMM_WORLD
  call MPI_COMM_RANK(MPI_COMM_WORLD,tProcInfo_G%rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,tProcInfo_G%size,ierr)
  tProcInfo_G%qroot=tProcInfo_G%rank==0
  call get_command_argument(1,filename)
  undtype='aprl-ramp77'
  qscaled_G=.true.
  qUndEnds_G=.true.
  zUndType_G=undtype
  ux=0.0_wp;uy=1.0_wp;taper=0.0_wp;kbx=0.0_wp;kby=0.0_wp
  call calcScaling(.005_wp,2.060074987322757_wp,298.26088438022981_wp,.028_wp,undtype,ux,uy,ctx)
  nsteps=560;dz=4*pi*.005_wp/80
  call setupMods(filename,taper,.005_wp,nsteps,dz,ux,uy,kbx,kby,ctx%frame)
  sz=0.0_wp
  call initUndulator(1,sz,szl,ctx%frame,ctx%und)
  call PopulateUndulatorFromGlobals(ctx%und)
  if (ctx%und%undulator_type/='aprl-ramp77') call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
  if (.not.ctx%und%model_undulator_ends) call MPI_ABORT(MPI_COMM_WORLD,3,ierr)
  x=0.0_wp;y=0.0_wp
  call getBFields(x,y,.007_wp/ctx%frame%gain_length,bx,by,bz,ctx%und,ctx%frame)
  if (abs(by(1)-.5_wp*(1-sqrt(.5_wp)))>1e-12_wp) call MPI_ABORT(MPI_COMM_WORLD,4,ierr)
  if (tProcInfo_G%qroot) then
    print*, 'ROUTE_TYPE ',trim(ctx%und%undulator_type)
    print*, 'ROUTE_MODEL_ENDS ',ctx%und%model_undulator_ends
    print*, 'ROUTE_FIELD_AT_QUARTER_PERIOD ',by(1)
    print*, 'ROUTE_PASS'
  end if
  call MPI_FINALIZE(ierr)
end program route_probe
