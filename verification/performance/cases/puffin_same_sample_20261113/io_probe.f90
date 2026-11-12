program io_probe
  use mpi
  use hdf5
  implicit none
  integer :: ierr, rank, np, nlocal, offset, j
  integer, parameter :: ntotal=1954743
  integer(hid_t) :: file_id, dset_id, file_space, mem_space, access, transfer
  integer(hsize_t) :: global_shape(2), local_shape(2), start(2)
  real(8), allocatable :: column(:)
  real(8) :: began
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
  offset=ntotal*rank/np
  nlocal=ntotal*(rank+1)/np-offset
  allocate(column(nlocal))
  global_shape=[7,ntotal]
  local_shape=[1,nlocal]
  call h5open_f(ierr)
  call h5pcreate_f(H5P_FILE_ACCESS_F,access,ierr)
  call h5pset_fapl_mpio_f(access,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)
  call h5fcreate_f('probe.h5',H5F_ACC_TRUNC_F,file_id,ierr,access_prp=access)
  call h5screate_simple_f(2,global_shape,file_space,ierr)
  call h5dcreate_f(file_id,'electrons',H5T_NATIVE_DOUBLE,file_space,dset_id,ierr)
  call h5screate_simple_f(2,local_shape,mem_space,ierr)
  call h5pcreate_f(H5P_DATASET_XFER_F,transfer,ierr)
  call h5pset_dxpl_mpio_f(transfer,H5FD_MPIO_COLLECTIVE_F,ierr)
  began=MPI_Wtime()
  do j=0,6
    start=[j,offset]
    column=real(j,8)
    if (rank==0) then
      write(*,*) 'IO_BEGIN_COLUMN',j,MPI_Wtime()-began
      flush(6)
    end if
    call h5sselect_hyperslab_f(file_space,H5S_SELECT_SET_F,start,local_shape,ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,column,local_shape,ierr, &
                   xfer_prp=transfer,file_space_id=file_space,mem_space_id=mem_space)
    if (ierr/=0) call MPI_ABORT(MPI_COMM_WORLD,2,ierr)
  end do
  call h5dclose_f(dset_id,ierr)
  call h5fclose_f(file_id,ierr)
  call h5close_f(ierr)
  if (rank==0) then
    write(*,*) 'IO_PROBE_PASS',MPI_Wtime()-began
    flush(6)
  end if
  call MPI_FINALIZE(ierr)
end program io_probe
