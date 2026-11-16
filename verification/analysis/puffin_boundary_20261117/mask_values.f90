! Read-only evaluation of the masks already compiled into the production library.
! No MPI initialization, field propagation, particle evolution, or source change.
program mask_values
  use puffin_kinds, only: wp, ip
  use puffin_constants, only: pi
  use masks, only: getZ2Mask, getMask
  implicit none
  integer(ip) :: nz, i, j, u, case_i
  real(wp) :: dz, transverse_dx, a(81*81), z(37), zrank(750)
  character(1024) :: out, arg
  call get_command_argument(1,out)
  call get_command_argument(2,arg)
  read(arg,*) dz
  call get_command_argument(3,arg)
  read(arg,*) transverse_dx
  open(newunit=u,file=trim(out)//'/native_tail_masks.csv',status='new')
  write(u,'(a)') 'nz,node_zero_based,dz2_bar,pi_native,mask_native'
  do case_i=-1,1
    nz=96000+20*case_i
    z=getZ2Mask(dz,nz,37_ip,37_ip,nz-37_ip)
    do i=1,37
      write(u,'(i0,a,i0,3(a,es26.17e3))') nz,',',nz-37+i-1,',',dz,',',pi,',',z(i)
    end do
  end do
  close(u)
  zrank=getZ2Mask(dz,96000_ip,750_ip,37_ip,95250_ip)
  open(newunit=u,file=trim(out)//'/native_tail_actual_rank.csv',status='new')
  write(u,'(a)') 'node_zero_based,mask_native'
  do i=714,750
    write(u,'(i0,a,es26.17e3)') 95250+i-1,',',zrank(i)
  end do
  close(u)
  call getMask(81_ip,81_ip,transverse_dx,transverse_dx,16_ip,16_ip,a)
  open(newunit=u,file=trim(out)//'/native_transverse_mask.csv',status='new')
  write(u,'(a)') 'ix_zero_based,iy_zero_based,mask_native'
  do j=1,81
    do i=1,81
      write(u,'(i0,a,i0,a,es26.17e3)') i-1,',',j-1,',',a(i+81*(j-1))
    end do
  end do
  close(u)
end program mask_values
