program field_probe
  use puffin_kinds
  use puffin_constants
  use BFields
  use setupcalcs, only: calcScaling
  use GlobalTypes, only: tSimulationContext
  implicit none
  type(tSimulationContext) :: ctx
  character(len=32) :: undtype
  real(kind=wp) :: x(5), y(5), bx(5), by(5), bz(5), z, yy(5)
  real(kind=wp) :: ux, uy, nb, lg, lc
  integer :: i, j
  ux=0.0_wp
  uy=1.0_wp
  undtype='aprl-ramp77'
  call calcScaling(0.005_wp,2.060074987322757_wp,298.26088438022981_wp, &
                   0.028_wp,undtype,ux,uy,ctx)
  lg=ctx%frame%gain_length
  lc=ctx%frame%cooperation_length
  nb=lg*lc**2*e_0*m_e/q_e**2*ctx%frame%gamma_ref**3*ctx%frame%rho**3* &
      (4.0_wp*c*2.0_wp*pi/ctx%frame%lambda_w/ctx%frame%aw)**2
  write(*,'(a,10es25.16)') 'SCALE ',pi,q_e,e_0,m_e,lg,lc,nb,ctx%frame%eta, &
       ctx%frame%kappa,ctx%frame%lambda_r
  yy=[-0.0008_wp,-0.000075_wp,0.0_wp,0.000075_wp,0.0008_wp]
  y=yy/sqrt(lg*lc)
  x=0.0_wp
  ctx%und%undulator_type='aprl-ramp77'
  do i=0,3000
    if (i<=1000) then
      z=0.056_wp*real(i,wp)/1000.0_wp
    else if (i<=2000) then
      z=2.100_wp+0.056_wp*real(i-1000,wp)/1000.0_wp
    else
      z=2.156_wp*real(i-2000,wp)/1000.0_wp
    end if
    call getBFields(x,y,z/lg,bx,by,bz,ctx%und,ctx%frame)
    do j=1,5
      write(*,'(a,5es25.16)') 'FIELD ',z,yy(j),bx(j),by(j),bz(j)
    end do
  end do
end program field_probe
