module mod_usr
  use mod_hd

  implicit none
  double precision :: rhodens,rholight=1.0d0,Atwoods=0.8,Reynolds=2000 
  double precision :: lamda=0.4  !lamda=xprobmax1-xprobmin1 
  double precision :: pint, pbottom

  ! the location of demarcation line  
  double precision :: y0=1.6d0

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr
    usr_gravity       => gravity
    usr_special_bc    => specialbound_usr
    
    unit_time=sqrt(lamda/Atwoods)
    unit_velocity=sqrt(lamda*Atwoods/(1.0+Atwoods))
    unit_numberdensity=1.0d3
  
    call set_coordinate_system("Cartesian")
    call hd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ix^D

    double precision:: epsilon,kx
    logical::          first
    data first/.true./

    pint = (one+Atwoods)/(one-Atwoods)*(-y0)
    pbottom = pint+w(ixOmin1, ixOmin2, rho_)*y0

    ! density of two types
    rhodens=rholight*(one+Atwoods)/(one-Atwoods)

    ! setup the perturbation
    epsilon=0.5d0
    ! kx=2 pi
    kx=8.0d0*atan(one)
    lamda=xprobmax1-xprobmin1
    ! print out the info
    if (first) then
       if (mype==0) then
          print *,'HD Rayleigh Taylor problem'
          print *,'  --assuming y ranging from ', xprobmin2,'-',xprobmax2
          print *,'  --interface y0-epsilon:',y0,epsilon
          print *,'  --density ratio:',rhodens/rholight
          print *,'  --kx:',kx
          print *,'unit_time',unit_time
          print *,'unit_velocity',unit_velocity
          print *,'unit_numberdensity',unit_numberdensity
          print *,'lamda',lamda
          print *,'vc_mu in the .par file should be', lamda*sqrt(Atwoods*1.0*lamda/(Atwoods+1))/(Reynolds*2/(rhodens+rholight))
       end if
       first=.false.
    end if

    w(ixO^S,rho_)=epsilon*(1.d0+ERF(170.d0*(x(ixO^S,2)-y0)-2*cos(kx*x(ixO^S,1))))*(rhodens-rholight)+rholight

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero

    ! pressure at interface

    if(hd_energy) then
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S,e_)=w(ixO^S,e_)/(hd_gamma-one)/pbottom
    end if

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    select case(iB)
    case(3)
      w(ixO^S,rho_)=rholight
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S, e_)=w(ixO^S, e_)/pbottom
      w(ixO^S,mom(2))=0.d0
      w(ixO^S,mom(1))=0.d0
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      w(ixO^S,rho_)=rhodens
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S, e_)=w(ixO^S, e_)/pbottom
      w(ixO^S,mom(2))=0.d0
      w(ixO^S,mom(1))=0.d0
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    end select
  end subroutine specialbound_usr

  ! Calculate gravitational acceleration in each dimension
  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:)=0.d0
    gravity_field(ixO^S,2)=-1.d0

  end subroutine gravity

end module mod_usr
