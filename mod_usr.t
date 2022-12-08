module mod_usr
  use mod_hd

  implicit none
  double precision :: rhodens,rholight=1.0d0,Atwoods=0.04d0,Reynolds=100.0d0 
  double precision :: lamda=0.4  
  double precision :: pint
 
  ! the location of demarcation line  
  double precision :: y0=2.2d0

contains

  subroutine usr_init()
    double precision :: g = 9.81

    usr_init_one_grid => initonegrid_usr
    usr_gravity       => gravity
    usr_special_bc    => specialbound_usr
    
    unit_time=sqrt(g*lamda/Atwoods)
    unit_velocity=sqrt(g*lamda*Atwoods/(1.0+Atwoods))
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

    ! density of two types
    rhodens=rholight*(one+Atwoods)/(one-Atwoods)

    ! setup the perturbation
    epsilon=0.5d0
    lamda=xprobmax1-xprobmin1
    ! kx=2 pi
    kx=8.0d0*atan(one)/lamda
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
          print *,'lamda: ',lamda
          print *,'Atwoods: ', Atwoods
          print *,'vc_mu from .par file', lamda*sqrt(Atwoods*1.0*lamda/(Atwoods+1.0d0))/(2.0*Reynolds/(rhodens + rholight))
       end if
       first=.false.
    end if

    w(ixO^S,rho_)=epsilon*(1.d0+ERF(170.d0*(x(ixO^S,2)-y0)-epsilon*cos(kx*x(ixO^S,1))))*(rhodens-rholight)+rholight

    pint = rholight*(one+Atwoods)/(one-Atwoods)*(xprobmax2-y0)

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero

    ! pressure at interface

    if(hd_energy) then
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S,e_)=w(ixO^S,e_)/(hd_gamma-one)
    end if

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer ::  ix^D,ixOs^L
    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    select case(iB)
    case(3)
      w(ixO^S,rho_)=rholight
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S,mom(2))=0.d0
      w(ixO^S,mom(1))=0.d0
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      w(ixO^S,mom(1))=0
      w(ixO^S,mom(2))=0  
      ixOs^L=ixO^L;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call hd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)-invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      enddo
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    end select
  end subroutine specialbound_usr


  subroutine specialbound_usrOld(qt,ixI^L,ixO^L,iB,w,x)
    use mod_global_parameters
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    select case(iB)
    case(3)
      w(ixO^S,rho_)=rholight
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S, e_)=w(ixO^S, e_)
      w(ixO^S,mom(2))=0.d0
      w(ixO^S,mom(1))=0.d0
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      w(ixO^S,rho_)=rhodens
      w(ixO^S,e_)=pint-w(ixO^S,rho_)*(x(ixO^S,2)-y0)
      w(ixO^S, e_)=w(ixO^S, e_)
      w(ixO^S,mom(2))=0.d0
      w(ixO^S,mom(1))=0.d0
      call hd_to_conserved(ixI^L,ixO^L,w,x)
    end select
  end subroutine specialbound_usrOld

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
