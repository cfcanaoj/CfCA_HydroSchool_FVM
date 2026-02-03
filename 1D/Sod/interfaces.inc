!=============================================================
! Interface block (do not edit in this tutorial)
!=============================================================
interface

  subroutine GenerateGrid(xf, xv)
  use params, only : x1min, x1max, nxtotf, nxtotv, ngh, nx
  implicit none 
  real(8), intent(out) :: xf(:), xv(:) 
  end subroutine GenerateGrid

  subroutine GenerateProblem(xv, Q)
  implicit none
  real(8), intent(in ) :: xv(:)
  real(8), intent(out) :: Q(:,:)
  end subroutine GenerateProblem

  subroutine Prim2Consv(Q, U)
  use params, only : IDN, IMX, IEN, IVX, IPR, is, ie, gam
  implicit none
  real(8), intent(in) :: Q(:,:)
  real(8), intent(out) :: U(:,:)
  end subroutine Prim2Consv

  subroutine Output( time, flag, xv, Q )
  use params, only: IDN, IVX, IPR, dirname, is, ie, unitsnap
  implicit none
  real(8), intent(in) :: time
  logical, intent(in) :: flag
  real(8), intent(in) :: xv(:), Q(:,:)
  end subroutine Output

  Real(8) Function TimestepControl(xf, Q) 
  use params, only : IDN, IVX, IPR, is, ie, gam 
  real(8), intent(in) :: xf(:), Q(:,:)
  end Function TimestepControl

  subroutine BoundaryCondition(Q)
  use params, only : IDN, IVX, IPR, is, ie, ngh
  implicit none
  real(8), intent(inout) :: Q(:,:)
  end subroutine BoundaryCondition

  subroutine Consv2Prim( U, Q )
  use params, only : IDN, IVX, IPR, IMX, IEN, gam, is, ie
  implicit none
  real(8), intent(in) :: U(:,:)
  real(8), intent(out) :: Q(:,:)
  end subroutine Consv2Prim

  subroutine NumericalFlux( dt, xv, Q, F )
  use params, only : nxtotf, nvar, is, ie
  implicit none
  real(8), intent(in) :: dt
  real(8), intent(in) :: xv(:)
  real(8), intent(in) :: Q(:,:)
  real(8), intent(out) :: F(:,:)
  end subroutine NumericalFlux

  subroutine UpdateConsv( dt, xf, F, U )
  use params, only : NVAR, is, ie
  implicit none
  real(8), intent(in)  :: F(:,:), dt, xf(:)
  real(8), intent(inout) :: U(:,:)
  end subroutine UpdateConsv

  subroutine RealtimeAnalysis(xv,Q,phys_evo)
  use params, only : IDN, IVX, IPR, NVAR, is, ie, gam, nevo, nx, nxtotv
  implicit none
  real(8), intent(in)  :: xv(:), Q(:,:)
  real(8), intent(out) :: phys_evo(:)
  end subroutine

end interface

