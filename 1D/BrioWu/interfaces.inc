!=============================================================
! Interface block (do not edit in this tutorial)
!=============================================================
interface

    subroutine GenerateGrid(xf, xv)
    use params, only : nxtot, ngh, nx, xmax, xmin
    implicit none
    real(8), intent(out) :: xf(:), xv(:)
    end subroutine GenerateGrid

    subroutine GenerateProblem(xv, Q)
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, is, ie 
    implicit none
    real(8), intent(in ) :: xv(:)
    real(8), intent(out) :: Q(:,:)
    end subroutine GenerateProblem

    subroutine Prim2Consv(Q, U)
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, is, ie, Bx, gam
    implicit none
    real(8), intent(in) :: Q(:,:)
    real(8), intent(out) :: U(:,:)
    end subroutine Prim2Consv

    subroutine Output( time, flag, xv, Q )
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, is, ie, Bx, unitsnap
    implicit none
    real(8),       intent(in) :: time
    logical,       intent(in) :: flag
    real(8),       intent(in) :: xv(:), Q(:,:)
    end subroutine Output

    Real(8) Function TimestepControl(xf, Q) 
    use params, only : IDN, IVX, IPR, IBY, IBZ, Bx, gam, is, ie
    implicit none
    real(8), intent(in) :: xf(:), Q(:,:)
    end Function TimestepControl

    subroutine BoundaryCondition(Q)
    use params, only : NVAR, ngh, is, ie
    implicit none
    real(8), intent(inout) :: Q(:,:)
    end subroutine BoundaryCondition

    subroutine NumericalFlux( dt, xv, Q, flux )
    use params, only : NVAR, is, ie, nxtot, NFLX
    implicit none
    real(8), intent(in) :: dt
    real(8), intent(in) :: xv(:)
    real(8), intent(in) :: Q(:,:)
    real(8), intent(out) :: flux(:,:)
    end subroutine NumericalFlux

    subroutine UpdateConsv( dt, xf, flux, U )
    use params, only : NVAR, is, ie
    implicit none
    real(8), intent(in)  :: flux(:,:), dt, xf(:)
    real(8), intent(out) :: U(:,:)
    end subroutine UpdateConsv

    subroutine Consv2Prim( U, Q )
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, is, ie, Bx, gam
    implicit none
    real(8), intent(in) :: U(:,:)
    real(8), intent(out) :: Q(:,:)
    end subroutine Consv2Prim

    subroutine RealtimeAnalysis(xv,Q,phys_evo)
    use params, only : is, ie, nx
    implicit none
    real(8), intent(in)  :: xv(:), Q(:,:)
    real(8), intent(out) :: phys_evo(:)
    end subroutine RealtimeAnalysis

end interface

