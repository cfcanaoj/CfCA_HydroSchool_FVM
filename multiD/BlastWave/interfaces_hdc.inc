!=============================================================
! Interface block (do not edit in this tutorial)
!=============================================================
interface

    subroutine GenerateGrid(xf, xv, yf, yv)
    use params, only : nxtot, nytot, ngh, nx, ny, xmax, xmin, ymax, ymin
    implicit none
    real(8), intent(out) :: xf(:), xv(:)
    real(8), intent(out) :: yf(:), yv(:)
    end subroutine GenerateGrid

    subroutine GenerateProblem(xv, yv, Q )
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam
    implicit none
    real(8), intent(in ) :: xv(:), yv(:)
    real(8), intent(out) :: Q(:,:,:)
    end subroutine GenerateProblem

    subroutine Prim2Consv(Q, U)
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam
    implicit none
    real(8), intent(in) :: Q(:,:,:)
    real(8), intent(out) :: U(:,:,:)
    end subroutine Prim2Consv

    subroutine BoundaryCondition(Q)
    use params, only : nxtot, nytot, ngh, is, ie, js, je
    implicit none
    real(8), intent(inout) :: Q(:,:,:)
    end subroutine BoundaryCondition

    subroutine Output( time, flag, xv, yv, Q )
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, IPS, NVAR, & 
                       nx, ny, is, ie, js, je, gam, &
                       flag_binary, dirname, dtsnap, unitsnap
    implicit none
    real(8), intent(in) :: time! false --> output per dtsnap, true --> force to output
    logical, intent(in) :: flag
    real(8), intent(in) :: xv(:), yv(:), Q(:,:,:)
    end subroutine Output

    real(8) function TimestepControl(xf, yf, Q)
    use params, only : IDN, IVX, IVY, IPR, IBX, IBY, IBZ, is, ie, js, je, Ccfl, gam
    implicit none
    real(8), intent(in) :: xf(:), yf(:), Q(:,:,:)
    end function TimestepControl

    subroutine NumericalFlux( dt, xf, yf, Q, F, G )
    use params, only : nxtot, nytot, NVAR, NFLX, is, ie, js, je, Ccfl, flag_flux
    implicit none
    real(8), intent(in) :: dt
    real(8), intent(in) :: xf(:), yf(:)
    real(8), intent(in) :: Q(:,:,:)
    real(8), intent(out) :: F(:,:,:)
    real(8), intent(out) :: G(:,:,:)
    end subroutine NumericalFlux

    subroutine UpdateConsv( dt1, xf, yf, F, G, Uo, U)
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam, alpha, Ccfl
    implicit none
    real(8), intent(in) :: dt1
    real(8), intent(in)  :: xf(:), yf(:)
    real(8), intent(in)  :: F(:,:,:), G(:,:,:)
    real(8), intent(in)  :: Uo(:,:,:)
    real(8), intent(inout) :: U(:,:,:)
    end subroutine UpdateConsv

    subroutine SrcTerms( dt1, dt0, Q, U )
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam, alpha, Ccfl
    implicit none
    real(8), intent(in) :: dt1, dt0
    real(8), intent(in)  :: Q(:,:,:)
    real(8), intent(inout) :: U(:,:,:)
    end subroutine SrcTerms

    subroutine Consv2Prim( U, Q )
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam
    implicit none
    real(8), intent(in) :: U(:,:,:)
    real(8), intent(out) :: Q(:,:,:)
    end subroutine Consv2Prim

    subroutine RealtimeAnalysis(xv,yv,Q,phys_evo)
    use params, only : IDN, IVX, IVY, IVZ, IPR, IBX, IBY, IBZ, &
                       IMX, IMY, IMZ, IEN, IPS, NVAR, is, ie, js, je, gam, nevo
    implicit none
    real(8), intent(in) :: xv(:), yv(:), Q(:,:,:)
    real(8), intent(out) :: phys_evo(:)
    end subroutine RealtimeAnalysis

end interface
