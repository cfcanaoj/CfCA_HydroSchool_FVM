!=============================================================
! Interface block (do not edit in this tutorial)
!=============================================================
interface

    subroutine GenerateGrid(xf, xv, yf, yv)
    use params, only : nx, ny, ngh, xmax, xmin, ymax, ymin 
    implicit none
    real(8), intent(out) :: xf(:), xv(:)
    real(8), intent(out) :: yf(:), yv(:)
    end subroutine GenerateGrid

    subroutine GenerateProblem(xv, yv, Q, Bs, Bc )
    use params, only : is, ie, js, je, &
                       IDN, IVX, IVY, IVZ, IPR
    implicit none
    real(8), intent(in ) :: xv(:), yv(:)
    real(8), intent(out) :: Q(:,:,:)
    real(8), intent(out) :: Bs(:,:,:)
    real(8), intent(out) :: Bc(:,:,:)
    end subroutine GenerateProblem


    subroutine BoundaryCondition( Q, Bs, Bc )
    use params, only : is, ie, js, je, ngh, &
                       IDN, IVX, IVY, IVZ, IPR
    implicit none
    real(8), intent(inout) :: Q(:,:,:)
    real(8), intent(inout) :: Bs(:,:,:)
    real(8), intent(inout) :: Bc(:,:,:)
    end subroutine BoundaryCondition


    subroutine Prim2Consv(Q, Bc, U)
    use params, only : is, ie, js, je, &
                           IDN, IVX, IVY, IVZ, IPR, &
                           IMX, IMY, IMZ, IEN, gam
    implicit none
    real(8), intent(in) :: Q(:,:,:)
    real(8), intent(in) :: Bc(:,:,:)
    real(8), intent(out) :: U(:,:,:)
    end subroutine Prim2Consv


    subroutine Consv2Prim( U, Bs, Q, Bc )
    use params, only : is, ie, js, je, &
                       IDN, IVX, IVY, IVZ, IPR, &
                       IMX, IMY, IMZ, IEN, gam
    implicit none
    real(8), intent(in) :: U(:,:,:),Bs(:,:,:)
    real(8), intent(out) :: Q(:,:,:),Bc(:,:,:)
    end subroutine Consv2Prim

    real(8) function TimestepControl(xf, yf, Q, Bc )
    use params, only : is, ie, js, je, &
                       IDN, IVX, IVY, IVZ, IPR, gam, Ccfl
    implicit none
    real(8), intent(in) :: xf(:), yf(:), Q(:,:,:), Bc(:,:,:)
    end function TimestepControl

    subroutine NumericalFlux( dt, xf, yf, Q, Bc, Bs, F, G, E)
    use params, only : nxtot, nytot, NFLX, NVAR, is, ie, js, je, &
                       IBX, IBY, IBZ, flag_flux
    implicit none
    real(8), intent(in) :: dt, xf(:), yf(:)
    real(8), intent(in) :: Q(:,:,:)
    real(8), intent(in) :: Bc(:,:,:)
    real(8), intent(in) :: Bs(:,:,:)
    real(8), intent(out) :: F(:,:,:)
    real(8), intent(out) :: G(:,:,:)
    real(8), intent(out) :: E(:,:,:)
    end subroutine Numericalflux

    subroutine UpdateConsv( dt1, xf, xv, yf, yv, F, G, E, Q, Uo, Bso, U, Bs)
    use params, only : is, ie, js, je
    implicit none
    real(8), intent(in) :: dt1
    real(8), intent(in)  :: xv(:), yv(:)
    real(8), intent(in)  :: xf(:), yf(:)
    real(8), intent(in)  :: F(:,:,:), G(:,:,:)
    real(8), intent(in)  :: Uo(:,:,:), Q(:,:,:), Bso(:,:,:)
    real(8), intent(out) :: U(:,:,:), Bs(:,:,:), E(:,:,:)
    end subroutine UpdateConsv

    subroutine Output( time, flag, dirname, xv, yv, Q, Bc )
    use params, only : nx, ny, is, ie, js, je, NFLX, NVAR, unitsnap, unitbin, &
                      IDN, IVX, IVY, IVZ, IPR, dtsnap, flag_binary
    implicit none
    real(8), intent(in) :: time
    logical, intent(in) :: flag
    character(20), intent(in) :: dirname 
    real(8), intent(in) :: xv(:), yv(:)
    real(8), intent(in) :: Q(:,:,:), Bc(:,:,:)
    end subroutine Output

    subroutine RealtimeAnalysis(xv,yv,Q,Bc,Bs,phys_evo)
    use params, only : is, ie, js, je, gam
    implicit none
    real(8), intent(in) :: xv(:), yv(:), Q(:,:,:), Bc(:,:,:), Bs(:,:,:)
    real(8), intent(out) :: phys_evo(:)
    end subroutine RealtimeAnalysis


end interface
