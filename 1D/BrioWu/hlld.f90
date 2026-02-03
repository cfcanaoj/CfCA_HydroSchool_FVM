!=============================================================
! HLLD
! Description:
!   Compute the interface flux using the HLLD approximate Riemann solver for
!   the 1D Euler equations.
!
! Inputs:
!   Ql(:)  Left primitive state  (rho, v, p)
!   Qr(:)  Right primitive state (rho, v, p)
!
! Output:
!   flx(:) Conservative flux (mass, momentum, energy)
! 
! References
!   Miyoshi and Kusano, J. Comput. Phys., 208, 315 (2005)
!=============================================================
subroutine HLLD(Ql,Qr,flx)
implicit none
real(8),intent(in)  ::Ql(:), Qr(:)
real(8),intent(out) :: flx(:)
real(8):: b1
real(8):: Ul(NFLX), Ur(NFLX)
real(8):: Ulst(NFLX), Urst(NFLX)
real(8):: Uldst(NFLX), Urdst(NFLX)
real(8):: Fl(NFLX), Fr(NFLX)
real(8):: cfl,cfr
real(8):: S0, S1, S2, S3, S4
real(8):: pbl, pbr, ptotl, ptotr
real(8) :: sqrtdl, sqrtdr, v_dot_B_stl, v_dot_B_str
real(8) :: Ulst_d_inv, Urst_d_inv, sum_sqrtd_inv, tmp
real(8) :: ptot_stl, ptot_str,ptot_st, Cl, Cr, Cml, Cmr, Cml_inv, Cmr_inv, bxsgn

    pbl = 0.5d0*(Bx**2 + Ql(IBY)**2 + Ql(IBZ)**2)
    pbr = 0.5d0*(Bx**2 + Qr(IBY)**2 + Qr(IBZ)**2)
    ptotl = Ql(IPR) + pbl
    ptotr = Qr(IPR) + pbr

    cfl = sqrt( 0.5d0*( 2.0d0*pbl + gam*Ql(IPR) &
                     + dsqrt( (2.0d0*pbl - gam*Ql(IPR))**2 &
                     + 4.0d0*gam*Ql(IPR)*( Ql(IBY)**2 + Ql(IBZ)**2 ) ) )/Ql(IDN) )
    cfr = sqrt( 0.5d0*( 2.0d0*pbr + gam*Qr(IPR) &
                     + dsqrt( (2.0d0*pbr - gam*Qr(IPR))**2 &
                     + 4.0d0*gam*Qr(IPR)*( Qr(IBY)**2 + Qr(IBZ)**2 ) ) )/Qr(IDN) )

    S0 = min( Ql(IVX) - cfl, Qr(IVX) - cfr)
    S4 = max( Ql(IVX) + cfl, Qr(IVX) + cfr)

    ! conserved variables in the left and right states
    Ul(IDN) = Ql(IDN)
    Ul(IVX) = Ql(IDN)*Ql(IVX)
    Ul(IVY) = Ql(IDN)*Ql(IVY)
    Ul(IVZ) = Ql(IDN)*Ql(IVZ)
    Ul(IEN) = 0.5d0*Ql(IDN)*( Ql(IVX)**2 + Ql(IVY)**2 + Ql(IVZ)**2) & 
            + pbl + Ql(IPR)/(gam - 1.0d0)
    Ul(IBY) = Ql(IBY)
    Ul(IBZ) = Ql(IBZ)

    Ur(IDN) = Qr(IDN)
    Ur(IVX) = Qr(IDN)*Qr(IVX)
    Ur(IVY) = Qr(IDN)*Qr(IVY)
    Ur(IVZ) = Qr(IDN)*Qr(IVZ)
    Ur(IEN) = 0.5d0*Qr(IDN)*( Qr(IVX)**2 + Qr(IVY)**2 + Qr(IVZ)**2) & 
            + pbr + Qr(IPR)/(gam - 1.0d0)
    Ur(IBY) = Qr(IBY)
    Ur(IBZ) = Qr(IBZ)

    !--- Step 3.  Compute L/R fluxes
    Fl(IDN) = Ul(IVX)
    Fl(IVX) = Ul(IVX)*Ql(IVX) + ptotl - Bx**2
    Fl(IVY) = Ul(IVY)*Ql(IVX) - Bx*Ql(IBY)
    Fl(IVZ) = Ul(IVZ)*Ql(IVX) - Bx*Ql(IBZ)
    Fl(IEN) = ( Ul(IEN) + ptotl - Bx**2 )*Ql(IVX) &
            - Bx*( Ql(IBY)*Ql(IVY) + Ql(IBZ)*Ql(IVZ) )
    Fl(IBY) = Ql(IBY)*Ql(IVX) - Bx*Ql(IVY)
    Fl(IBZ) = Ql(IBZ)*Ql(IVX) - Bx*Ql(IVZ)

    Fr(IDN) = Ur(IVX)
    Fr(IVX) = Ur(IVX)*Qr(IVX) + ptotr - Bx**2
    Fr(IVY) = Ur(IVY)*Qr(IVX) - Bx*Qr(IBY)
    Fr(IVZ) = Ur(IVZ)*Qr(IVX) - Bx*Qr(IBZ)
    Fr(IEN) = ( Ur(IEN) + ptotr - Bx**2 )*Qr(IVX) &
            - Bx*( Qr(IBY)*Qr(IVY) + Qr(IBZ)*Qr(IVZ) )
    Fr(IBY) = Qr(IBY)*Qr(IVX) - Bx*Qr(IVY)
    Fr(IBZ) = Qr(IBZ)*Qr(IVX) - Bx*Qr(IVZ)

    !--- Step 4.  Compute middle and Alfven wave speeds
    Cl = S0 - Ql(IVX)
    Cr = S4 - Qr(IVX)

    S2 = ( Cr*Ur(IVX) - Cl*Ul(IVX) + (ptotl - ptotr) ) &
           /( Cr*Ur(IDN) - Cl*Ul(IDN) )

    Cml = S0 - S2
    Cmr = S4 - S2
    Cml_inv = 1.0d0/Cml
    Cmr_inv = 1.0d0/Cmr

    Ulst(IDN) = Ul(IDN)*Cl*Cml_inv
    Urst(IDN) = Ur(IDN)*Cr*Cmr_inv
    Ulst_d_inv = 1.0d0/Ulst(IDN)
    Urst_d_inv = 1.0d0/Urst(IDN)
    sqrtdl = dsqrt(Ulst(IDN))
    sqrtdr = dsqrt(Urst(IDN))

    S1 = S2 - dabs(Bx)/sqrtdl
    S3 = S2 + dabs(Bx)/sqrtdr

    !--- Step 5.  Compute intermediate states
   ptot_stl = ptotl + Ul(IDN)*Cl*(S2 - Ql(IVX))
   ptot_str = ptotr + Ur(IDN)*Cr*(S2 - Qr(IVX))

   ptot_st = 0.5d0*(ptot_stl + ptot_str)

   Ulst(IVX) = Ulst(IDN)*S2
   if( dabs( Ul(IDN)*Cl*Cml-Bx**2) < 1.0d-8*ptot_st ) then
       Ulst(IVY) = Ulst(IDN)*Ql(IVY)
       Ulst(IVZ) = Ulst(IDN)*Ql(IVZ)

       Ulst(IBY) = Ul(IBY)
       Ulst(IBZ) = Ul(IBZ)
   else 
       tmp = Bx*( Cl - Cml )/(Ul(IDN)*Cl*Cml - Bx**2)
       Ulst(IVY) = Ulst(IDN)*( Ql(IVY) - Ul(IBY)*tmp )
       Ulst(IVZ) = Ulst(IDN)*( Ql(IVZ) - Ul(IBZ)*tmp )

       tmp = (Ul(IDN)*Cl**2 - Bx**2)/( Ul(IDN)*Cl*Cml - Bx**2)
       Ulst(IBY) = Ul(IBY)*tmp
       Ulst(IBZ) = Ul(IBZ)*tmp
   endif

   v_dot_B_stl = ( Ulst(IVX)*Bx + Ulst(IVY)*Ulst(IBY) + Ulst(IVZ)*Ulst(IBZ) )*Ulst_d_inv
   Ulst(IEN) = ( Cl*Ul(IEN) - ptotl*Ql(IVX) + ptot_st*S2 &
               + Bx*( Ql(IVX)*Bx + Ql(IVY)*Ul(IBY) + Ql(IVZ)*Ul(IBZ) - v_dot_B_stl) )*Cml_inv

   Urst(IVX) = Urst(IDN)*S2
   if( dabs( Ur(IDN)*Cr*Cmr-Bx**2) < 1.0d-8*ptot_st ) then
       Urst(IVY) = Urst(IDN)*Qr(IVY)
       Urst(IVZ) = Urst(IDN)*Qr(IVZ)

       Urst(IBY) = Ur(IBY)
       Urst(IBZ) = Ur(IBZ)
   else 
       tmp = Bx*( Cr - Cmr )/(Ur(IDN)*Cr*Cmr - Bx**2)
       Urst(IVY) = Urst(IDN)*( Qr(IVY) - Ur(IBY)*tmp )
       Urst(IVZ) = Urst(IDN)*( Qr(IVZ) - Ur(IBZ)*tmp )

       tmp = (Ur(IDN)*Cr**2 - Bx**2)/( Ur(IDN)*Cr*Cmr - Bx**2)
       Urst(IBY) = Ur(IBY)*tmp
       Urst(IBZ) = Ur(IBZ)*tmp
   endif

   v_dot_B_str = ( Urst(IVX)*Bx + Urst(IVY)*Urst(IBY) + Urst(IVZ)*Urst(IBZ) )*Urst_d_inv
   Urst(IEN) = ( Cr*Ur(IEN) - ptotr*Qr(IVX) + ptot_st*S2 &
               + Bx*( Qr(IVX)*Bx + Qr(IVY)*Ur(IBY) + Qr(IVZ)*Ur(IBZ) - v_dot_B_str) )*Cmr_inv
 
   if( 0.5d0*Bx**2 < 1.0d-8*ptot_st )then
       Uldst(:) = Ulst(:)
       Urdst(:) = Urst(:)
   else 
       sum_sqrtd_inv = 1.0d0/(sqrtdl + sqrtdr) 
       if (Bx > 0.0d0 ) then 
           bxsgn = 1.0d0
       else 
           bxsgn = -1.0d0
       endif

       Uldst(IDN) = Ulst(IDN)
       Urdst(IDN) = Urst(IDN)

       Uldst(IVX) = Ulst(IVX)
       Urdst(IVX) = Urst(IVX)

       tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVY)*Ulst_d_inv) + sqrtdr*(Urst(IVY)*Urst_d_inv) &
                            + bxsgn*(Urst(IBY) - Ulst(IBY)) )
       Uldst(IVY) = Uldst(IDN)*tmp
       Urdst(IVY) = Urdst(IDN)*tmp
!

       tmp = sum_sqrtd_inv*(  sqrtdl*(Ulst(IVZ)*Ulst_d_inv) + sqrtdr*(Urst(IVZ)*Urst_d_inv) &
                            + bxsgn*(Urst(IBZ) - Ulst(IBZ)) )
       Uldst(IVZ) = Uldst(IDN)*tmp
       Urdst(IVZ) = Urdst(IDN)*tmp

       tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBY) + sqrtdr*Ulst(IBY) &
                 + bxsgn*sqrtdl*sqrtdr*( (Urst(IVY)*Urst_d_inv) - (Ulst(IVY)*Ulst_d_inv) ) )
       Uldst(IBY) = tmp
       Urdst(IBY) = tmp

       tmp = sum_sqrtd_inv*(  sqrtdl*Urst(IBZ) + sqrtdr*Ulst(IBZ) &
                 + bxsgn*sqrtdl*sqrtdr*( (Urst(IVZ)*Urst_d_inv) - (Ulst(IVZ)*Ulst_d_inv) ) )
       Uldst(IBZ) = tmp
       Urdst(IBZ) = tmp
!
       tmp = S2*Bx + (Uldst(IVY)*Uldst(IBY) + Uldst(IVZ)*Uldst(IBZ))/Uldst(IDN)
       Uldst(IEN) = Ulst(IEN) - sqrtdl*bxsgn*(v_dot_B_stl - tmp)
       Urdst(IEN) = Urst(IEN) + sqrtdr*bxsgn*(v_dot_B_str - tmp)
   endif

    !--- Step 6.  Compute flux
    if( S0 >= 0.0d0 ) then
         flx(:) = Fl(:)
    else if (S4 <= 0.0d0 ) then
         flx(:) = Fr(:)
    else  if  (S1 >= 0.0d0) then
         flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:))
     else if (S2 >= 0.0d0) then
         flx(:) = Fl(:) + S0*(Ulst(:) - Ul(:)) + S1*(Uldst(:) - Ulst(:))
     else if (S3 > 0.0d0 ) then
         flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) + S3*(Urdst(:) - Urst(:))
     else 
         flx(:) = Fr(:) + S4*(Urst(:) - Ur(:)) 
     endif

return
end subroutine HLLD
