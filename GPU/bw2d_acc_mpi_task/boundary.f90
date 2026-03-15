!=============================================================
! BoundaryCondition
! Description:
!   Apply boundary conditions by filling ghost cells
!   of the primitive array Q.
!=============================================================
subroutine BoundaryCondition(Q)
  use params, only : nxtot, nytot, NVAR, ngh, is, ie, js, je
  use mpimod
  implicit none
  real(8), intent(inout) :: Q(NVAR,nxtot,nytot)
  integer::i,j,ihy
  integer::sizex,sizey
  real(8),dimension(NVAR,ngh,nytot-1):: QsendXm,QsendXp
  real(8),dimension(NVAR,ngh,nytot-1):: QrecvXm,QrecvXp
  real(8),dimension(NVAR,nxtot-1,ngh):: QsendYm,QsendYp
  real(8),dimension(NVAR,nxtot-1,ngh):: QrecvYm,QrecvYp
    
  if(ntiles(1) ==1) then
     !$acc parallel loop collapse(3) present(Q)
     do j=1,nytot-1
        do i=1,ngh
           do ihy=1,NVAR
              Q(ihy,is-i,j) = Q(ihy,ie+1-i,j)
              Q(ihy,ie+i,j) = Q(ihy,is+i-1,j)
           enddo
        enddo
     enddo
  endif
  
  if(ntiles(2) ==1) then
     !$acc parallel loop collapse(3) present(Q)
     do j=1,ngh
        do i=1,nxtot-1
           do ihy=1,NVAR
              Q(ihy,i,js-j)  = Q(ihy,i,je+1-j)
              Q(ihy,i,je+j)  = Q(ihy,i,js+j-1)
           enddo
        enddo
     enddo
  endif

  !$acc data create(QsendXp,QsendXm,QrecvXp,QrecvXm)
  mpixdir: if(ntiles(1) /=1) then
     !$acc parallel loop collapse(3) present(Q)
     do j=1,nytot-1
        do i=1,ngh
           do ihy=1,NVAR
              QsendXp(ihy,i,j) = Q(ihy,ie-ngh+i,j)
              QsendXm(ihy,i,j) = Q(ihy,is  +i-1,j)
           enddo
        enddo
     enddo
     
     n1mdir: if (n1m /= MPI_PROC_NULL) then
        sizex = NVAR*ngh*(nytot-1)
        nreq = nreq + 1
        call MPI_IRECV(QrecvXm,sizex &
             & , MPI_DOUBLE &
             & , n1m,1100, comm3d, req(nreq), ierr)
     
        nreq = nreq + 1
        call MPI_ISEND(QsendXm,sizex &
             & , MPI_DOUBLE &
             & , n1m, 1200, comm3d, req(nreq), ierr)
     endif n1mdir

     n1pdir: if (n1p /= MPI_PROC_NULL) then
        sizex = NVAR*ngh*(nytot-1)
        nreq = nreq + 1
        call MPI_IRECV(QrecvXp,sizex &
             & , MPI_DOUBLE &
             & , n1p,1200, comm3d, req(nreq), ierr)
     
        nreq = nreq + 1
        call MPI_ISEND(QsendXp,sizex &
             & , MPI_DOUBLE &
             & , n1p, 1100, comm3d, req(nreq), ierr)
     endif n1pdir
     
     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0

     !$acc parallel loop collapse(3) present(Q)
     do j=1,nytot-1
        do i=1,ngh
           do ihy=1,NVAR
              Q(ihy,is-ngh-1+i,j) = QrecvXm(ihy,i,j)
              Q(ihy,ie+i      ,j) = QrecvXp(ihy,i,j)
           enddo
        enddo
     enddo
  endif mpixdir
  !$acc end data
  
  !$acc data create(QsendYp,QsendYm,QrecvYp,QrecvYm)
  mpiydir: if(ntiles(2) /=1) then
     !$acc parallel loop collapse(3) present(Q)
     do i=1,nxtot-1
        do j=1,ngh
           do ihy=1,NVAR
              QsendYp(ihy,i,j) = Q(ihy,i,je-ngh+j)
              QsendYm(ihy,i,j) = Q(ihy,i,js+j-1)
           enddo
        enddo
     enddo
     
     n2mdir: if (n2m /= MPI_PROC_NULL) then
        sizey = NVAR*ngh*(nxtot-1)
        nreq = nreq + 1
        call MPI_IRECV(QrecvXm,sizey &
             & , MPI_DOUBLE &
             & , n2m,2100, comm3d, req(nreq), ierr)
     
        nreq = nreq + 1
        call MPI_ISEND(QsendXm,sizey &
             & , MPI_DOUBLE &
             & , n2m, 2200, comm3d, req(nreq), ierr)
     endif n2mdir

     n2pdir: if (n2p /= MPI_PROC_NULL) then
        sizey = NVAR*ngh*(nxtot-1)
        nreq = nreq + 1
        call MPI_IRECV(QrecvXp,sizey &
             & , MPI_DOUBLE &
             & , n2p,2200, comm3d, req(nreq), ierr)
     
        nreq = nreq + 1
        call MPI_ISEND(QsendXp,sizey &
             & , MPI_DOUBLE &
             & , n2p, 2100, comm3d, req(nreq), ierr)
     endif n2pdir
     
     if(nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
     nreq = 0

     !$acc parallel loop collapse(3) present(Q)
     do i=1,nxtot-1
        do j=1,ngh
           do ihy=1,NVAR
              Q(ihy,i,js-ngh-1+j) = QrecvYm(ihy,i,j)
              Q(ihy,i,je+j      ) = QrecvYp(ihy,i,j)
           enddo
        enddo
     enddo
  endif mpiydir
  !$acc end data
  
  
  
return
end subroutine BoundaryCondition
