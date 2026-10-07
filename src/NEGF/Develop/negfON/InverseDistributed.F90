!
! Distributed block-tridiagonal inversion over the inverse_comm process group.
!
! Contract (shared with InvertSparseONv3): gfsparse holds E*S-H-Sigma in CRS form on the
! group master; on return its diagonal blocks (and the first off-diagonal blocks for
! opindex 1,3,5) contain the corresponding blocks of G on the stored pattern;
! gfout(:,1:nl) = G(:,1:nl) and gfout(:,nl+1:nl+nr) = G(:,N1-nr+1:N1) for opindex 2,3;
! gfout(1:nr,1:nl) = G(N1-nr+1:N1,1:nl) for opindex 4,5. gfsparse and gfout are
! referenced on mynode_inverse == 0 only; every rank of inverse_comm must call.
!
module mInverseDistributed
  use mConstants
  use mTypes
  use mMatrixUtil
  use mPartition
  use mMPI_NEGF
  use mNegfOutput, only: negf_abort, negf_log_unit
  use negfmod, only: outinfo
  implicit none
  private

  public :: InvertSparseONDistributed
  public :: DistributedInversionActive

  integer, parameter :: tag_blocks = 3101
  integer, parameter :: tag_results = 3102
  logical, save :: first_call = .true.

contains

  logical function DistributedInversionActive(solver)
    integer, intent(in) :: solver
    DistributedInversionActive = (solver == 2 .and. nnodes_inverse > 1)
  end function DistributedInversionActive

  subroutine InvertSparseONDistributed(N1,gfsparse,nl,nr,gfout,opindex)
    integer, intent(in) :: N1,nl,nr,opindex
    type(matrixTypeGeneral), intent(inout) :: gfsparse,gfout

    character(len=*), parameter :: sMyName="InvertSparseONDistributed"
    type(ioType) :: io
    type(matrixType), allocatable :: h0(:),h1(:),hm1(:)
    type(matrixType), allocatable :: lh0(:),lh1(:),lhm1(:)
    type(matrixType), allocatable :: rD(:),rU(:),rL(:),rSigL(:),rSigR(:),rG(:),rC1(:),rCK(:),rMM(:),rML(:)
    type(matrixType), allocatable :: sL(:),M1(:),M2(:),g0(:),g1(:),gm1(:),c1loc(:),cKloc(:)
    type(matrixType) :: Saa,Sab,Sba,Sbb,sR,M,tmp,blk
    integer, allocatable :: nb(:),off(:),ca(:),cb(:),kept(:),ka(:),kb(:)
    complex(kdp), allocatable :: buf(:)
    integer :: me,np,nblk,nchunks,a,b,bb,i,p,k,nred,cnt,pos,mpierror
    logical :: need_offdiag,need_col,need_corner,have_chunk
    real(kdp) :: total
#ifdef MPI
    integer :: istatus(MPI_STATUS_SIZE)
#endif

    io%isDebug=.false.
    me=mynode_inverse
    np=nnodes_inverse
    need_offdiag=(opindex==1.or.opindex==3.or.opindex==5)
    need_col=(opindex==2.or.opindex==3)
    need_corner=(opindex==4.or.opindex==5)
    nblk=0

    if(me==0)then
      if(gfsparse%mattype/=2) call negf_abort(sMyName//": the distributed inverter needs the sparse (EM.OrderN) Green function matrix")
      call PartitionMatrix(h0,h1,hm1,nblk,gfsparse,nl,nr,io)
      call FillBlocksFromMatrixSparse(h0,h1,hm1,nblk,gfsparse)
    endif
#ifdef MPI
    call MPI_Bcast(nblk,1,MPI_integer,0,inverse_comm,mpierror)
#endif
    allocate(nb(nblk),off(nblk))
    if(me==0)then
      do i=1,nblk
        nb(i)=h0(i)%iRows
        off(i)=h0(i)%iVert
        if(h0(i)%iCols/=nb(i).or.h0(i)%iHorz/=off(i)) &
          call negf_abort(sMyName//": non-square diagonal block in the block-tridiagonal partition")
      enddo
      do i=1,nblk-1
        if(off(i+1)/=off(i)+nb(i).or.h1(i)%iRows/=nb(i).or.h1(i)%iCols/=nb(i+1).or.h1(i)%iHorz/=off(i+1).or. &
           h1(i)%iVert/=off(i).or.hm1(i)%iRows/=nb(i+1).or.hm1(i)%iCols/=nb(i).or.hm1(i)%iHorz/=off(i).or. &
           hm1(i)%iVert/=off(i+1)) &
          call negf_abort(sMyName//": off-diagonal block layout does not match the diagonal blocks")
      enddo
      if(off(1)/=1.or.off(nblk)+nb(nblk)-1/=N1.or.nb(1)/=nl.or.nb(nblk)/=nr) &
        call negf_abort(sMyName//": the first and last blocks must be the lead blocks (nl, nr) and cover 1..N1")
    endif
#ifdef MPI
    call MPI_Bcast(nb(1),nblk,MPI_integer,0,inverse_comm,mpierror)
    call MPI_Bcast(off(1),nblk,MPI_integer,0,inverse_comm,mpierror)
#endif

    allocate(ca(0:np-1),cb(0:np-1),kept(0:np-1),ka(0:np-1),kb(0:np-1))
    call AssignChunks(nb,nblk,np,nchunks,ca,cb)
    a=ca(me)
    b=cb(me)
    have_chunk=(a<=b)
    bb=min(b,nblk-1)
    kept=0
    nred=0
    do p=0,nchunks-1
      kept(p)=2
      if(ca(p)==cb(p)) kept(p)=1
      ka(p)=nred+1
      kb(p)=nred+kept(p)
      nred=nred+kept(p)
    enddo

    if(first_call.and.me==0)then
      total=0.0_kdp
      do i=1,nblk
        total=total+real(nb(i),kdp)**3
      enddo
      write(negf_log_unit,'(a,i0,a,i0,a,i0,a,i0)') 'InvertSparseONDistributed: N1=',N1,' blocks=',nblk, &
        ' inverse ranks=',np,' chunks=',nchunks
      do p=0,nchunks-1
        write(negf_log_unit,'(a,i0,a,i0,a,i0,a,f7.4)') '  rank ',p,': blocks ',ca(p),'..',cb(p), &
          ' cost share ',ChunkCost(nb,ca(p),cb(p))/total
      enddo
      first_call=.false.
    endif

    ! scatter the chunks; the master keeps the storage of its own blocks
    if(me==0)then
      do p=1,nchunks-1
        cnt=ChunkBlockCount(nb,nblk,ca(p),cb(p))
        allocate(buf(cnt))
        pos=0
        do i=ca(p),cb(p)
          call PackBlock(buf,pos,h0(i))
        enddo
        do i=ca(p),min(cb(p),nblk-1)
          call PackBlock(buf,pos,h1(i))
          call PackBlock(buf,pos,hm1(i))
        enddo
#ifdef MPI
        call MPI_Send(buf(1),cnt,DAT_dcomplex,p,tag_blocks,inverse_comm,mpierror)
#endif
        deallocate(buf)
      enddo
      allocate(lh0(a:b),lh1(a:bb),lhm1(a:bb))
      do i=a,b
        call move_alloc(h0(i)%a,lh0(i)%a)
        call SetBlock(lh0(i),nb(i),nb(i),off(i),off(i))
      enddo
      do i=a,bb
        call move_alloc(h1(i)%a,lh1(i)%a)
        call SetBlock(lh1(i),nb(i),nb(i+1),off(i+1),off(i))
        call move_alloc(hm1(i)%a,lhm1(i)%a)
        call SetBlock(lhm1(i),nb(i+1),nb(i),off(i),off(i+1))
      enddo
      do i=1,nblk
        call FreeBlock(h0(i))
      enddo
      do i=1,nblk-1
        call FreeBlock(h1(i))
        call FreeBlock(hm1(i))
      enddo
      call DestroyArray(h0,sMyName,io)
      call DestroyArray(h1,sMyName,io)
      call DestroyArray(hm1,sMyName,io)
    elseif(have_chunk)then
      cnt=ChunkBlockCount(nb,nblk,a,b)
      allocate(buf(cnt))
#ifdef MPI
      call MPI_Recv(buf(1),cnt,DAT_dcomplex,0,tag_blocks,inverse_comm,istatus,mpierror)
#endif
      allocate(lh0(a:b),lh1(a:bb),lhm1(a:bb))
      pos=0
      do i=a,b
        call UnpackBlock(buf,pos,lh0(i),nb(i),nb(i),off(i),off(i))
      enddo
      do i=a,bb
        call UnpackBlock(buf,pos,lh1(i),nb(i),nb(i+1),off(i+1),off(i))
        call UnpackBlock(buf,pos,lhm1(i),nb(i+1),nb(i),off(i),off(i+1))
      enddo
      deallocate(buf)
    endif

    ! Schur complement of every chunk onto its boundary blocks
    if(have_chunk) call SchurChunk(lh0,lh1,lhm1,a,b,Saa,Sab,Sba,Sbb,io)

    ! reduced block-tridiagonal system, replicated on every rank
    allocate(rD(nred),rU(max(nred-1,0)),rL(max(nred-1,0)))
    do p=0,nchunks-1
      cnt=nb(ca(p))**2
      if(kept(p)==2) cnt=cnt+2*nb(ca(p))*nb(cb(p))+nb(cb(p))**2
      if(cb(p)<nblk) cnt=cnt+2*nb(cb(p))*nb(cb(p)+1)
      allocate(buf(cnt))
      if(p==me)then
        pos=0
        call PackBlock(buf,pos,Saa)
        if(kept(p)==2)then
          call PackBlock(buf,pos,Sab)
          call PackBlock(buf,pos,Sba)
          call PackBlock(buf,pos,Sbb)
        endif
        if(cb(p)<nblk)then
          call PackBlock(buf,pos,lh1(b))
          call PackBlock(buf,pos,lhm1(b))
        endif
      endif
#ifdef MPI
      call MPI_Bcast(buf(1),cnt,DAT_dcomplex,p,inverse_comm,mpierror)
#endif
      pos=0
      k=ka(p)
      call UnpackBlock(buf,pos,rD(k),nb(ca(p)),nb(ca(p)),0,0)
      if(kept(p)==2)then
        call UnpackBlock(buf,pos,rU(k),nb(ca(p)),nb(cb(p)),0,0)
        call UnpackBlock(buf,pos,rL(k),nb(cb(p)),nb(ca(p)),0,0)
        call UnpackBlock(buf,pos,rD(k+1),nb(cb(p)),nb(cb(p)),0,0)
      endif
      if(cb(p)<nblk)then
        call UnpackBlock(buf,pos,rU(kb(p)),nb(cb(p)),nb(cb(p)+1),0,0)
        call UnpackBlock(buf,pos,rL(kb(p)),nb(cb(p)+1),nb(cb(p)),0,0)
      endif
      deallocate(buf)
    enddo

    if(have_chunk)then
      allocate(rSigL(nred),rSigR(nred),rG(nred),rMM(max(nred-1,0)),rML(max(nred-1,0)))
      call ZeroBlock(rSigL(1),rD(1)%iRows)
      do k=1,nred-1
        call InvMinus(M,rD(k),rSigL(k),io)
        call Mul(rMM(k),-kcone,M,rU(k),io)
        call Mul(rSigL(k+1),-kcone,rL(k),rMM(k),io)
        call FreeBlock(M)
      enddo
      call ZeroBlock(rSigR(nred),rD(nred)%iRows)
      do k=nred,2,-1
        call InvMinus(M,rD(k),rSigR(k),io)
        call Mul(rML(k-1),-kcone,M,rL(k-1),io)
        call Mul(rSigR(k-1),-kcone,rU(k-1),rML(k-1),io)
        call FreeBlock(M)
      enddo
      do k=1,nred
        call InvMinus2(rG(k),rD(k),rSigL(k),rSigR(k),io)
      enddo
      if(need_col.or.need_corner)then
        allocate(rC1(nred))
        call CopyBlock(rC1(1),rG(1),io)
        do k=2,nred
          call Mul(rC1(k),kcone,rML(k-1),rC1(k-1),io)
        enddo
      endif
      if(need_col)then
        allocate(rCK(nred))
        call CopyBlock(rCK(nred),rG(nred),io)
        do k=nred-1,1,-1
          call Mul(rCK(k),kcone,rMM(k),rCK(k+1),io)
        enddo
      endif

      ! chunk solve seeded with the exact boundary self-energies
      allocate(sL(a:b),M1(a:b),M2(a:b),g0(a:b))
      call CopyBlock(sL(a),rSigL(ka(me)),io)
      do i=a,b-1
        call InvMinus(M1(i),lh0(i),sL(i),io)
        call Mul3(sL(i+1),kcone,lhm1(i),M1(i),lh1(i),io)
      enddo
      if(need_offdiag.and.b<nblk) call InvMinus(M1(b),lh0(b),sL(b),io)
      call CopyBlock(sR,rSigR(kb(me)),io)
      do i=b,a,-1
        call InvMinus2(g0(i),lh0(i),sL(i),sR,io)
        if(i>a)then
          call InvMinus(M2(i),lh0(i),sR,io)
          call FreeBlock(sR)
          call Mul3(sR,kcone,lh1(i-1),M2(i),lhm1(i-1),io)
        endif
      enddo
      call FreeBlock(sR)
      if(outinfo)then
        write(negf_log_unit,'(a,i0,2es12.3)') 'distributed inversion boundary check, rank ',me, &
          maxval(abs(rG(ka(me))%a-g0(a)%a)),maxval(abs(rG(kb(me))%a-g0(b)%a))
      endif
      if(need_offdiag)then
        allocate(g1(a:bb),gm1(a:bb))
        do i=a,b-1
          call Mul3(g1(i),-kcone,M1(i),lh1(i),g0(i+1),io)
          call Mul3(gm1(i),-kcone,g0(i+1),lhm1(i),M1(i),io)
        enddo
        if(b<nblk)then
          call Mul3(g1(b),-kcone,M1(b),lh1(b),rG(ka(me+1)),io)
          call Mul3(gm1(b),-kcone,rG(ka(me+1)),lhm1(b),M1(b),io)
        endif
      endif
      if(need_col)then
        allocate(c1loc(a:b),cKloc(a:b))
        call CopyBlock(c1loc(a),rC1(ka(me)),io)
        do i=a+1,b
          call Mul3(c1loc(i),-kcone,M2(i),lhm1(i-1),c1loc(i-1),io)
        enddo
        call CopyBlock(cKloc(b),rCK(kb(me)),io)
        do i=b-1,a,-1
          call Mul3(cKloc(i),-kcone,M1(i),lh1(i),cKloc(i+1),io)
        enddo
      endif
    endif

    ! results to the master
    if(me==0)then
      if(need_col.or.need_corner) gfout%matdense%a=kczero
      if(have_chunk)then
        do i=a,b
          call SetBlock(g0(i),nb(i),nb(i),off(i),off(i))
          call CopySparseBlocksSingle(gfsparse,g0(i))
          if(need_col)then
            gfout%matdense%a(off(i):off(i)+nb(i)-1,1:nl)=c1loc(i)%a
            gfout%matdense%a(off(i):off(i)+nb(i)-1,nl+1:nl+nr)=cKloc(i)%a
          endif
        enddo
        if(need_offdiag)then
          do i=a,bb
            call SetBlock(g1(i),nb(i),nb(i+1),off(i+1),off(i))
            call CopySparseBlocksSingle(gfsparse,g1(i))
            call SetBlock(gm1(i),nb(i+1),nb(i),off(i),off(i+1))
            call CopySparseBlocksSingle(gfsparse,gm1(i))
          enddo
        endif
      endif
      do p=1,nchunks-1
        cnt=ResultCount(nb,nblk,ca(p),cb(p),nl,nr,need_offdiag,need_col)
        allocate(buf(cnt))
#ifdef MPI
        call MPI_Recv(buf(1),cnt,DAT_dcomplex,p,tag_results,inverse_comm,istatus,mpierror)
#endif
        pos=0
        do i=ca(p),cb(p)
          call UnpackBlock(buf,pos,blk,nb(i),nb(i),off(i),off(i))
          call CopySparseBlocksSingle(gfsparse,blk)
          call FreeBlock(blk)
        enddo
        if(need_offdiag)then
          do i=ca(p),min(cb(p),nblk-1)
            call UnpackBlock(buf,pos,blk,nb(i),nb(i+1),off(i+1),off(i))
            call CopySparseBlocksSingle(gfsparse,blk)
            call FreeBlock(blk)
            call UnpackBlock(buf,pos,blk,nb(i+1),nb(i),off(i),off(i+1))
            call CopySparseBlocksSingle(gfsparse,blk)
            call FreeBlock(blk)
          enddo
        endif
        if(need_col)then
          do i=ca(p),cb(p)
            call UnpackBlock(buf,pos,blk,nb(i),nl,0,0)
            gfout%matdense%a(off(i):off(i)+nb(i)-1,1:nl)=blk%a
            call FreeBlock(blk)
            call UnpackBlock(buf,pos,blk,nb(i),nr,0,0)
            gfout%matdense%a(off(i):off(i)+nb(i)-1,nl+1:nl+nr)=blk%a
            call FreeBlock(blk)
          enddo
        endif
        deallocate(buf)
      enddo
      if(need_corner) gfout%matdense%a(1:nr,1:nl)=rC1(nred)%a
    elseif(have_chunk)then
      cnt=ResultCount(nb,nblk,a,b,nl,nr,need_offdiag,need_col)
      allocate(buf(cnt))
      pos=0
      do i=a,b
        call PackBlock(buf,pos,g0(i))
      enddo
      if(need_offdiag)then
        do i=a,bb
          call PackBlock(buf,pos,g1(i))
          call PackBlock(buf,pos,gm1(i))
        enddo
      endif
      if(need_col)then
        do i=a,b
          call PackBlock(buf,pos,c1loc(i))
          call PackBlock(buf,pos,cKloc(i))
        enddo
      endif
#ifdef MPI
      call MPI_Send(buf(1),cnt,DAT_dcomplex,0,tag_results,inverse_comm,mpierror)
#endif
      deallocate(buf)
    endif

    ! release everything
    if(have_chunk)then
      do i=a,b
        call FreeBlock(lh0(i))
        call FreeBlock(sL(i))
        call FreeBlock(M1(i))
        call FreeBlock(M2(i))
        call FreeBlock(g0(i))
      enddo
      do i=a,bb
        call FreeBlock(lh1(i))
        call FreeBlock(lhm1(i))
      enddo
      deallocate(lh0,lh1,lhm1,sL,M1,M2,g0)
      if(need_offdiag)then
        do i=a,bb
          call FreeBlock(g1(i))
          call FreeBlock(gm1(i))
        enddo
        deallocate(g1,gm1)
      endif
      if(need_col)then
        do i=a,b
          call FreeBlock(c1loc(i))
          call FreeBlock(cKloc(i))
        enddo
        deallocate(c1loc,cKloc)
      endif
      do k=1,nred
        call FreeBlock(rSigL(k))
        call FreeBlock(rSigR(k))
        call FreeBlock(rG(k))
        if(allocated(rC1)) call FreeBlock(rC1(k))
        if(allocated(rCK)) call FreeBlock(rCK(k))
      enddo
      do k=1,nred-1
        call FreeBlock(rMM(k))
        call FreeBlock(rML(k))
      enddo
      deallocate(rSigL,rSigR,rG,rMM,rML)
      if(allocated(rC1)) deallocate(rC1)
      if(allocated(rCK)) deallocate(rCK)
      call FreeBlock(Saa)
      call FreeBlock(Sab)
      call FreeBlock(Sba)
      call FreeBlock(Sbb)
    endif
    do k=1,nred
      call FreeBlock(rD(k))
    enddo
    do k=1,nred-1
      call FreeBlock(rU(k))
      call FreeBlock(rL(k))
    enddo
    deallocate(rD,rU,rL)
    deallocate(nb,off,ca,cb,kept,ka,kb)

  end subroutine InvertSparseONDistributed

  !> contiguous chunk assignment balancing sum(n_i^3); every chunk gets at least one block,
  !> ranks beyond the number of blocks get an empty chunk (ca > cb)
  subroutine AssignChunks(nb,nblk,np,nchunks,ca,cb)
    integer, intent(in) :: nb(:),nblk,np
    integer, intent(out) :: nchunks
    integer, intent(out) :: ca(0:),cb(0:)
    real(kdp) :: total,acc,target
    integer :: i,p

    nchunks=min(np,nblk)
    ca=1
    cb=0
    total=0.0_kdp
    do i=1,nblk
      total=total+real(nb(i),kdp)**3
    enddo
    p=0
    ca(0)=1
    acc=0.0_kdp
    target=total/real(nchunks,kdp)
    do i=1,nblk
      acc=acc+real(nb(i),kdp)**3
      if(p<nchunks-1)then
        if(acc>=target.or.nblk-i==nchunks-1-p)then
          cb(p)=i
          total=total-acc
          p=p+1
          ca(p)=i+1
          acc=0.0_kdp
          target=total/real(nchunks-p,kdp)
        endif
      endif
    enddo
    cb(nchunks-1)=nblk
  end subroutine AssignChunks

  real(kdp) function ChunkCost(nb,a,b)
    integer, intent(in) :: nb(:),a,b
    integer :: i
    ChunkCost=0.0_kdp
    do i=a,b
      ChunkCost=ChunkCost+real(nb(i),kdp)**3
    enddo
  end function ChunkCost

  integer function ChunkBlockCount(nb,nblk,a,b)
    integer, intent(in) :: nb(:),nblk,a,b
    integer :: i
    ChunkBlockCount=0
    do i=a,b
      ChunkBlockCount=ChunkBlockCount+nb(i)**2
    enddo
    do i=a,min(b,nblk-1)
      ChunkBlockCount=ChunkBlockCount+2*nb(i)*nb(i+1)
    enddo
  end function ChunkBlockCount

  integer function ResultCount(nb,nblk,a,b,nl,nr,need_offdiag,need_col)
    integer, intent(in) :: nb(:),nblk,a,b,nl,nr
    logical, intent(in) :: need_offdiag,need_col
    integer :: i
    ResultCount=0
    do i=a,b
      ResultCount=ResultCount+nb(i)**2
      if(need_col) ResultCount=ResultCount+nb(i)*(nl+nr)
    enddo
    if(need_offdiag)then
      do i=a,min(b,nblk-1)
        ResultCount=ResultCount+2*nb(i)*nb(i+1)
      enddo
    endif
  end function ResultCount

  !> Schur complement of the chunk a..b onto its boundary blocks {a,b}; Sab, Sba, Sbb are set only for b > a
  subroutine SchurChunk(lh0,lh1,lhm1,a,b,Saa,Sab,Sba,Sbb,io)
    integer, intent(in) :: a,b
    type(matrixType), intent(inout) :: lh0(a:),lh1(a:),lhm1(a:)
    type(matrixType), intent(inout) :: Saa,Sab,Sba,Sbb
    type(ioType), intent(inout) :: io

    type(matrixType) :: M,sig,gL,gR,col,fc,tmp
    type(matrixType), allocatable :: mm(:),ml(:)
    integer :: i

    call CopyBlock(Saa,lh0(a),io)
    if(b==a) return
    if(b==a+1)then
      call CopyBlock(Sab,lh1(a),io)
      call CopyBlock(Sba,lhm1(a),io)
      call CopyBlock(Sbb,lh0(b),io)
      return
    endif

    allocate(mm(a+1:b-2),ml(a+1:b-2))
    call ZeroBlock(sig,lh0(a+1)%iRows)
    do i=a+1,b-2
      call InvMinus(M,lh0(i),sig,io)
      call Mul(mm(i),-kcone,M,lh1(i),io)
      call FreeBlock(sig)
      call Mul(sig,-kcone,lhm1(i),mm(i),io)
      call FreeBlock(M)
    enddo
    call InvMinus(gL,lh0(b-1),sig,io)
    call FreeBlock(sig)
    call Mul3(tmp,kcone,lhm1(b-1),gL,lh1(b-1),io)
    call AllocLike(Sbb,lh0(b),io)
    call MatrixAdd(Sbb,kcone,lh0(b),-kcone,tmp,io)
    call FreeBlock(tmp)
    call move_alloc(gL%a,col%a)
    col%iRows=gL%iRows
    col%iCols=gL%iCols
    do i=b-2,a+1,-1
      call Mul(tmp,kcone,mm(i),col,io)
      call FreeBlock(col)
      call move_alloc(tmp%a,col%a)
      col%iRows=tmp%iRows
      col%iCols=tmp%iCols
    enddo
    call Mul3(Sab,-kcone,lh1(a),col,lh1(b-1),io)
    call FreeBlock(col)

    call ZeroBlock(sig,lh0(b-1)%iRows)
    do i=b-1,a+2,-1
      call InvMinus(M,lh0(i),sig,io)
      call Mul(ml(i-1),-kcone,M,lhm1(i-1),io)
      call FreeBlock(sig)
      call Mul(sig,-kcone,lh1(i-1),ml(i-1),io)
      call FreeBlock(M)
    enddo
    call InvMinus(gR,lh0(a+1),sig,io)
    call FreeBlock(sig)
    call Mul3(tmp,kcone,lh1(a),gR,lhm1(a),io)
    call FreeBlock(Saa)
    call AllocLike(Saa,lh0(a),io)
    call MatrixAdd(Saa,kcone,lh0(a),-kcone,tmp,io)
    call FreeBlock(tmp)
    call move_alloc(gR%a,fc%a)
    fc%iRows=gR%iRows
    fc%iCols=gR%iCols
    do i=a+2,b-1
      call Mul(tmp,kcone,ml(i-1),fc,io)
      call FreeBlock(fc)
      call move_alloc(tmp%a,fc%a)
      fc%iRows=tmp%iRows
      fc%iCols=tmp%iCols
    enddo
    call Mul3(Sba,-kcone,lhm1(b-1),fc,lhm1(a),io)
    call FreeBlock(fc)

    do i=a+1,b-2
      call FreeBlock(mm(i))
      call FreeBlock(ml(i))
    enddo
    deallocate(mm,ml)
  end subroutine SchurChunk

  !> res = (A - S)^-1
  subroutine InvMinus(res,A,S,io)
    type(matrixType), intent(inout) :: res,A,S
    type(ioType), intent(inout) :: io
    call AllocLike(res,A,io)
    call MatrixAdd(res,kcone,A,-kcone,S,io)
    call InverseMatrix(res,io)
  end subroutine InvMinus

  !> res = (A - S1 - S2)^-1
  subroutine InvMinus2(res,A,S1,S2,io)
    type(matrixType), intent(inout) :: res,A,S1,S2
    type(ioType), intent(inout) :: io
    call AllocLike(res,A,io)
    call MatrixAdd(res,kcone,A,-kcone,S1,-kcone,S2,io)
    call InverseMatrix(res,io)
  end subroutine InvMinus2

  !> res = alpha*A*B
  subroutine Mul(res,alpha,A,B,io)
    type(matrixType), intent(inout) :: res,A,B
    complex(kdp), intent(in) :: alpha
    type(ioType), intent(inout) :: io
    call AllocateMatrix(A%iRows,B%iCols,res,"mInverseDistributed",io)
    call ProductCeAxB(res,alpha,A,B,io)
  end subroutine Mul

  !> res = alpha*A*B*C
  subroutine Mul3(res,alpha,A,B,C,io)
    type(matrixType), intent(inout) :: res,A,B,C
    complex(kdp), intent(in) :: alpha
    type(ioType), intent(inout) :: io
    type(matrixType) :: t
    call Mul(t,kcone,A,B,io)
    call Mul(res,alpha,t,C,io)
    call FreeBlock(t)
  end subroutine Mul3

  subroutine AllocLike(dst,src,io)
    type(matrixType), intent(inout) :: dst,src
    type(ioType), intent(inout) :: io
    call AllocateMatrix(src%iRows,src%iCols,src%iHorz,src%iVert,dst,"mInverseDistributed",io)
  end subroutine AllocLike

  subroutine CopyBlock(dst,src,io)
    type(matrixType), intent(inout) :: dst,src
    type(ioType), intent(inout) :: io
    call AllocLike(dst,src,io)
    dst%a=src%a
  end subroutine CopyBlock

  subroutine ZeroBlock(blk,n)
    type(matrixType), intent(inout) :: blk
    integer, intent(in) :: n
    allocate(blk%a(n,n))
    blk%a=kczero
    blk%iRows=n
    blk%iCols=n
    blk%iHorz=0
    blk%iVert=0
  end subroutine ZeroBlock

  subroutine SetBlock(blk,rows,cols,horz,vert)
    type(matrixType), intent(inout) :: blk
    integer, intent(in) :: rows,cols,horz,vert
    blk%iRows=rows
    blk%iCols=cols
    blk%iHorz=horz
    blk%iVert=vert
  end subroutine SetBlock

  subroutine FreeBlock(blk)
    type(matrixType), intent(inout) :: blk
    if(allocated(blk%a)) deallocate(blk%a)
  end subroutine FreeBlock

  subroutine PackBlock(buf,pos,blk)
    complex(kdp), intent(inout) :: buf(:)
    integer, intent(inout) :: pos
    type(matrixType), intent(in) :: blk
    integer :: n
    n=blk%iRows*blk%iCols
    buf(pos+1:pos+n)=reshape(blk%a,(/n/))
    pos=pos+n
  end subroutine PackBlock

  subroutine UnpackBlock(buf,pos,blk,rows,cols,horz,vert)
    complex(kdp), intent(in) :: buf(:)
    integer, intent(inout) :: pos
    type(matrixType), intent(inout) :: blk
    integer, intent(in) :: rows,cols,horz,vert
    allocate(blk%a(rows,cols))
    blk%a=reshape(buf(pos+1:pos+rows*cols),(/rows,cols/))
    call SetBlock(blk,rows,cols,horz,vert)
    pos=pos+rows*cols
  end subroutine UnpackBlock

end module mInverseDistributed
