!
! Distributed block-tridiagonal inversion over the inverse_comm process group.
!
! Contract (shared with InvertSparseONv3): gfsparse holds E*S-H-Sigma in CRS form; on return its
! diagonal blocks (and the first off-diagonal blocks for opindex 1,3,5) contain the corresponding
! blocks of G on the stored pattern; gfout(:,1:nl) = G(:,1:nl) and gfout(:,nl+1:nl+nr) =
! G(:,N1-nr+1:N1) for opindex 2,3; gfout(1:nr,1:nl) = G(N1-nr+1:N1,1:nl) for opindex 4,5.
! Every rank of inverse_comm must call.
!   replicated = .false.: gfsparse exists on mynode_inverse == 0 only; it scatters the blocks and
!                         gathers every result back (gfout on the master only).
!   replicated = .true. : every rank holds the same gfsparse (pattern and values) and fills only
!                         the blocks of its own chunk; it receives the G blocks of its chunk in its
!                         own gfsparse and, for opindex 2,3, the lead-column rows of its chunk and of
!                         the next chunk's first block in its own gfout (gfout must exist on every
!                         rank); gather = .true. additionally delivers everything to the master.
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
  public :: DistributedEntryOwnerMask
  public :: DistributedLocalRowRange
  public :: ReduceEnergySliceToMaster

  integer, parameter :: tag_blocks = 3101
  integer, parameter :: tag_results = 3102
  integer, parameter :: tag_columns = 3103
  logical, save :: first_call = .true.

  ! block layout of the pattern inverted last (pattern arrays kept only when derived from a local copy)
  integer, save :: c_n1 = -1, c_nnz = -1, c_nblk = 0, c_nchunks = 0
  integer, allocatable, save :: c_q(:), c_j(:), c_nb(:), c_off(:), c_ca(:), c_cb(:), c_blkofrow(:), c_ownerofblk(:)

contains

  logical function DistributedInversionActive(solver)
    integer, intent(in) :: solver
    DistributedInversionActive = (solver == 2 .and. nnodes_inverse > 1)
  end function DistributedInversionActive

  subroutine InvertSparseONDistributed(N1,gfsparse,nl,nr,gfout,opindex,replicated,gather)
    integer, intent(in) :: N1,nl,nr,opindex
    type(matrixTypeGeneral), intent(inout) :: gfsparse,gfout
    logical, intent(in) :: replicated,gather

    character(len=*), parameter :: sMyName="InvertSparseONDistributed"
    type(ioType) :: io
    type(matrixType), allocatable :: lh0(:),lh1(:),lhm1(:),th0(:),th1(:),thm1(:)
    type(matrixType), allocatable :: rD(:),rU(:),rL(:),rSigL(:),rSigR(:),rG(:),rC1(:),rCK(:),rMM(:),rML(:)
    type(matrixType), allocatable :: sL(:),M1(:),M2(:),g0(:),g1(:),gm1(:),c1loc(:),cKloc(:)
    type(matrixType) :: Saa,Sab,Sba,Sbb,sR,M,blk
    complex(kdp), allocatable :: buf(:)
    integer :: me,np,nblk,nchunks,a,b,bb,i,p,k,nred,cnt,pos,mpierror,an
    logical :: need_offdiag,need_col,need_corner,have_chunk,gath,write_local
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
    gath=gather.or..not.replicated

    call AcquireLayout(gfsparse,nl,nr,N1,replicated,io)
    nblk=c_nblk
    nchunks=c_nchunks
    a=c_ca(me)
    b=c_cb(me)
    have_chunk=(a<=b)
    bb=min(b,nblk-1)
    write_local=(me==0).or.(replicated.and..not.gath)
    nred=0
    do p=0,nchunks-1
      nred=nred+c_kept(p)
    enddo

    if(first_call.and.me==0)then
      total=0.0_kdp
      do i=1,nblk
        total=total+real(c_nb(i),kdp)**3
      enddo
      write(negf_log_unit,'(a,i0,a,i0,a,i0,a,i0,a,l1)') 'InvertSparseONDistributed: N1=',N1,' blocks=',nblk, &
        ' inverse ranks=',np,' chunks=',nchunks,' replicated matrix=',replicated
      do p=0,nchunks-1
        write(negf_log_unit,'(a,i0,a,i0,a,i0,a,f7.4)') '  rank ',p,': blocks ',c_ca(p),'..',c_cb(p), &
          ' cost share ',ChunkCost(c_nb,c_ca(p),c_cb(p))/total
      enddo
      first_call=.false.
    endif

    ! blocks of the own chunk: from the local copy, or scattered by the master
    if(have_chunk) allocate(lh0(a:b),lh1(a:bb),lhm1(a:bb))
    if(replicated)then
      if(have_chunk) call FillChunkBlocks(gfsparse,a,b,lh0,lh1,lhm1,io)
    else
      if(me==0)then
        do p=1,nchunks-1
          allocate(th0(c_ca(p):c_cb(p)),th1(c_ca(p):min(c_cb(p),nblk-1)),thm1(c_ca(p):min(c_cb(p),nblk-1)))
          call FillChunkBlocks(gfsparse,c_ca(p),c_cb(p),th0,th1,thm1,io)
          cnt=ChunkBlockCount(c_nb,nblk,c_ca(p),c_cb(p))
          allocate(buf(cnt))
          pos=0
          do i=c_ca(p),c_cb(p)
            call PackBlock(buf,pos,th0(i))
            call FreeBlock(th0(i))
          enddo
          do i=c_ca(p),min(c_cb(p),nblk-1)
            call PackBlock(buf,pos,th1(i))
            call PackBlock(buf,pos,thm1(i))
            call FreeBlock(th1(i))
            call FreeBlock(thm1(i))
          enddo
          deallocate(th0,th1,thm1)
#ifdef MPI
          call MPI_Send(buf(1),cnt,DAT_dcomplex,p,tag_blocks,inverse_comm,mpierror)
#endif
          deallocate(buf)
        enddo
        call FillChunkBlocks(gfsparse,a,b,lh0,lh1,lhm1,io)
      elseif(have_chunk)then
        cnt=ChunkBlockCount(c_nb,nblk,a,b)
        allocate(buf(cnt))
#ifdef MPI
        call MPI_Recv(buf(1),cnt,DAT_dcomplex,0,tag_blocks,inverse_comm,istatus,mpierror)
#endif
        pos=0
        do i=a,b
          call UnpackBlock(buf,pos,lh0(i),c_nb(i),c_nb(i),c_off(i),c_off(i))
        enddo
        do i=a,bb
          call UnpackBlock(buf,pos,lh1(i),c_nb(i),c_nb(i+1),c_off(i+1),c_off(i))
          call UnpackBlock(buf,pos,lhm1(i),c_nb(i+1),c_nb(i),c_off(i),c_off(i+1))
        enddo
        deallocate(buf)
      endif
    endif

    ! Schur complement of every chunk onto its boundary blocks
    if(have_chunk) call SchurChunk(lh0,lh1,lhm1,a,b,Saa,Sab,Sba,Sbb,io)

    ! reduced block-tridiagonal system, replicated on every rank
    allocate(rD(nred),rU(max(nred-1,0)),rL(max(nred-1,0)))
    do p=0,nchunks-1
      cnt=c_nb(c_ca(p))**2
      if(c_kept(p)==2) cnt=cnt+2*c_nb(c_ca(p))*c_nb(c_cb(p))+c_nb(c_cb(p))**2
      if(c_cb(p)<nblk) cnt=cnt+2*c_nb(c_cb(p))*c_nb(c_cb(p)+1)
      allocate(buf(cnt))
      if(p==me)then
        pos=0
        call PackBlock(buf,pos,Saa)
        if(c_kept(p)==2)then
          call PackBlock(buf,pos,Sab)
          call PackBlock(buf,pos,Sba)
          call PackBlock(buf,pos,Sbb)
        endif
        if(c_cb(p)<nblk)then
          call PackBlock(buf,pos,lh1(b))
          call PackBlock(buf,pos,lhm1(b))
        endif
      endif
#ifdef MPI
      call MPI_Bcast(buf(1),cnt,DAT_dcomplex,p,inverse_comm,mpierror)
#endif
      pos=0
      k=c_ka(p)
      call UnpackBlock(buf,pos,rD(k),c_nb(c_ca(p)),c_nb(c_ca(p)),0,0)
      if(c_kept(p)==2)then
        call UnpackBlock(buf,pos,rU(k),c_nb(c_ca(p)),c_nb(c_cb(p)),0,0)
        call UnpackBlock(buf,pos,rL(k),c_nb(c_cb(p)),c_nb(c_ca(p)),0,0)
        call UnpackBlock(buf,pos,rD(k+1),c_nb(c_cb(p)),c_nb(c_cb(p)),0,0)
      endif
      if(c_cb(p)<nblk)then
        call UnpackBlock(buf,pos,rU(c_kb(p)),c_nb(c_cb(p)),c_nb(c_cb(p)+1),0,0)
        call UnpackBlock(buf,pos,rL(c_kb(p)),c_nb(c_cb(p)+1),c_nb(c_cb(p)),0,0)
      endif
      deallocate(buf)
    enddo

    allocate(rSigL(nred),rSigR(nred),rG(nred),rMM(max(nred-1,0)),rML(max(nred-1,0)))
    if(.true.)then
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
    endif

    if(have_chunk)then
      ! chunk solve seeded with the exact boundary self-energies
      allocate(sL(a:b),M1(a:b),M2(a:b),g0(a:b))
      call CopyBlock(sL(a),rSigL(c_ka(me)),io)
      do i=a,b-1
        call InvMinus(M1(i),lh0(i),sL(i),io)
        call Mul3(sL(i+1),kcone,lhm1(i),M1(i),lh1(i),io)
      enddo
      if(need_offdiag.and.b<nblk) call InvMinus(M1(b),lh0(b),sL(b),io)
      call CopyBlock(sR,rSigR(c_kb(me)),io)
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
          maxval(abs(rG(c_ka(me))%a-g0(a)%a)),maxval(abs(rG(c_kb(me))%a-g0(b)%a))
      endif
      if(need_offdiag)then
        allocate(g1(a:bb),gm1(a:bb))
        do i=a,b-1
          call Mul3(g1(i),-kcone,M1(i),lh1(i),g0(i+1),io)
          call Mul3(gm1(i),-kcone,g0(i+1),lhm1(i),M1(i),io)
        enddo
        if(b<nblk)then
          call Mul3(g1(b),-kcone,M1(b),lh1(b),rG(c_ka(me+1)),io)
          call Mul3(gm1(b),-kcone,rG(c_ka(me+1)),lhm1(b),M1(b),io)
        endif
      endif
      if(need_col)then
        allocate(c1loc(a:b),cKloc(a:b))
        call CopyBlock(c1loc(a),rC1(c_ka(me)),io)
        do i=a+1,b
          call Mul3(c1loc(i),-kcone,M2(i),lhm1(i-1),c1loc(i-1),io)
        enddo
        call CopyBlock(cKloc(b),rCK(c_kb(me)),io)
        do i=b-1,a,-1
          call Mul3(cKloc(i),-kcone,M1(i),lh1(i),cKloc(i+1),io)
        enddo
      endif
    endif

    ! own results into the local containers
    if(write_local.and.(need_col.or.need_corner)) gfout%matdense%a=kczero
    if(write_local.and.have_chunk)then
      do i=a,b
        call SetBlock(g0(i),c_nb(i),c_nb(i),c_off(i),c_off(i))
        call CopySparseBlocksSingle(gfsparse,g0(i))
        if(need_col)then
          gfout%matdense%a(c_off(i):c_off(i)+c_nb(i)-1,1:nl)=c1loc(i)%a
          gfout%matdense%a(c_off(i):c_off(i)+c_nb(i)-1,nl+1:nl+nr)=cKloc(i)%a
        endif
      enddo
      if(need_offdiag)then
        do i=a,bb
          call SetBlock(g1(i),c_nb(i),c_nb(i+1),c_off(i+1),c_off(i))
          call CopySparseBlocksSingle(gfsparse,g1(i))
          call SetBlock(gm1(i),c_nb(i+1),c_nb(i),c_off(i),c_off(i+1))
          call CopySparseBlocksSingle(gfsparse,gm1(i))
        enddo
      endif
    endif
    if(write_local.and.need_corner) gfout%matdense%a(1:nr,1:nl)=rC1(nred)%a

    ! lead-column rows of the next chunk's first block, needed by the owner of the inter-chunk pair
    if(replicated.and..not.gath.and.need_col.and.have_chunk)then
      if(me>0)then
        cnt=c_nb(a)*(nl+nr)
        allocate(buf(cnt))
        pos=0
        call PackBlock(buf,pos,c1loc(a))
        call PackBlock(buf,pos,cKloc(a))
#ifdef MPI
        call MPI_Send(buf(1),cnt,DAT_dcomplex,me-1,tag_columns,inverse_comm,mpierror)
#endif
        deallocate(buf)
      endif
      if(b<nblk)then
        an=c_ca(me+1)
        cnt=c_nb(an)*(nl+nr)
        allocate(buf(cnt))
#ifdef MPI
        call MPI_Recv(buf(1),cnt,DAT_dcomplex,me+1,tag_columns,inverse_comm,istatus,mpierror)
#endif
        pos=0
        call UnpackBlock(buf,pos,blk,c_nb(an),nl,0,0)
        gfout%matdense%a(c_off(an):c_off(an)+c_nb(an)-1,1:nl)=blk%a
        call FreeBlock(blk)
        call UnpackBlock(buf,pos,blk,c_nb(an),nr,0,0)
        gfout%matdense%a(c_off(an):c_off(an)+c_nb(an)-1,nl+1:nl+nr)=blk%a
        call FreeBlock(blk)
        deallocate(buf)
      endif
    endif

    ! all results to the master
    if(gath)then
      if(me==0)then
        do p=1,nchunks-1
          cnt=ResultCount(c_nb,nblk,c_ca(p),c_cb(p),nl,nr,need_offdiag,need_col)
          allocate(buf(cnt))
#ifdef MPI
          call MPI_Recv(buf(1),cnt,DAT_dcomplex,p,tag_results,inverse_comm,istatus,mpierror)
#endif
          pos=0
          do i=c_ca(p),c_cb(p)
            call UnpackBlock(buf,pos,blk,c_nb(i),c_nb(i),c_off(i),c_off(i))
            call CopySparseBlocksSingle(gfsparse,blk)
            call FreeBlock(blk)
          enddo
          if(need_offdiag)then
            do i=c_ca(p),min(c_cb(p),nblk-1)
              call UnpackBlock(buf,pos,blk,c_nb(i),c_nb(i+1),c_off(i+1),c_off(i))
              call CopySparseBlocksSingle(gfsparse,blk)
              call FreeBlock(blk)
              call UnpackBlock(buf,pos,blk,c_nb(i+1),c_nb(i),c_off(i),c_off(i+1))
              call CopySparseBlocksSingle(gfsparse,blk)
              call FreeBlock(blk)
            enddo
          endif
          if(need_col)then
            do i=c_ca(p),c_cb(p)
              call UnpackBlock(buf,pos,blk,c_nb(i),nl,0,0)
              gfout%matdense%a(c_off(i):c_off(i)+c_nb(i)-1,1:nl)=blk%a
              call FreeBlock(blk)
              call UnpackBlock(buf,pos,blk,c_nb(i),nr,0,0)
              gfout%matdense%a(c_off(i):c_off(i)+c_nb(i)-1,nl+1:nl+nr)=blk%a
              call FreeBlock(blk)
            enddo
          endif
          deallocate(buf)
        enddo
      elseif(have_chunk)then
        cnt=ResultCount(c_nb,nblk,a,b,nl,nr,need_offdiag,need_col)
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
      call FreeBlock(Saa)
      call FreeBlock(Sab)
      call FreeBlock(Sba)
      call FreeBlock(Sbb)
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
    do k=1,nred
      call FreeBlock(rD(k))
    enddo
    do k=1,nred-1
      call FreeBlock(rU(k))
      call FreeBlock(rL(k))
    enddo
    deallocate(rD,rU,rL)

  end subroutine InvertSparseONDistributed

  !> ownership of the stored entries of a matrix with the cached block layout: entry (i,j) belongs to the
  !> owner of the lower of the two block indices (the chunk that computed that diagonal block or the
  !> (k,k+1)/(k+1,k) pair); valid after an inversion of a replicated matrix with the same row count
  subroutine DistributedEntryOwnerMask(n1,q,j,nnz,mask)
    integer, intent(in) :: n1,nnz
    integer, intent(in) :: q(n1+1),j(nnz)
    logical, intent(out) :: mask(nnz)
    integer :: i,ind,bi,bj

    if(c_n1/=n1) call negf_abort("DistributedEntryOwnerMask: no block layout of a replicated matrix with this row count")
    do i=1,n1
      bi=c_blkofrow(i)
      do ind=q(i),q(i+1)-1
        bj=c_blkofrow(j(ind))
        mask(ind)=(c_ownerofblk(min(bi,bj))==mynode_inverse)
      enddo
    enddo
  end subroutine DistributedEntryOwnerMask

  !> rows whose lead-column entries the local mode fills on this rank: its own chunk and the first block
  !> of the next chunk (r0 > r1 for a rank without a chunk); valid after an inversion of a replicated matrix
  subroutine DistributedLocalRowRange(r0,r1)
    integer, intent(out) :: r0,r1
    integer :: a,b,last

    r0=1
    r1=0
    if(c_nblk<=0) return
    a=c_ca(mynode_inverse)
    b=c_cb(mynode_inverse)
    if(a>b) return
    last=b
    if(b<c_nblk) last=b+1
    r0=c_off(a)
    r1=c_off(last)+c_nb(last)-1
  end subroutine DistributedLocalRowRange

  !> sums buf(ispin,ie,:) over inverse_comm onto the master (the other ranks keep their partial slice)
  subroutine ReduceEnergySliceToMaster(buf,nspin,ne,nnz,ispin,ie)
    integer, intent(in) :: nspin,ne,nnz,ispin,ie
    complex(kdp), intent(inout) :: buf(nspin,ne,nnz)
    complex(kdp), allocatable :: tmp(:),res(:)
    integer :: mpierror

#ifdef MPI
    if(nnodes_inverse<=1.or.nnz<=0) return
    allocate(tmp(nnz),res(nnz))
    tmp=buf(ispin,ie,:)
    res=kczero
    call MPI_Reduce(tmp(1),res(1),nnz,DAT_dcomplex,MPI_SUM,0,inverse_comm,mpierror)
    if(mynode_inverse==0) buf(ispin,ie,:)=res
    deallocate(tmp,res)
#endif
  end subroutine ReduceEnergySliceToMaster

  !> block layout for this call: derived from the local copy (replicated) or by the master and broadcast
  subroutine AcquireLayout(gfsparse,nl,nr,N1,replicated,io)
    type(matrixTypeGeneral), intent(in) :: gfsparse
    integer, intent(in) :: nl,nr,N1
    logical, intent(in) :: replicated
    type(ioType), intent(inout) :: io
    character(len=*), parameter :: sMyName="InvertSparseONDistributed"
    integer, allocatable :: nb(:),off(:)
    integer :: nblk,mpierror
    logical :: same

    if(mynode_inverse==0.or.replicated)then
      if(gfsparse%mattype/=2) call negf_abort(sMyName//": the distributed inverter needs the sparse (EM.OrderN) Green function matrix")
      same=(c_n1==gfsparse%iRows.and.c_nnz==gfsparse%matSparse%nnz.and.allocated(c_q))
      if(same) same=all(c_q==gfsparse%matSparse%q(1:c_n1+1)).and.all(c_j==gfsparse%matSparse%j(1:c_nnz))
      if(.not.same)then
        call PartitionBlockLayout(gfsparse,nl,nr,nblk,nb,off,io)
        if(off(1)/=1.or.off(nblk)+nb(nblk)-1/=N1.or.nb(1)/=nl.or.nb(nblk)/=nr) &
          call negf_abort(sMyName//": the first and last blocks must be the lead blocks (nl, nr) and cover 1..N1")
        call StoreLayout(nblk,nb,off,N1)
        c_n1=gfsparse%iRows
        c_nnz=gfsparse%matSparse%nnz
        if(allocated(c_q)) deallocate(c_q)
        if(allocated(c_j)) deallocate(c_j)
        allocate(c_q(c_n1+1),c_j(c_nnz))
        c_q=gfsparse%matSparse%q(1:c_n1+1)
        c_j=gfsparse%matSparse%j(1:c_nnz)
      endif
    endif
    if(.not.replicated)then
#ifdef MPI
      nblk=c_nblk
      call MPI_Bcast(nblk,1,MPI_integer,0,inverse_comm,mpierror)
      if(allocated(nb)) deallocate(nb,off)
      allocate(nb(nblk),off(nblk))
      if(mynode_inverse==0)then
        nb=c_nb
        off=c_off
      endif
      call MPI_Bcast(nb(1),nblk,MPI_integer,0,inverse_comm,mpierror)
      call MPI_Bcast(off(1),nblk,MPI_integer,0,inverse_comm,mpierror)
      if(mynode_inverse/=0)then
        call StoreLayout(nblk,nb,off,N1)
        c_n1=-1
        c_nnz=-1
      endif
#endif
    endif
  end subroutine AcquireLayout

  subroutine StoreLayout(nblk,nb,off,N1)
    integer, intent(in) :: nblk,N1
    integer, intent(in) :: nb(nblk),off(nblk)
    integer :: i,p,np

    np=nnodes_inverse
    if(allocated(c_nb)) deallocate(c_nb,c_off,c_ca,c_cb,c_blkofrow,c_ownerofblk)
    allocate(c_nb(nblk),c_off(nblk),c_ca(0:np-1),c_cb(0:np-1),c_blkofrow(N1),c_ownerofblk(nblk))
    c_nblk=nblk
    c_nb=nb
    c_off=off
    call AssignChunks(nb,nblk,np,c_nchunks,c_ca,c_cb)
    c_ownerofblk=-1
    do p=0,c_nchunks-1
      do i=c_ca(p),c_cb(p)
        c_ownerofblk(i)=p
      enddo
    enddo
    do i=1,nblk
      c_blkofrow(off(i):off(i)+nb(i)-1)=i
    enddo
  end subroutine StoreLayout

  integer function c_kept(p)
    integer, intent(in) :: p
    c_kept=2
    if(c_ca(p)==c_cb(p)) c_kept=1
  end function c_kept

  !> index of the first reduced block of chunk p
  integer function c_ka(p)
    integer, intent(in) :: p
    integer :: ip
    c_ka=1
    do ip=0,p-1
      c_ka=c_ka+c_kept(ip)
    enddo
  end function c_ka

  !> index of the last reduced block of chunk p
  integer function c_kb(p)
    integer, intent(in) :: p
    c_kb=c_ka(p)+c_kept(p)-1
  end function c_kb

  !> dense blocks of the chunk a..b (h0 for a..b, h1/hm1 for a..min(b,nblk-1)) filled from a CRS matrix
  subroutine FillChunkBlocks(gfsparse,a,b,h0,h1,hm1,io)
    type(matrixTypeGeneral), intent(in) :: gfsparse
    integer, intent(in) :: a,b
    type(matrixType), intent(inout) :: h0(a:),h1(a:),hm1(a:)
    type(ioType), intent(inout) :: io
    integer :: i

    do i=a,b
      call NewBlock(h0(i),c_nb(i),c_nb(i),c_off(i),c_off(i),io)
      call FillBlock(gfsparse,h0(i))
    enddo
    do i=a,min(b,c_nblk-1)
      call NewBlock(h1(i),c_nb(i),c_nb(i+1),c_off(i+1),c_off(i),io)
      call FillBlock(gfsparse,h1(i))
      call NewBlock(hm1(i),c_nb(i+1),c_nb(i),c_off(i),c_off(i+1),io)
      call FillBlock(gfsparse,hm1(i))
    enddo
  end subroutine FillChunkBlocks

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

  subroutine NewBlock(blk,rows,cols,horz,vert,io)
    type(matrixType), intent(inout) :: blk
    integer, intent(in) :: rows,cols,horz,vert
    type(ioType), intent(inout) :: io
    call AllocateMatrix(rows,cols,horz,vert,blk,"mInverseDistributed",io)
    blk%a=kczero
  end subroutine NewBlock

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
