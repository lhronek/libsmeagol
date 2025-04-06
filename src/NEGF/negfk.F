! 
! Copyright (c) Authors:
! Ivan Rungger and Stefano Sanvito
! Trinity College Dublin, Ireland
! October 2008 
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! THE SUBROUTINE
!                   NEGFK  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

      SUBROUTINE NEGFK(N1,NspinBlocks,NspinComplexMatrix,nsc,Nenerg,
     &    Nenerg1, Nenerg2,
     &    NPOLES,istep,inicoor,IV,ITER,ik,ICASE,Vini,Vfinal,T,Delta,V,
     &    EB,CONVERGED,NSLICES,CALCTRANSM, CalcIETS,
     &    nk,kpoint,slabel,wk,
     &    hgeneralp,sgeneralp,rhogeneralp,ematgeneralp,ldos,
     .    idyn, tmdskip, tmdsampling)

! *****************************************************************
! Calculates the charge density of a nanoscopic system coupled to
! charge reservoirs. It uses the non-equilibrium Green s
! function approach.
! The equilibrium charge density is calculated using a contour
! integral in the complex plane.
!
! Collinear spins, k points version
!
! On convergence the current is calculated (and possibly the 
! transmission coefficients)
!
! Written by Alexandre Reily Rocha and Victor Garcia-Suarez, June 2003 
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: rochaa@tcd.ie
! Departamento de Fisica
! Universidad de Oviedo
! e-mail: victor@condmat.uniovi.es
! *************************** HISTORY **********************************
! Original version: June 2003
! Altered by Alexandre Reily Rocha and Ivan Rungger, October 2008:
!      Extensive changes to join negfk and negf into one subroutine
!      Modularisation into separate subroutines, mose of these are 
!      found in the file negfk_subs.F
! ***************************** INPUT ***************************** 

      use precision
      use negfmod
      use mTypes
      use mMatrixUtil
      use mBoundstates
      use mMPI_NEGF
      use sigma
      use mEnergyGrid
      use mCurrentDistribution, only : CurrentDistributionMatrix_RhoK
      use ScissorOperator, only : sco,scoSCF,SCOSetHamiltonianBlock,
     .    SCOApplyK_nc

      INCLUDE "negfk_include.h"



      integer, intent(in) :: NspinComplexMatrix,NspinBlocks
      logical, external :: leqi
      integer ii,jj,ind,nl2,nr2,n12
      double complex, allocatable :: sdense2(:,:)
      type(matrixTypeGeneral) :: rhobs_general(NspinComplexMatrix)
      type(matrixTypeGeneral) :: ematbs_general(NspinComplexMatrix)
      type(matrixTypeGeneral) :: hgeneralp(NspinComplexMatrix),sgeneralp
      type(matrixTypeGeneral) :: ematgeneralp(NspinComplexMatrix)
      type(matrixTypeGeneral) :: ematgeneralp_l_real(NspinComplexMatrix)
      type(matrixTypeGeneral) :: ematgeneralp_r_real(NspinComplexMatrix)
      type(matrixTypeGeneral) :: ematgeneralp_l_im(NspinComplexMatrix)
      type(matrixTypeGeneral) :: ematgeneralp_r_im(NspinComplexMatrix)
      type(matrixTypeGeneral) :: rhogeneralp(NspinComplexMatrix)
      type(matrixTypeGeneral) :: rhogeneralp_l_real(NspinComplexMatrix)
      type(matrixTypeGeneral) :: rhogeneralp_r_real(NspinComplexMatrix)
      type(matrixTypeGeneral) :: rhogeneralp_l_im(NspinComplexMatrix)
      type(matrixTypeGeneral) :: rhogeneralp_r_im(NspinComplexMatrix)
      type(ioType) :: io
      integer*4:: sc_0,sc_1,sc_r,sc_m
      integer*4:: sc_0b,sc_1b,sc_rb,sc_mb
      integer*4:: sc_0a,sc_1a
      integer :: nleadslr=2
      integer, allocatable :: leadsdim(:)
      double precision, allocatable :: LeadsVoltageShift(:)
      logical :: ldos, CalcIETS
      integer nspinMin
      DOUBLE PRECISION, DIMENSION (NspinBlocks) :: Icb,tk
      DOUBLE PRECISION, DIMENSION (NspinBlocks) :: Icfinal
      double precision, allocatable :: Rhopdostotk(:,:)
      double precision :: emdostotk(4), weight
      double complex, allocatable :: rho_real_l(:,:,:)
      double complex, allocatable :: rho_real_r(:,:,:)
      double complex, allocatable :: rho_im_l(:,:,:)
      double complex, allocatable :: rho_im_r(:,:,:)
      double complex, allocatable :: rho_f_im_l(:,:,:)
      double complex, allocatable :: rho_f_im_r(:,:,:)


#ifdef MPI
      INTEGER :: MPIerror
#endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
        CALL SYSTEM_CLOCK(sc_0a,sc_rb,sc_mb)
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
        write(12347,'(A,f12.6)')
     $       't_negf1',(sc_1b-sc_0b)*1.0d0/sc_rb
        CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
      endif




 
#ifdef MPI
      CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
      CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
#else
      Nnodes=1
      MyNode=0
#endif
      
      io%isDebug=.false.
      CALL TIMER('NEGFK',1)
      CALL TIMER('NEGFK_su',1)

! Create density matrices for real and imaginary components
            do ispin=1,NspinComplexMatrix

              call AllocateMatrixGeneral(n1,n1,
     .             rhogeneralp(ispin)%matSparse%nnz,2,
     .        rhogeneralp_l_real(ispin),"negfk",io)
              rhogeneralp_l_real(ispin)%matSparse%q(:)=
     .            rhogeneralp(ispin)%matSparse%q(:)
              rhogeneralp_l_real(ispin)%matSparse%j(:)=
     .            rhogeneralp(ispin)%matSparse%j(:)
              rhogeneralp_l_real(ispin)%matSparse%b(:)=0D0

              call AllocateMatrixGeneral(n1,n1,
     .             rhogeneralp(ispin)%matSparse%nnz,2,
     .        rhogeneralp_r_real(ispin),"negfk",io)
              rhogeneralp_r_real(ispin)%matSparse%q(:)=
     .            rhogeneralp(ispin)%matSparse%q(:)
              rhogeneralp_r_real(ispin)%matSparse%j(:)=
     .            rhogeneralp(ispin)%matSparse%j(:)
              rhogeneralp_r_real(ispin)%matSparse%b(:)=0D0

              call AllocateMatrixGeneral(n1,n1,
     .             rhogeneralp(ispin)%matSparse%nnz,2,
     .        rhogeneralp_l_im(ispin),"negfk",io)
              rhogeneralp_l_im(ispin)%matSparse%q(:)=
     .            rhogeneralp(ispin)%matSparse%q(:)
              rhogeneralp_l_im(ispin)%matSparse%j(:)=
     .            rhogeneralp(ispin)%matSparse%j(:)
              rhogeneralp_l_im(ispin)%matSparse%b(:)=0D0

              call AllocateMatrixGeneral(n1,n1,
     .             rhogeneralp(ispin)%matSparse%nnz,2,
     .        rhogeneralp_r_im(ispin),"negfk",io)
              rhogeneralp_r_im(ispin)%matSparse%q(:)=
     .            rhogeneralp(ispin)%matSparse%q(:)
              rhogeneralp_r_im(ispin)%matSparse%j(:)=
     .            rhogeneralp(ispin)%matSparse%j(:)
              rhogeneralp_r_im(ispin)%matSparse%b(:)=0D0

             if(emforces)then

              call AllocateMatrixGeneral(n1,n1,
     .             ematgeneralp(ispin)%matSparse%nnz,2,
     .        ematgeneralp_l_real(ispin),"negfk",io)
              ematgeneralp_l_real(ispin)%matSparse%q(:)=
     .            ematgeneralp(ispin)%matSparse%q(:)
              ematgeneralp_l_real(ispin)%matSparse%j(:)=
     .            ematgeneralp(ispin)%matSparse%j(:)
              ematgeneralp_l_real(ispin)%matSparse%b(:)=0D0

              call AllocateMatrixGeneral(n1,n1,
     .             ematgeneralp(ispin)%matSparse%nnz,2,
     .        ematgeneralp_r_real(ispin),"negfk",io)
              ematgeneralp_r_real(ispin)%matSparse%q(:)=
     .            ematgeneralp(ispin)%matSparse%q(:)
              ematgeneralp_r_real(ispin)%matSparse%j(:)=
     .            ematgeneralp(ispin)%matSparse%j(:)
              ematgeneralp_r_real(ispin)%matSparse%b(:)=0D0

              call AllocateMatrixGeneral(n1,n1,
     .             ematgeneralp(ispin)%matSparse%nnz,2,
     .        ematgeneralp_l_im(ispin),"negfk",io)
              ematgeneralp_l_im(ispin)%matSparse%q(:)=
     .            ematgeneralp(ispin)%matSparse%q(:)
              ematgeneralp_l_im(ispin)%matSparse%j(:)=
     .            ematgeneralp(ispin)%matSparse%j(:)
              ematgeneralp_l_im(ispin)%matSparse%b(:)=0D0

              call AllocateMatrixGeneral(n1,n1,
     .             ematgeneralp(ispin)%matSparse%nnz,2,
     .        ematgeneralp_r_im(ispin),"negfk",io)
              ematgeneralp_r_im(ispin)%matSparse%q(:)=
     .            ematgeneralp(ispin)%matSparse%q(:)
              ematgeneralp_r_im(ispin)%matSparse%j(:)=
     .            ematgeneralp(ispin)%matSparse%j(:)
              ematgeneralp_r_im(ispin)%matSparse%b(:)=0D0

             endif  ! if(emforces)then
            enddo ! do ispin=1,NspinComplexMatrix

!---set negfmod values-------------------------------
      itermod=iter
      ivmod=IV
!----------------------------------------------------

      if(WriteNk.and.iv.eq.0.and.Mynode.eq.0)
     .    write(*,*)"i_kp=",ik,mynode_groupK
     
!  ****** Allocation of Memory for the electrodes*******************
      If ((ITER .EQ. 1) .AND. (istep .EQ. inicoor) .AND. (ik .EQ. 1)
     &    .AND. (IV .EQ. 0)) THEN

! Read some data for the electrodes
        call leads_info(slabelL, NL, nspinL, maxnhL, Ef_LeadL,
     &      slabelR, NR, nspinR, maxnhR, Ef_LeadR,Ef_Lead,NspinBlocks)

      EndIf        
! ****************************************************************

! Read and calculate the Hamiltonians, overlaps  and density matrices
! of the leads


      if (Mynode .eq. 0) Then
!        write(*,*)"hsleads",NspinBlocks, NspinBlocks, NL, NR, nsc, ITER, istep, IV,
!     .      gamma_negf, ik, kpoint, T

        call hsleads(NspinBlocks, NspinBlocks, NL, NR, nsc, ITER,
     .      istep,inicoor,
     .      IV, gamma_negf, ik, kpoint, T,
     .      H0_L2, H1_L2, S0_L2, S1_L2, H0_R2, H1_R2, S0_R2, S1_R2)

        if (ABS(Ef_LeadL-Ef_LeadR).GT.1.0d-5) then
          if (iv.eq.0.and.Mynode.eq.0.and.ITER.EQ.1.and.ik.eq.1)
     .         write(6,'(a63)')
     .        'smeagol: different Fermi energies for left and right elec
     .trodes'
          different_eflr=.true.
          nspinMin=NspinBlocks
          if(NspinBlocks>1)nspinMin=2
          DO ISPIN=1,nspinMin
            H0_R2(:,:,ISPIN)=H0_R2(:,:,ISPIN)+(Ef_Lead-Ef_LeadR)*S0_R2
            H1_R2(:,:,ISPIN)=H1_R2(:,:,ISPIN)+(Ef_Lead-Ef_LeadR)*S1_R2
          ENDDO
        else
          different_eflr=.false.
        endif
      endif

#ifdef MPI
! If parallel broacasts the information about the electrods to all the other
! compute nodes
      call broadcast_leads(NspinBlocks,NL,H0_L2,S0_L2,H1_L2,S1_L2,
     &                     NR,H0_R2,S0_R2,H1_R2,S1_R2)
#endif

      call convertcopy_hsleads_sigma(NspinBlocks,NspinComplexMatrix,
     &    NL,NR)

      nl2=nl
      nr2=nr
      n12=n1

      if(NspinBlocks>2)then
        nl2=2*nl
        nr2=2*nr
        n12=2*n1
      endif


      if((N1-NL-NR).le.0)negfon=.false.

! Sets the values of the Hamiltonian and Overlap matrices at the boundaries 
      if(mynode_inverse.eq.0)then
        sgeneralp%mattype=2
        do ispin=1,NspinComplexMatrix
          hgeneralp(ispin)%mattype=2
        enddo

      if(.false.)then
        call PrintSparse2DenseReorderedNC(hgeneralp,n1,
     .      NspinComplexMatrix, NspinBlocks, NspinComplexMatrix,1,
     .      "hn_1")
        call PrintSparse2DenseReorderedNC(sgeneralp,n1, 1, NspinBlocks,
     .      NspinComplexMatrix,2,"sn_1")
      endif



        call set_boundaryelementsgeneral(n_replace_L,n_replace_R,
     .      NspinBlocks, NspinComplexMatrix,
     &        V,N1,NL,H0_L2,S0_L2,H1_L2,S1_L2,
     &        NR,H0_R2,S0_R2,H1_R2,S1_R2,hgeneralp,sgeneralp,
     &        Set_HBoundary_Leads,Set_HLR_Zero,HStol)

      if(.false.)then
        call PrintSparse2DenseReorderedNC(hgeneralp,n1,
     .      NspinComplexMatrix, NspinBlocks, NspinComplexMatrix,1,
     .      "hn_2")
        call PrintSparse2DenseReorderedNC(sgeneralp,n1, 1, NspinBlocks,
     .      NspinComplexMatrix,2,"sn_2")
      endif


      endif

      if(nprocs_hs.ne.1)then
        sgeneralp%mattype=3
        do ispin=1,NspinComplexMatrix
          hgeneralp(ispin)%mattype=3
        enddo
        call set_boundaryelementsgeneral(n_replace_L,n_replace_R,
     .      NspinBlocks, NspinComplexMatrix,
     &        V,N1,NL,H0_L2,S0_L2,H1_L2,S1_L2,
     &        NR,H0_R2,S0_R2,H1_R2,S1_R2,hgeneralp,sgeneralp,
     &        Set_HBoundary_Leads,Set_HLR_Zero)
        sgeneralp%mattype=2
        do ispin=1,NspinComplexMatrix
          hgeneralp(ispin)%mattype=2
        enddo
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
          write(12347,'(A,f12.6)')
     $         't_negf2',(sc_1b-sc_0b)*1.0d0/sc_rb
          CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
        endif



        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

        call MaxDifferenceMatrixCRS_CRSP(sgeneralp,sgeneralp,
     .      mynode_inverse,"(sgeneralap) ",io)

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_matrixdifference_s',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

        do ispin=1,NspinComplexMatrix
          call MaxDifferenceMatrixCRS_CRSP(hgeneralp(ispin),
     .        hgeneralp(ispin),mynode_inverse,"(hgeneraludap) ",io)

          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           't_matrixdifference_h',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif

        enddo
      endif
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
        write(12347,'(A,f12.6)')
     $       't_negf3',(sc_1b-sc_0b)*1.0d0/sc_rb
        CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
      endif

      if(NspinComplexMatrix>2)then
        if (sco.and.(converged.or.scoSCF)) then
          SCOSetHamiltonianBlock=.true.
          call SCOApplyK_nc(hgeneralp,sgeneralp,n1,NspinComplexMatrix,
     .        nl,nr)
        else
          SCOSetHamiltonianBlock=.false.
        endif
      endif



!calcualte the energy grid on the real axis and the self-energies for each
!energy point.

      allocate(leadsdim(nleadslr),LeadsVoltageShift(nleadslr))

      leadsdim(1)=nl
      leadsdim(2)=nr
      LeadsVoltageShift(1)=0.5D0 * V
      LeadsVoltageShift(2)=-0.5D0 * V

      call energygrid_selfenergies_real2(slabel,ef_lead,v,kb*T,iter,
     .    myhead,nheads,n1,NspinBlocks,delta,ik,hgeneralp,sgeneralp,
     .    istep,inicoor,vini,vfinal,nenerg,nk,
     .    nleadslr,leadsdim, storesigma, LeadsVoltageShift,sigmatodisk)

!#ifdef MPI
!      call MPI_Finalize( MPIerror )
!#endif
!      stop



      deallocate(leadsdim,LeadsVoltageShift)
      
!      write(*,*)"exiting ",mynode
!      call MPI_Finalize( MPIerror )
!      stop
!
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
        write(12347,'(A,f12.6)')
     $       't_negf4',(sc_1b-sc_0b)*1.0d0/sc_rb
        CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
      endif



      IF (.NOT. CONVERGED) THEN

        if(emtimings)CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)

        if((.not.ldos))then

          allocate(leadsdim(nleadslr),LeadsVoltageShift(nleadslr))

          leadsdim(1)=nl
          leadsdim(2)=nr
          LeadsVoltageShift(1)=0.5D0 * V
          LeadsVoltageShift(2)=-0.5D0 * V

          call energygrid_selfenergies_imag2(slabel,iter,ik,
     &        nenerg1,nenerg2, npoles, istep, inicoor,
     &        NspinBlocks,V,T,Ef_Lead,r0,eb,
     &        nleadslr,leadsdim,nk,storesigma,LeadsVoltageShift,
     &        sigmatodisk)

          deallocate(leadsdim,LeadsVoltageShift)
        endif
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'seneib',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif


#ifdef MPI
        CALL MPI_BARRIER(negf_comm, MPIerror)
#endif

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'seneia',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

!CCCCCCC Allocating Memory CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        CALL TIMER('NEGFK_su',2)
! ########## Charge in equilibrium loop ###################

        if(.not.ldos.and..not.SkipEquilibriumRho)then

           ALLOCATE(rho_f_im_l(NspinComplexMatrix,
     &                        EImagGrid%nEnergies,
     &                        rhogeneralp(1)%matSparse%nnz))

           ALLOCATE(rho_f_im_r(NspinComplexMatrix,
     &                        EImagGrid%nEnergies,
     &                        rhogeneralp(1)%matSparse%nnz))

            ALLOCATE(rho_im_l(NspinComplexMatrix,
     &                        EImagGrid%nEnergies,
     &                        rhogeneralp(1)%matSparse%nnz))

           ALLOCATE(rho_im_r(NspinComplexMatrix,
     &                        EImagGrid%nEnergies,
     &                        rhogeneralp(1)%matSparse%nnz))

          rho_im_l=0D0
          rho_im_r=0D0
          rho_f_im_l=0D0
          rho_f_im_r=0D0

          CALL TIMER('NEGFK_imag',1)
          call keldyshImag(N12,NspinBlocks,NspinComplexMatrix,NL2,NR2,
     .        ik,iter, NEnerg1,NEnerg2,NPOLES,R0,
     .        V,T,Ef_Lead,
     .        hgeneralp,sgeneralp,rhogeneralp,ematgeneralp,weightrho,
     .        rho_im_l, rho_im_r, rho_f_im_l, rho_f_im_r)
          CALL TIMER('NEGFK_imag',2)

      ! equilibrium contribution to density evaluated on imaginary space contour
      if(mynode_inverse.eq.0)then
        DO ISPIN=1,EImagGrid%nSpin
         DO I=1,EImagGrid%nEnergies

              rhogeneralp_l_im(ispin)%matSparse%b(:)=
     &            rhogeneralp_l_im(ispin)%matSparse%b(:)+
     &            rho_im_l(ISPIN, I, :)

              rhogeneralp_r_im(ispin)%matSparse%b(:)=
     &            rhogeneralp_r_im(ispin)%matSparse%b(:)+
     &             rho_im_r(ISPIN, I, :)

               if(emforces)then

                ematgeneralp_l_im(ispin)%matSparse%b(:)=
     &            ematgeneralp_l_im(ispin)%matSparse%b(:)+
     &            rho_f_im_l(ISPIN, I, :)

                ematgeneralp_r_im(ispin)%matSparse%b(:)=
     &            ematgeneralp_r_im(ispin)%matSparse%b(:)+
     &            rho_f_im_r(ISPIN, I, :)

                 endif ! if(emforces)then

             enddo ! DO ISPIN=1,ERealGrid%nSpin
        enddo ! DO I=1,ERealGrid%nEnergies

          DEALLOCATE(rho_im_l)
          DEALLOCATE(rho_im_r)
          DEALLOCATE(rho_f_im_l)
          DEALLOCATE(rho_f_im_r)

        endif  !if(mynode_inverse.eq.0)then
       endif  ! if(.not.ldos.and..not.SkipEquilibriumRho)then

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
          write(12347,'(A,f12.6)')
     $         't_negf5',(sc_1b-sc_0b)*1.0d0/sc_rb
          CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
        endif



        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'keldyshib',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

#ifdef MPI
        CALL MPI_BARRIER(negf_comm, MPIerror)
#endif
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'keldyshia',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

        if(.not.SkipNonEquilibriumRho.and.v.ne.0.0D0)then
          CALL TIMER('NEGFK_real',1)

          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
            write(12347,'(A,f12.6)')
     $           't_negf5b',(sc_1b-sc_0b)*1.0d0/sc_rb
            CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
          endif

          ALLOCATE(rho_real_l(NspinComplexMatrix,
     &                        ERealGrid%nEnergies,
     &                        rhogeneralp(1)%matSparse%nnz))
          ALLOCATE(rho_real_r(NspinComplexMatrix,
     &                        ERealGrid%nEnergies,
     &                        rhogeneralp(1)%matSparse%nnz))

          rho_real_l=0D0
          rho_real_r=0D0

          call keldyshReal(N12,NspinBlocks,NspinComplexMatrix,
     &        NL2,NR2,nk,ik,
     &        V,T,Delta,Ef_Lead,rhobs_general,ematbs_general,
     &        hgeneralp,sgeneralp,rhogeneralp,ematgeneralp,
     &        iter,ldos,weightrho, rho_real_l, rho_real_r)

          ! non-equilibrium contribution to density evaluated on real space contour
      DO ISPIN=1,ERealGrid%nSpin
       DO I=1,ERealGrid%nEnergies

              rhogeneralp_l_real(ispin)%matSparse%b(:)=
     &            rhogeneralp_l_real(ispin)%matSparse%b(:)+
     &            rho_real_l(ISPIN, I, :)

              rhogeneralp_r_real(ispin)%matSparse%b(:)=
     &            rhogeneralp_r_real(ispin)%matSparse%b(:)+
     &             rho_real_r(ISPIN, I, :)

               if(emforces)then
                Ei=ERealGrid%e(i)

                ematgeneralp_l_real(ispin)%matSparse%b(:)=
     &            ematgeneralp_l_real(ispin)%matSparse%b(:)+
     &            Ei*rho_real_l(ISPIN, I, :)

                ematgeneralp_r_real(ispin)%matSparse%b(:)=
     &            ematgeneralp_r_real(ispin)%matSparse%b(:)+
     &            Ei*rho_real_r(ISPIN, I, :)

                 endif ! if(emforces)then

             enddo ! DO ISPIN=1,ERealGrid%nSpin
      enddo ! DO I=1,ERealGrid%nEnergies

          DEALLOCATE(rho_real_l)
          DEALLOCATE(rho_real_r)


          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
            write(12347,'(A,f12.6)')
     $           't_negf5c',(sc_1b-sc_0b)*1.0d0/sc_rb
            CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
          endif




        endif

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'keldyshrb',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

#ifdef MPI
        CALL MPI_BARRIER(negf_comm, MPIerror)
#endif

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'keldyshra',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
          write(12347,'(A,f12.6)')
     $         't_negf6',(sc_1b-sc_0b)*1.0d0/sc_rb
          CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
        endif



! ########## End of Energy Loop out of equilibrium ###################################


!C ************************************************************************************
!C Calculate the total charge density by summing equilibrium and out-of-equilibrium
!C contributions


        CALL TIMER('NEGFK_col',1)


        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'beforebarrier',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

#ifdef MPI
        CALL MPI_BARRIER(negf_comm, MPIerror)
#endif
        

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'afterbarrier',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
          write(12347,'(A,f12.6)')
     $         't_negf6b',(sc_1b-sc_0b)*1.0d0/sc_rb
          CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
        endif

! Collect real and imaginary components across all MPI ranks
      if(emforces)then
        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        ematgeneralp_l_real, ematgeneralp_l_real)
        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        ematgeneralp_r_real, ematgeneralp_r_real)
        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        ematgeneralp_l_im, ematgeneralp_l_im)
        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        ematgeneralp_r_im, ematgeneralp_r_im)
      endif ! if(emforces)then

        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        rhogeneralp_l_real, rhogeneralp_l_real)
        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        rhogeneralp_r_real, rhogeneralp_r_real)
        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        rhogeneralp_l_im, rhogeneralp_l_im)
        call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &        rhogeneralp_r_im, rhogeneralp_r_im)

! Sum density matrices real and imaginary components
       weight=weightrho
       do ispin=1,NspinComplexMatrix
           do ii=1,rhogeneralp(ispin)%matSparse%nnz

        if(.not.SkipNonEquilibriumRho.and.v.ne.0.0D0)then

! Calculate weight from non-equilibrium density
         weight=(rhogeneralp_l_real(ispin)%matSparse%b(ii))**2.0D0/
     &           (rhogeneralp_l_real(ispin)%matSparse%b(ii)**2.0D0+
     &            rhogeneralp_r_real(ispin)%matSparse%b(ii)**2.0D0)

        endif !if(.not.SkipNonEquilibriumRho.and.v.ne.0.0D0)then

         if(emforces)then
         ematgeneralp(ispin)%matSparse%b(ii) =
     &     ematgeneralp(ispin)%matSparse%b(ii) +
     &      (weight)*(ematgeneralp_l_im(ispin)%matSparse%b(ii) +
     &                ematgeneralp_r_real(ispin)%matSparse%b(ii)) +
     &      (1D0 - weight)*(ematgeneralp_r_im(ispin)%matSparse%b(ii) +
     &                       ematgeneralp_l_real(ispin)%matSparse%b(ii))
          endif

         rhogeneralp(ispin)%matSparse%b(ii) =
     &     rhogeneralp(ispin)%matSparse%b(ii) +
     &      (weight)*(rhogeneralp_l_im(ispin)%matSparse%b(ii) +
     &                rhogeneralp_r_real(ispin)%matSparse%b(ii)) +
     &      (1D0 - weight)*(rhogeneralp_r_im(ispin)%matSparse%b(ii) +
     &                       rhogeneralp_l_real(ispin)%matSparse%b(ii))

           enddo ! do ii=1,rhogeneralp(ispin)%matSparse%nnz
        enddo ! do ispin=1,NspinComplexMatrix

#ifdef MPI
         CALL MPI_BARRIER(negf_comm, MPIerror )
#endif

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'rhogreduce',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
          write(12347,'(A,f12.6)')
     $         't_negf7',(sc_1b-sc_0b)*1.0d0/sc_rb
          CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
        endif


        ! Destroy density matrices for real and imaginary components
      do ispin=1,NspinComplexMatrix
        call DestroyMatrixGeneral(rhogeneralp_l_real(ispin),"negfk",io)
        call DestroyMatrixGeneral(rhogeneralp_r_real(ispin),"negfk",io)
        call DestroyMatrixGeneral(rhogeneralp_l_im(ispin),"negfk",io)
        call DestroyMatrixGeneral(rhogeneralp_r_im(ispin),"negfk",io)

        if(emforces)then
        call DestroyMatrixGeneral(ematgeneralp_l_real(ispin),"negfk",io)
        call DestroyMatrixGeneral(ematgeneralp_r_real(ispin),"negfk",io)
        call DestroyMatrixGeneral(ematgeneralp_l_im(ispin),"negfk",io)
        call DestroyMatrixGeneral(ematgeneralp_r_im(ispin),"negfk",io)
        endif  ! if(emforces)then

      enddo !  ispin=1,NspinComplexMatrix


        CALL TIMER('NEGFK_col',2)

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'negfkcollect',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
          write(12347,'(A,f12.6)')
     $         't_negf8',(sc_1b-sc_0b)*1.0d0/sc_rb
          CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
        endif


      ELSE
! ######################################################################################
!       If the density matrix has converged, we calculate the current
!
!TO   Either of IETS or normal TrCoef

        IF (CALCIETS) THEN

          CALL TIMER('LOE',1)

          call loe_wrap(ik,nk,N12,nl2,nr2,NspinBlocks,V,IV,
     &        slabel,kpoint,wk, Ef_Lead,kb,t,
     &        hgeneralp,sgeneralp,rhogeneralp,ematgeneralp)
          CALL TIMER('LOE',2)

        ELSEIF (CALCTRANSM .and. 
     .     (idyn .lt. 1 .or. idyn .gt. 6 .or.    ! not MD or FC
     .      idyn .eq. 6 .and. istep .eq. 0 .or.  ! FC, only undisplaced coordinate
     .      idyn .gt. 0 .and. idyn .lt. 6 .and.  
     .      istep.gt.tmdskip .and. 
     .      mod(istep-tmdskip-1, tmdsampling).eq.0) ! MD
     .  ) THEN

          CALL TIMER('TRANSM',1)

          call transm_wrap(ik,nk,N12,nl2,nr2,NspinBlocks,
     &        NspinComplexMatrix,V,IV,
     &        slabel,kpoint,wk, Ef_Lead,kb,t,
     &        hgeneralp,sgeneralp,rhogeneralp,ematgeneralp,
     .        istep, inicoor, idyn) !!
          CALL TIMER('TRANSM',2)

        ENDIF

        if(NspinBlocks<=2)then
        CALL TIMER('CURRENT',1)
       
        If (ik .EQ. 1) then
          Ic=0.D0
          Icfinal=0.d0
        endif

        Icb=0D0
        if(.not.negfon.and.mynode_inverse.eq.0.and.NspinBlocks<=2)then
          DO I=1,ERealGrid%nEnergies
            Ei=ERealGrid%e(i)

            DO ISPIN=1,NspinBlocks

              IF ((ERealGrid%GridType.eq.1))then
                iksigma=1
              else
                iksigma=ik
              endif
              call current(I,nk,N1,NL,NR,NspinBlocks,ISPIN,
     $            iksigma,Ei,
     &            V,T,Ef_Lead,
     &            ERealGrid%w(I),Delta,
     $            wk,Icb,tk,
     &            hgeneralp,sgeneralp)

            ENDDO

            if(writetk)then
              if(NspinBlocks.eq.1)then
                write(12347,*)"tk=",ei,tk(1),ivmod,ik
              else
                write(12347,*)"tk=",ei,tk(1)+tk(2),tk(1),tk(2),ivmod,ik
              endif
            endif

          ENDDO
        endif

        DO ISPIN=1,NspinBlocks
          Ic(ISPIN)=Ic(ISPIN)+Icb(ISPIN)
        ENDDO
        If (ik .EQ. nk) then

#ifdef MPI
          CALL MPI_BARRIER(negfo_comm, MPIerror )
#ifdef NODAT
          CALL MPI_REDUCE(Ic(1:NspinBlocks),Icfinal,NspinBlocks,
     .        MPI_DOUBLE_PRECISION,
     .        MPI_SUM,0,negfo_comm,MPIerror)
#else
          CALL MPI_REDUCE(Ic(1:NspinBlocks),Icfinal,NspinBlocks,
     .        DAT_double, MPI_SUM,0,negfo_comm,MPIerror)
#endif
#else
          Icfinal=Ic
#endif
          If (mynode_negfo.EQ.0) Then
            call io_assign(iuc)
            curfile = paste(slabel,'.CUR')
!            If (IV .EQ. 0) Then
!              OPEN(UNIT=iuc,FILE=curfile)
!            Else
            OPEN(UNIT=iuc,FILE=curfile,POSITION='append')
!            EndIf
            If (NspinBlocks .EQ. 1) Icfinal=2.0D0*Icfinal
            IF (NspinBlocks .EQ. 1) THEN
              ! Meilin Bai, Dec 2012
              if (idyn .gt. 0 .and. idyn .le. 6) then  ! MD or FC
                WRITE(iuc,'(2d17.8,i7)') 13.6057D0 * V,Icfinal, istep 
              else
              ! end Meilin Bai
                WRITE(iuc,'(2d17.8)') 13.6057D0 * V,Icfinal
              endif

            ELSE
              ! Meilin Bai, Dec 2012
              if (idyn .gt. 0 .and. idyn .lt. 6) then
                WRITE(iuc,'(4d17.8,i7)')
     .               13.6057D0 * V,SUM(Icfinal),Icfinal,istep
              else
              ! end Meilin Bai
                WRITE(iuc,'(4d17.8)') 13.6057D0 * V,SUM(Icfinal),Icfinal
              endif

            ENDIF
            call io_close(iuc)
          EndIf
        EndIf

        CALL TIMER('CURRENT',2)
        endif

!        IF (ALLOCATED(EImagGrid%e).and.ik.eq.nk) THEN
!          call deallocate_energygrid_selfenergies(EImagGrid)
!        ENDIF


        if(emldos2.or.curr_dist.or.TransmissionChannels)then
          if(mynode_inverse.eq.0)then
            call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &          rhobs_general, rhogeneralp)
            if(emforces)then
              call rhoG_reduce(N1,NL,NR,NspinBlocks,NspinComplexMatrix,
     &            ematbs_general, ematgeneralp)
            endif

!             call WriteMatrixSparse(rhogeneralp(1),"rhot1")
!             call WriteMatrixSparse(rhogeneralp(2),"rhot2")

            if(emSTTk.and.mynode.eq.0)then
              allocate(Rhopdostotk(n1,NspinBlocks))
              call GetRhoPDOS_nc_noON(rhogeneralp,sgeneralp,NspinBlocks,
     .            .true.,emdostotk,Rhopdostotk)
              Rhopdostotk=Rhopdostotk*emSTTkScale
              call writePDOS(iv,ikpmod,'.TRC.RHOPDOS',slabel,
     .            NspinBlocks,n1,1,0.0D0,0.0D0,Rhopdostotk)
              deallocate(Rhopdostotk)
            endif


          endif
        endif

      ENDIF

!      IF ((ERealGrid%GridType.eq.1).or.
!     .    (allocated(ERealGrid%e).and.ik.eq.nk.and.
!     *    ((converged.and.(.not.emldos)).or.ldos))) then
!        call deallocate_energygrid_selfenergies(ERealGrid)
!      ENDIF

      CALL TIMER('NEGFK',2)      


      if(curr_distK)then
        call CurrentDistributionMatrix_RhoK(rhogeneralp,ematgeneralp,
     .      hgeneralp,sgeneralp,NspinComplexMatrix,      
     .      n1, v,wk,ik)
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
        write(12347,'(A,f12.6)')
     $       't_negf9',(sc_1b-sc_0b)*1.0d0/sc_rb
        CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1a,sc_rb,sc_mb)
        write(12347,'(A,f12.6)')
     $       't_negf_tot',(sc_1a-sc_0a)*1.0d0/sc_rb
      endif

      END SUBROUTINE NEGFK


