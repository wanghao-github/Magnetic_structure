* Sample Fortran source code for reading the data file, PIR_SSG_data.txt
* This data file contains the 2011 version of the matrices for the
* physically irreducible representations of the 230 crystallographic space
* groups for irrational k vectors in (3+d)-dimensional superspace.
*****************************************************************************
* Routines
*  pir_ssg_data_read: read data from file, PIR_SSG_data.txt
*  pir_ssg_data_unpack_operatormatrix: get matrix for given operator in space
*    group
*  pir_ssg_data_unpack_irmatrix_pointoperator: get point operator part of
*    IR matrix for given operator
*  pir_ssg_data_unpack_irmatrix_translation: get translation part of
*    IR matrix for given translation
*  pir_ssg_data_get_irmatrix:  get IR matrix for given space group operator
*  pir_ssg_data_get_transform2complex: get matrix that transforms PIR to
*    complex block-diagonal form
* Routines for testing data
*  pir_ssg_data_write: write data into file, PIR_SSG_data2.txt.  Compare
*    PIR_SSG_data2.txt with PIR_SSG_data.txt to determine if the data file was
*    read correctly
*  pir_ssg_data_test_multiplication_table: check that the IR matrices have the
*    same multiplication table as the space group operators
* Additional routines used by the above routines
*  pir_ssg_data_constant: store 16-byte double precision numbers as
*    1-byte integers
*  pir_ssg_data_vadd: add two vectors with rational components
*  pir_ssg_data_vsub: subtract two vectors with rational components
*  pir_ssg_data_vmlt: multiply a vector of rational numbers by a matrix of
*    integers
*  pir_ssg_data_opmlt: multiply two operator matrices
*  pir_ssg_data_dmatmlt: multiply two matrices of double precision numbers
*  pir_ssg_data_factor: remove greatest common factor from a set of integers
*****************************************************************************
      module pir_ssg_data_module
* storage area for data from table
* number of IRs
      integer ircount
      parameter (ircount=10294)
* for the ith IR:
* spacegroupnumber(i), space group number (1-230)
      integer spacegroupnumber(ircount)
* spacegroupsymbol(i), space group symbol
      character spacegroupsymbol(ircount)*10
* irsymbol(i), IR symbol (Miller-Love convention)
      character irsymbol(ircount)*8
* irdimension(i), dimension of IR
      integer irdimension(ircount)
* irtype(i), type of IR (1,2,3)
      integer irtype(ircount)
* kcount(i), number of k vectors in the star of k
      integer kcount(ircount)
* pmkcount(i), number of k vectors in the star of +/-k
      integer pmkcount(ircount)
* operatorcount(i), number of representative operators in space group
      integer operatorcount(ircount)
* kvectorpointer(i), pointer to k vectors in kvector
      integer kvectorpointer(ircount)
* kvector(m), k vectors stored sequentially beginning at
*   kvectorpointer(i).  Each k vector is stored as kvec(4,4), where
*   kvec(1:3,1) are the x,y,z components of the constant term,
*   kvec(1:3,2) are the x,y,z components of the alpha term,
*   kvec(1:3,3) are the x,y,z components of the beta term,
*   kvec(1:3,4) are the x,y,z components of the gamma term,
*   kvec(4,m) are the common denominators.
*   alpha,beta,gamma are free parameters for non-special k vectors
      integer*1 kvector(240000)

* for the jth space group operator:
* operatormatrixpointer(j,i), pointer to operator matrix in
*   operatormatrix
      integer operatormatrixpointer(48,ircount)
* operatormatrix(m), operator matrix stored sequentially beginning
*   at operatormatrixpointer(j,i).  This matrix is augmented with
*   the point operation in A(1:3,1:3), the translation in A(1:3,4),
*   and a common denominator in A(4,4).
      integer*1 operatormatrix(4800000)
* irmatrixpointer(j,i), pointer to IR matrix in irmatrix
      integer irmatrixpointer(48,ircount)
* irmatrix(m), IR matrix stored sequentially beginning at
*   irmatrixpointer(j,i).  The elements of this matrix are stored
*   as integers which are mapped onto double precision values using
*   the routine, pir_ssg_data_constant.
      integer*1 irmatrix(8400000)
      end module pir_ssg_data_module

****************************************************************************
* program to show how to use the routines
      program pir_ssg_data

      use pir_ssg_data_module

      implicit none
      integer i

* read data base from PIR_SSG_data.txt
      call pir_ssg_data_read

* write data base to PIR_SSG_data2.txt.  If you then compare it with 
* PIR_SSG_data.txt, you can determine whether the data was read correctly.
*      call pir_ssg_data_write

* test IR matrices: Do they have the same multiplcation table as the operators?
*      do i=1,ircount
*        if(spacegroupnumber(i).eq.0)cycle
*        write(6,*)i
*        call pir_ssg_data_test_multiplication_table(i)
*      enddo

      end

**************************************************************************

      subroutine pir_ssg_data_read
* Read data file, PIR_SSG_data.txt

      use pir_ssg_data_module

      implicit none
      integer i,j,k,m,n4,kp,op,ip,irnumber
      double precision irtemp(2304)

* initialize data
      spacegroupnumber=0
      spacegroupsymbol=' '
      irsymbol=' '
      irdimension=0
      pmkcount=0
      operatorcount=0
      kvectorpointer=0
      kvector=0
      operatormatrix=0
      irmatrix=0
* initialize pointers
      kp=0
      op=0
      ip=0
* open data file
      open(30,file='PIR_SSG_data.txt')
* skip title line
      read(30,*)
      read(30,*)
      read(30,*)
* read next IR
      do while(.true.)
* read IR number, space group number, space group symbol, IR label,
* IR dimension, number of k vectors in the star of +/-k,
* number of operators in space group
        read(30,*,end=1)irnumber,spacegroupnumber(irnumber),
     $       spacegroupsymbol(irnumber),irsymbol(irnumber),
     $       irdimension(irnumber),irtype(irnumber),kcount(irnumber),
     $       pmkcount(irnumber),operatorcount(irnumber)
* read each k vector in the star of +/-k
        kvectorpointer(irnumber)=kp+1
        read(30,*)kvector(kp+1:kp+16*pmkcount(irnumber))
        kp=kp+16*pmkcount(irnumber)
* dimension of operator in extended matrix form
        n4=pmkcount(irnumber)+4
* do each operator in space group
        do i=1,operatorcount(irnumber)
* read operator in extended matrix form
          operatormatrixpointer(i,irnumber)=op+1
          read(30,*)operatormatrix(op+1:op+n4**2)
          op=op+n4**2
* read IR matrix
          irmatrixpointer(i,irnumber)=ip+1
          read(30,*)irtemp(1:irdimension(irnumber)**2)
          do j=1,irdimension(irnumber)**2
            call pir_ssg_data_constant(2,irmatrix(ip+j),irtemp(j))
          enddo
          ip=ip+irdimension(irnumber)**2
* next operator
        enddo
* next IR
      enddo
 1    close(30)
      end
**************************************************************************
      subroutine pir_ssg_data_constant(choice,n,x)
* store 16-byte double precision numbers as 1-byte integers
* arguments:
*   choice (input): choice=1, input n, output x
*                   choice=2, input x, output n
*   n (input/output): 1-byte integer
*   x (input/output): 16-byte double precision
      implicit none
      integer choice,i
      integer*1 n
      double precision x,constants(17)
      data constants/0,1,-1,0.5,-0.5,0.25,-0.25,
     $     0.866025403784439,-0.866025403784439,  ! sqrt(3)/2
     $     0.707106781186548,-0.707106781186548,  ! sqrt(2)/2
     $     0.433012701892219,-0.433012701892219,  ! sqrt(3)/4
     $     0.683012701892219,-0.683012701892219,  ! cos(15)/sqrt(2)
     $     0.183012701892219,-0.183012701892219/  ! sin(15)/sqrt(2)
* choice=1: input integer and output double precision
      if(choice.eq.1)then
        x=constants(n)
* choice=2: input double precision, output integer
      else
        do i=1,17
          if(dabs(x-constants(i)).lt.1d-4)then
            n=i
            return
          endif
        enddo
        stop 'Error in pir_ssg_data_constant: value of x not found.'
      endif
      end
*****************************************************************************
      subroutine pir_ssg_data_write
* check whether the data base was read correctly bu writing the data to
* PIR_SSG_data2.txt and then comparing it with PIR_SSG_data.txt.

      use pir_ssg_data_module

      implicit none
      integer i,j,k,m,n,n4,kp,op,ip,i1,i2,nd,irnumber
      double precision xmat(48,48),x
      character lineout*10000

* open file
      open(20,file='PIR_SSG_data2.txt')
* title line
      write(20,'(a)')'ISO-IR: Physically Irreducible Representations '
     $     //'of the 230 Crystallographic Space Groups'
      write(20,'(a)')'2011 version of matrices, Operators in '
     $     //'(3+d)-dimensional superspace'
      write(20,'(a)')'Harold T. Stokes and Branton J. Campbell, 2013'
* do each IR
      do irnumber=1,ircount
        if(spacegroupnumber(irnumber).eq.0)cycle
* write IR number, space group number, space group symbol, IR label,
* IR dimension, number of k vectors in the star of +/-k,
* number of operators in space group
        write(20,'(i5,i4,a,5i3)')irnumber,spacegroupnumber(irnumber),
     $       ' "'//spacegroupsymbol(irnumber)//'" "'
     $       //irsymbol(irnumber)//'"',
     $       irdimension(irnumber),irtype(irnumber),kcount(irnumber),
     $       pmkcount(irnumber),operatorcount(irnumber)
* write each k vector in the star of +/-k
        kp=kvectorpointer(irnumber)-1
        write(20,'(24i3)')kvector(kp+1:kp+16*pmkcount(irnumber))
* dimension of operator in extended matrix form
        n4=pmkcount(irnumber)+4
* do each operator in space group
        do i=1,operatorcount(irnumber)
* write operator in extended matrix form
          op=operatormatrixpointer(i,irnumber)-1
          write(20,'(20i4)')operatormatrix(op+1:op+n4**2)
* write IR matrix
          ip=irmatrixpointer(i,irnumber)-1
          n=0
          do j=1,irdimension(irnumber)
          do k=1,irdimension(irnumber)
            n=n+1
            call pir_ssg_data_constant(1,irmatrix(ip+n),xmat(j,k))
          enddo
          enddo
        lineout=' '
        m=0
        nd=irdimension(irnumber)
        do j=1,nd
        do k=1,nd
          x=xmat(j,k)
          if(dabs(x).lt.1d-4)then
            lineout(m+1:m+2)=' 0'
            m=m+2
          else if(dabs(x-1).lt.1d-4)then
            lineout(m+1:m+2)=' 1'
            m=m+2
          else if(dabs(x+1).lt.1d-4)then
            lineout(m+1:m+3)=' -1'
            m=m+3
          else
            write(lineout(m+1:m+10),'(f10.5)')x
            m=m+1
            do n=m+1,m+9
              if(lineout(n:n).ne.' ')then
                lineout(m+1:m+1)=lineout(n:n)
                m=m+1
              endif
            enddo
            lineout(m+1:m+10)=' '
          endif
        enddo
        enddo
        i1=1
        i2=i1+79
        if(i2.gt.m)i2=m
        do while(i1.lt.m)
          do j=i2+1,i1,-1
            if(lineout(j:j).eq.' ')then
              write(20,'(a)')lineout(i1:j-1)
              i1=j
              i2=i1+79
              if(i2.gt.m)i2=m
              exit
            endif
          enddo
        enddo
* next operator
        enddo
* next IR
      enddo
 1    close(20)
      end
*****************************************************************************
      subroutine pir_ssg_data_unpack_operatormatrix(irnumber,opnumber,
     $     operatormatrix_out,nr)
* get operators and IR matrices for a selected IR
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   opnumber (input): selected operator number
*   operatormatrix_out(k,j) (output): operator matrix
*   nr (input): number of rows in the array operatormatrix_out

      use pir_ssg_data_module

      implicit none
      integer i,j,n,n4,nr,op,irnumber,opnumber,operatormatrix_out(nr,nr)

* check for valid selection
      if(irnumber.lt.1.or.irnumber.gt.ircount)then
        write(6,'(a)')'Error in pir_ssg_data_unpack_operatormatrix: '
     $       //'invalid value of irnumber'
        stop
      endif
      if(opnumber.lt.1.or.opnumber.gt.operatorcount(irnumber))then
        write(6,'(a)')'Error in pir_ssg_data_unpack_operatormatrix: '
     $       //'invalid value of opnumber'
        stop
      endif
* pointer to matrix
      op=operatormatrixpointer(opnumber,irnumber)
* unpack matrix
      n4=pmkcount(irnumber)+4
      do i=1,n4
      do j=1,n4
        operatormatrix_out(i,j)=operatormatrix(op)
        op=op+1
      enddo
      enddo
      end
*****************************************************************************
      subroutine pir_ssg_data_unpack_irmatrix_pointoperator(irnumber,
     $     opnumber,irmatrix_out,nr)
* get point operator part of IR matrix for a selected IR and operator
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   opnumber (input): selected operator number
*   irmatrix_out(k,j) (output): IR matrix
*   nr (input): number of rows in the array, irmatrix_out

      use pir_ssg_data_module

      implicit none
      integer nr
      integer i,j,n,ip,irnumber,opnumber
      double precision irmatrix_out(nr,nr)

* check for valid selection
      if(irnumber.lt.1.or.irnumber.gt.ircount)then
        write(6,'(a)')'Error in irr_data_unpack_irmatrix_pointoperator:'
     $       //' invalid value of irnumber'
        stop
      endif
      if(opnumber.lt.1.or.opnumber.gt.operatorcount(irnumber))then
        write(6,'(a)')'Error in irr_data_unpack_irmatrix_pointoperator:'
     $       //' invalid value of opnumber'
        stop
      endif
      if(nr.lt.irdimension(irnumber))then
        write(6,'(a)')'Error in irr_data_unpack_irmatrix_pointoperator:'
     $       //' invalid value of nr'
        stop
      endif
* pointer to matrix
      ip=irmatrixpointer(opnumber,irnumber)
* unpack matrix
      do i=1,irdimension(irnumber)
      do j=1,irdimension(irnumber)
        call pir_ssg_data_constant(1,irmatrix(ip),irmatrix_out(i,j))
        ip=ip+1
      enddo
      enddo
      end
*****************************************************************************
      subroutine pir_ssg_data_get_irmatrix_phase(irnumber,
     $     phasevector,phasevector_denom,irmatrix_out,nr)
* get translation part of IR matrix for a selected IR and translation
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   phasevector(i) (input), phase shift of ith modulation
*   phasevector_denom (input), common denominator in phase shifts
*   irmatrix_out(k,j) (output): IR matrix
*   nr (input): number of rows in the array, irmatrix_out

      use pir_ssg_data_module

      implicit none
      integer i,j,k,n,ik,nr,irnumber,nd,nk,nb,kp,kvec(4,4),
     $     phasevector(24),phasevector_denom
      double precision irmatrix_out(nr,nr),cosval,sinval
      double precision, parameter :: pi = 3.1415926535897932384d0

* check for valid selection
      if(irnumber.lt.1.or.irnumber.gt.ircount)then
        write(6,'(a)')'Error in pir_ssg_data_get_irmatrix_phase: '
     $       //'invalid value of irnumber'
        stop
      endif
      if(phasevector_denom.eq.0)then
        write(6,'(a)')'Error in pir_ssg_data_get_irmatrix_phase: '
     $       //'invalid value of phasevector_denom'
        stop
      endif
      if(nr.lt.irdimension(irnumber))then
        write(6,'(a)')'Error in pir_ssg_data_get_irmatrix_phase: '
     $       //'invalid value of nr'
        stop
      endif

      irmatrix_out=0
* dimension of IR
      nd=irdimension(irnumber)
* number of k vectors in star of +-k
      nk=pmkcount(irnumber)
* dimension of block matrix in IR
      nb=nd/nk
* do each block
      do i=1,nk
* cos and sin values
        cosval=dcos(2*pi*phasevector(i)/phasevector_denom)
        sinval=dsin(2*pi*phasevector(i)/phasevector_denom)
* construct block
        n=(i-1)*nb
* cos value on diagonal
        do j=1,nb
          irmatrix_out(n+j,n+j)=cosval
        enddo
* sin value off diagonal
        if(dabs(sinval).gt.1d-6)then
          do j=1,nb/2
            irmatrix_out(n+j,n+nb/2+j)=sinval
            irmatrix_out(n+nb/2+j,n+j)=-sinval
          enddo
        endif
      enddo
      end
*****************************************************************************
      subroutine pir_ssg_data_get_irmatrix(irnumber,
     $     operatormatrix_in,irmatrix_out,nr1,nr2)
* get IR matrix for a selected IR and operator
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   operatormatrix_in(j,i) (input), matrix of selected operator
*   irmatrix_out(k,j) (output): IR matrix
*   nr1 (input): number of rows in the array, operatormatrix_in
*   nr2 (input): number of rows in the array, irmatrix_out

      use pir_ssg_data_module

      implicit none
      integer nrop,nrir
      integer i,j,k,n4,nd,nr1,nr2,nmod,
     $     irnumber,operatormatrix_in(nr1,nr1),
     $     operatormatrix_stored(28,28),pvector(25),
     $     opnumber
      double precision irmatrix_out(nr2,nr2),
     $     irmatrix_pointoperation(48,48),irmatrix_phase(48,48)

* number of modulations
      nmod=pmkcount(irnumber)
* dimension of operator in extended matrix form
      n4=nmod+4
* identify operator
* try every operator in data base
      iloop: do i=1,operatorcount(irnumber)
        call pir_ssg_data_unpack_operatormatrix(irnumber,i,
     $     operatormatrix_stored,28)
* compare point operator parts
        do j=1,n4-1
        do k=1,n4-1
* not the same
          if(operatormatrix_in(k,j)/operatormatrix_in(n4,n4)
     $         .ne.operatormatrix_stored(k,j)
     $         /operatormatrix_stored(n4,n4))then
* last operator: operator not found
            if(i.eq.operatorcount(irnumber))then
              write(6,'(a)')'Error in pir_ssg_data_get_irmatrix: '
     $             //'operator not found'
              stop
* try next operator in data base
            else
              cycle iloop
            endif
          endif
        enddo
        enddo
* found it
        opnumber=i
        exit
      enddo iloop
* get phase vector
      pvector(1:nmod+1)=operatormatrix_in(4:n4,n4)
      call pir_ssg_data_factor(nmod+1,pvector)
* get phase part of IR
      call pir_ssg_data_get_irmatrix_phase(irnumber,
     $     pvector,pvector(nmod+1),irmatrix_phase,48)
* get point operation part of IR
      call pir_ssg_data_unpack_irmatrix_pointoperator(irnumber,
     $     opnumber,irmatrix_pointoperation,48)
* multiply them to obtain final IR matrix
      irmatrix_out=0
      nd=irdimension(irnumber)
      do i=1,nd
      do j=1,nd
      do k=1,nd
        irmatrix_out(i,j)=irmatrix_out(i,j)+irmatrix_phase(i,k)
     $       *irmatrix_pointoperation(k,j)
      enddo
      enddo
      enddo
      end
*****************************************************************************
      subroutine pir_ssg_data_test_multiplication_table(irnumber)
* get the IR matrices and check that they have the same multiplication
* table as the operators
      use pir_ssg_data_module
      implicit none
      integer i,j,k,m,n,i1,i2,j1,j2,n4,
     $     irnumber,opmat(28,28,48),ndim,opmat2(28,28)
      double precision irmat(48,48,48),irmat2(48,48),irmat3(48,48)
* no data for this irrep
      if(spacegroupnumber(irnumber).eq.0)return
* get operator matrix for each operator
      opmat=0
      do i=1,operatorcount(irnumber)
        call pir_ssg_data_unpack_operatormatrix(irnumber,i,
     $       opmat(1,1,i),28)
      enddo
      n=operatorcount(irnumber)
* get IR matrices
      do i=1,n
        call pir_ssg_data_get_irmatrix(irnumber,
     $       opmat(1,1,i),irmat(1,1,i),28,48)
      enddo
* dimension of IR
      ndim=irdimension(irnumber)
* dimension of operator in extended matrix form
      n4=pmkcount(irnumber)+4
* do each pair of operators
      do i1=1,n
      do i2=1,n
* multiply operators
        call pir_ssg_data_opmlt(opmat(1,1,i1),
     $       opmat(1,1,i2),opmat2,n4,28,28,28)
* get IR for operator
        call pir_ssg_data_get_irmatrix(irnumber,
     $       opmat2,irmat2,28,48)
* multiply IR matrices
        call pir_ssg_data_dmatmlt(irmat(1,1,i1),irmat(1,1,i2),
     $       irmat3,ndim,ndim,ndim,48,48,48)
* compare IR matrices
        do j1=1,ndim
        do j2=1,ndim
          if(dabs(irmat3(j2,j1)-irmat2(j2,j1)).gt.1d-6)then
            write(6,*)'Error in pir_ssg_data_test_multiplication_table:'
     $           //' IR matrices do not have the same multiplication '
     $           //'table as the operators'
            stop
          endif
        enddo
        enddo
      enddo
      enddo
      end
***************************************************************************
      subroutine pir_ssg_data_get_transform2complex(irnumber,smat,
     $     smatinv,nr)
* construct matrix that transforms PIR to complex block-diagonal form
* arguments:
*   irnumber (input): selected IR number (1-10294)
*   smat(i,j) (output): transformation matrix
*   smatinv(i,j) (output): inverse of smat
*   nr (input): number of rows in the arrays, smat and smatinv

      use pir_ssg_data_module

      implicit none
      integer i,j,k,m,n,nr,irnumber,nmod,nd,nd2,n1,n2,
     $     ipermute(4,4),ndk
      logical minuskvec
      double precision sqrt2
      complex*16 smat(nr,nr),zmat(nr,nr),c2real(2,2),smatinv(nr,nr)
      data ipermute/1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0/
      data c2real/1,(0,1),1,(0,-1)/

      sqrt2=dsqrt(2d0)
* determine if -k is in the star of k
      if(pmkcount(irnumber).eq.kcount(irnumber))then
        minuskvec=.false.
      else
        minuskvec=.true.
      endif
* dimension of irrep matrix
      nd=irdimension(irnumber)
* number of blocks
      nmod=pmkcount(irnumber)
* size of block
      nd2=nd/nmod
* matrix to undo permutation if -k is in the star of k
      if(minuskvec.and.irtype(irnumber).ne.1)then
        smat=0
        do i=1,nmod
          do j=1,4
          do k=1,4
            n1=(j-1)*nd2/4+(i-1)*nd2
            n2=(k-1)*nd2/4+(i-1)*nd2
            do m=1,nd2/4
              smat(n1+m,n2+m)=ipermute(j,k)
            enddo
          enddo
          enddo
        enddo
      else
        smat=0
        do i=1,nd
          smat(i,i)=1
        enddo
      endif
* matrix to undo transformation to real form
      if(irtype(irnumber).ne.1.or.minuskvec)then
      zmat=0
      do i=1,nmod
        do j=1,2
        do k=1,2
          n1=(j-1)*nd2/2+(i-1)*nd2
          n2=(k-1)*nd2/2+(i-1)*nd2
          do m=1,nd2/2
            zmat(n1+m,n2+m)=dconjg(c2real(k,j))/sqrt2
          enddo
        enddo
        enddo
      enddo
      call pir_ssg_data_zmatmlt(zmat,smat,smat,nd,nd,nd,nr,nr,nr)
      endif
* matrix to permute blocks
      if(irtype(irnumber).ne.1)then
      zmat=0
      j=1
      do i=1,nmod*2
        n1=(i-1)*nd2/2
        n2=(j-1)*nd2/2
        do k=1,nd2/2
          zmat(n1+k,n2+k)=1
        enddo
        j=j+2
        if(j.gt.nmod*2)j=2
      enddo
      call pir_ssg_data_zmatmlt(zmat,smat,smat,nd,nd,nd,nr,nr,nr)
      endif
* inverse
      do i=1,nd
      do j=1,nd
        smatinv(j,i)=dconjg(smat(i,j))
      enddo
      enddo
      end
***************************************************************************
      subroutine pir_ssg_data_vadd(nv1,nv2,nv3)
* add two vectors of rational numbers: nv3=nv1+nv2
      implicit none
      integer nv1(4),nv2(4),nv3(4),j
      do j=1,3
        nv3(j)=nv1(j)*nv2(4)+nv2(j)*nv1(4)
      enddo
      nv3(4)=nv1(4)*nv2(4)
      call pir_ssg_data_factor(4,nv3)
      return
      end
****************************************************************************
      subroutine pir_ssg_data_vsub(nv1,nv2,nv3)
* subtract two vectors of rational numbers: nv3=nv1-nv2
      implicit none
      integer nv1(4),nv2(4),nv3(4),nv(4),j
      do j=1,3
        nv(j)=-nv2(j)
      enddo
      nv(4)=nv2(4)
      call pir_ssg_data_vadd(nv1,nv,nv3)
      end
*****************************************************************************
      subroutine pir_ssg_data_vmlt(mat,nv1,nv2)
* multiply a 3x3 matrix times a vector: nv2=mat*nv1
      implicit none
      integer mat(3,3),nv1(4),nv2(4),nv3(4),j,k
      nv3=0
      do j=1,3
      do k=1,3
        nv3(j)=nv3(j)+mat(j,k)*nv1(k)
      enddo
      enddo
      nv3(4)=nv1(4)
      call pir_ssg_data_factor(4,nv3)
      nv2=nv3
      end
*****************************************************************************
      subroutine pir_ssg_data_opmlt(op1,op2,op3,n4,nr1,nr2,nr3)
* multiply two operator matrices op3=op1*op2
      implicit none
      integer i,j,k,m,n,n4,nr1,nr2,nr3,op1(nr1,nr1),op2(nr2,nr2),
     $     op3(nr3,nr3),op(n4,n4)
      op=0
      do i=1,n4
      do j=1,n4
      do k=1,n4
        op(i,j)=op(i,j)+op1(i,k)*op2(k,j)
      enddo
      enddo
      enddo
* remove common factor
      call pir_ssg_data_factor(n4**2,op)
* copy into output
      op3(1:n4,1:n4)=op(1:n4,1:n4)
      end
*****************************************************************************
      subroutine pir_ssg_data_dmatmlt(x1,x2,x3,nrow1,ncol1,ncol2,nr1,
     $     nr2,nr3)
* multiply two matrices of double precision numbers, x3=x1*x2
* arguments:
*     x1,x2 (input), first and second matrix
*     x3 (output), product x1*x2
*     nrow1 (input), number of rows in x1, also the number of rows in x3
*     ncol1 (input), number of columns in x1, also the number of
*          rows in x2
*     ncol2 (input), number of columns in x2, also the number of
*          columns in x3
*     nr1 (input), number of rows in the physical array x1
*     nr2 (input), number of rows in the physical array x2
*     nr3 (input), number of rows in the physical array x3
      implicit none
      integer i,j,k,nrow1,ncol1,ncol2,nr1,nr2,nr3
      double precision x1(nr1,ncol1),x2(nr2,ncol2),x3(nr3,ncol2),
     $     x(nrow1,ncol2)
      x=0
      do i=1,ncol2
      do j=1,nrow1
      do k=1,ncol1
        x(j,i)=x(j,i)+x1(j,k)*x2(k,i)
      enddo
      enddo
      enddo
      x3(1:nrow1,1:ncol2)=x(1:nrow1,1:ncol2)
      end
*****************************************************************************
      subroutine pir_ssg_data_zmatmlt(x1,x2,x3,nrow1,ncol1,ncol2,nr1,
     $     nr2,nr3)
* multiply two matrices of complex*16 numbers, x3=x1*x2
* arguments:
*     x1,x2 (input), first and second matrix
*     x3 (output), product x1*x2
*     nrow1 (input), number of rows in x1, also the number of rows in x3
*     ncol1 (input), number of columns in x1, also the number of
*          rows in x2
*     ncol2 (input), number of columns in x2, also the number of
*          columns in x3
*     nr1 (input), number of rows in the physical array x1
*     nr2 (input), number of rows in the physical array x2
*     nr3 (input), number of rows in the physical array x3
      implicit none
      integer i,j,k,nrow1,ncol1,ncol2,nr1,nr2,nr3
      complex*16 x1(nr1,ncol1),x2(nr2,ncol2),x3(nr3,ncol2),
     $     x(nrow1,ncol2)
      x=0
      do i=1,ncol2
      do j=1,nrow1
      do k=1,ncol1
        x(j,i)=x(j,i)+x1(j,k)*x2(k,i)
      enddo
      enddo
      enddo
      x3(1:nrow1,1:ncol2)=x(1:nrow1,1:ncol2)
      end
*****************************************************************************
      subroutine pir_ssg_data_factor(n,numbers)
* remove the greatest common factor contained in n integers in numbers
      implicit none
      integer n,numbers(n),min,i,j,factor
      factor=1
* find a nonzero integer
      do i=1,n
        if(numbers(i).ne.0)goto 2
      enddo
* all zeros
      return
* find the minimum absolute nonzero value among the integers
2     min=iabs(numbers(i))
      do i=2,n
        if(numbers(i).ne.0.and.iabs(numbers(i)).lt.min)
     +      min=iabs(numbers(i))
      enddo
* try each number from the minimum on down to 2
      do i=min,2,-1
* is i a common factor?
        do j=1,n
          if(mod(numbers(j),i).ne.0)goto 1
        enddo
* yes, divide it out
        do j=1,n
          numbers(j)=numbers(j)/i
        enddo
* save it too
        factor=factor*i
* done
        return
* try next number
1       continue
      enddo
      end





