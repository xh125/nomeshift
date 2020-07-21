program imagfreqcorr
    !used to correction the imaginary frequencies by nomemode shift in VASP
    !need file POSCAR, OUTCAR with frequencies calculation
    !POSCAR(new)=POSCAR(old)+scal*natom*Σ[(j=1,njmode)((freq_j/i)/freq_j_max)*Q_j(R)]

	implicit none
        integer,parameter :: maxlen=80
        integer,parameter :: dp=kind(0.0d0)
        integer :: io_error,istatus
        integer:: natomstype,natom,nmode,nimagemode,initimagemode
        integer:: iatomstype,iatom,imode,iimagemode
        character(len=maxlen):: ctemp,POSCARfile
        character(len=maxlen) :: POSCAR_name,OUTCAR_name
        integer :: POSCAR_unit,OUTCAR_unit
        real(kind=dp) :: Latt(3,3)
        character(len=5),allocatable :: Symatom(:)
        integer,allocatable :: numatom(:)
        real(kind=dp),allocatable :: atomposition_0(:,:)
        real(kind=dp),allocatable :: atomposition_cart(:,:)
        real(kind=dp),allocatable :: phononvecter(:,:,:)
        real(kind=dp),allocatable :: phononfreq(:)
        real(kind=dp) :: scal = 0.01

        
    call readOUTCAR()
    
    scal = real(natom)*scal/phononfreq(nmode)
    
    if(initimagemode<=nmode) then
        do imode=initimagemode,nmode
            atomposition_cart=atomposition_cart+scal*phononfreq(imode)*phononvecter(:,:,imode)
        enddo
    endif
    
    call writePOSCAR()
    
    deallocate(atomposition_0,atomposition_cart,phononvecter,phononfreq)
    
    contains
    subroutine readPOSCAR()
        implicit none
        integer :: i,j
        character(len=maxlen) :: commentposcar
        real(kind=dp) :: scaling
        character(len=1) :: coordtype
        POSCAR_name='POSCAR'
        POSCAR_unit= io_file_unit()
        call open_file(POSCAR_name,POSCAR_unit)
        read(POSCAR_unit,*) commentposcar
        read(POSCAR_unit,*) scaling
        do i=1,3
            read(POSCAR_unit,*) (Latt(j,i),j=1,3)
        enddo
        Latt = Latt*scaling
        
        !natomstype=1
        !do while (io_error <= 0)
        !    allocate(Symatom(natomstype))
        !    read(POSCAR_unit,fmt=*,iostat=io_error) (Symatom(i),i=1,natomstype)
        !    natomstype = natomstype+1
        !    deallocate(Symatom)
        !    backspace(POSCAR_unit)
        !enddo
        
        !natomstype = natomstype - 2
        !allocate(Symatom(natomstype),stat=istatus)
        allocate(numatom(natomstype),stat=istatus)
        
        read(POSCAR_unit,fmt=*,iostat=io_error) (Symatom(i),i=1,natomstype)
        read(POSCAR_unit,fmt=*,iostat=io_error) (numatom(i),i=1,natomstype)
        natom=0
        do i=1,natomstype
           natom=natom+numatom(i)
        enddo
        allocate(atomposition_0(3,natom),stat=istatus)
        allocate(atomposition_cart(3,natom),stat=istatus)
        
        read(POSCAR_unit,*) coordtype
        if(coordtype=='S' .or. coordtype=='s') then
            read(POSCAR_unit,*) coordtype
        endif
        
        if ((coordtype=='D') .or. (coordtype=='d')) then
            do i=1,natom
                read(POSCAR_unit,*) atomposition_0(:,i)
                atomposition_cart(:,i)=matmul(atomposition_0(:,i),Latt)
            enddo
        elseif(coordtype=='c' .or. coordtype=='C' .or. coordtype=='K' .or. coordtype=='k') then
            do i=1,natom
                read(POSCAR_unit,*) atomposition_cart(:,i)
            enddo
        endif
        
        call close_file(POSCAR_name,POSCAR_unit)
        
        
    end subroutine
    
    subroutine writePOSCAR()
        implicit none
        character(len=maxlen) :: newposcarname
        integer:: newposcarunit
        integer :: i,j
        
        newposcarname = "POSCAR_shift"
        newposcarunit = io_file_unit()
        call open_file(newposcarname,newposcarunit)
        
        write(newposcarunit,*) "POSCAR shift with image noremel mode"
        write(newposcarunit,"(F12.6)") 1.00
        do i=1,3
            write(newposcarunit,"(3F20.12)") (Latt(j,i),j=1,3)
        enddo
        write(newposcarunit,"(*(A5))") Symatom(:)
        write(newposcarunit,"(*(I5))") numatom(:)
        
        write(newposcarunit,*) "Caracter"
        
        do iatom=1,natom
            write(newposcarunit,"(3F16.9)") (atomposition_cart(j,iatom),j=1,3)
        enddo
        
        call close_file(newposcarname,newposcarunit)
    
    end subroutine
    
    subroutine readOUTCAR()
        implicit none
        logical :: findincar= .false.,findpotcat=.false.,findeigvector = .false.
        character(len=maxlen) :: ctemp1,ctemp2,ctemp3,ctemp4
        integer :: i
        OUTCAR_name = 'OUTCAR'
        OUTCAR_unit = io_file_unit()
        call open_file(OUTCAR_name,OUTCAR_unit)
        do while(findincar /= .true.)
            read(OUTCAR_unit,*) ctemp
            if(trim(adjustl(ctemp))=="INCAR:") findincar=.true.
        enddo
        
        findpotcat = .true.
        natomstype = 0
        do while (findpotcat ==.true.)
            read(OUTCAR_unit,*) ctemp
            if(trim(adjustl(ctemp))=="POTCAR:") then
                natomstype = natomstype + 1
            else
                findpotcat = .false.
            endif
        enddo

        do i=1,natomstype
            backspace(OUTCAR_unit)
        enddo
        backspace(OUTCAR_unit)
        natomstype = natomstype -1    
        allocate(Symatom(natomstype))
        do i=1,natomstype
            read(OUTCAR_unit,*) ctemp1,ctemp2,ctemp3,ctemp4
            Symatom(i)=trim(adjustl(ctemp3))
        enddo
        
        call readPOSCAR()
        
        do while (findeigvector /= .true.)
            read(OUTCAR_unit,*) ctemp
            if(trim(adjustl(ctemp))=="Eigenvectors") then
                findeigvector = .true.
            endif
        enddo
        
        nmode = natom*3
        allocate(phononfreq(nmode))
        allocate(phononvecter(3,natom,nmode))
        read(OUTCAR_unit,"(A80,/)") ctemp

        
        initimagemode =1
        do imode=1,nmode
            read(OUTCAR_unit,"(/,A5,A3,A1,A12)") ctemp1,ctemp2,ctemp3,ctemp4
            read(ctemp4,"(F12.6)") phononfreq(imode)
            
            if(trim(adjustl(ctemp2))=='f') then
                initimagemode = initimagemode +1
            endif
            !if(trim(adjustl(ctemp2))=='f/i' .and. phononfreq(imode)<= imagethro) then
            !    initimagemode = initimagemode+1
            !endif
            
            read(OUTCAR_unit,*) ctemp
            do iatom=1,natom
                read(OUTCAR_unit,*) atomposition_cart(:,iatom),phononvecter(:,iatom,imode)
            enddo

        enddo
        
        call close_file(OUTCAR_name,OUTCAR_unit)
        
        
    end subroutine
  
  function io_file_unit() !得到一个当前未使用的unit，用于打开文件
  !==========================================                                     
  !! Returns an unused unit number
  !! so we can later open a file on that unit.                                       
  !==========================================
  implicit none

    integer :: io_file_unit,unit_index
    logical :: file_open

    unit_index = 9
    file_open  = .true.
    do while ( file_open )
      unit_index = unit_index + 1
      inquire( unit=unit_index, OPENED = file_open ) !用于检查文件状态，参考P.536
    end do
    
    io_file_unit = unit_index

    return
  
  end function io_file_unit
  
  subroutine open_file(file_name,file_unit)
    implicit none
    
    character(len=*),intent(in) :: file_name
    integer,intent(in)          :: file_unit
    integer :: ierr
    character(len=maxlen) :: msg
    open(unit=file_unit, file=file_name,iostat=ierr,iomsg=msg)
    !if(ierr /= 0 ) then
    !  call io_error('Error: Problem opening "'//trim(adjustl(file_name))//' " file')
    !  call io_error(msg)
    !endif
  end subroutine open_file
  
  subroutine close_file(file_name,file_unit)
    implicit none
    
    integer,intent(in)          :: file_unit
    character(len=*),intent(in) :: file_name
    integer :: ierr
    character(len=maxlen) :: msg
    close(file_unit,iostat=ierr,iomsg=msg)
    !if(ierr /= 0 ) then
    ! call io_error('Error: Problem close "'//trim(adjustl(file_name))//' " file')
    !  call io_error(msg)
    !endif   
  end subroutine close_file

end program	
    

