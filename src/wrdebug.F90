#include "redef.h"
Module wrdebug
  Implicit None
  Integer, Parameter::    tfile = 25

  Interface diag_array
     Module Procedure diag_array_1d, diag_array_2d
  End Interface diag_array

  Interface diag_nob
     Module Procedure&
          diag_2d_nb_real, diag_2d_nb_complex,&
          diag_3d_nb_real, diag_3d_nb_complex,&
          diag_6d_nb_real, diag_6d_nb_complex
  End Interface diag_nob

  Interface diag_onlyb
     Module Procedure diag_3d_onlyb
  End Interface diag_onlyb

  Interface diag_onlyb
  End Interface diag_onlyb

  !---------------- Hilfsroutinen zum debuggen --------------------
  Contains
  Function getoffx(xsize)
    Use par_mod
    Implicit None

    Integer:: getoffx, xsize
    If (xsize == li0) Then
       getoffx = nib
    Else If (xsize == li0) Then
       getoffx = 0
    Else
       Stop "bad x size"
    Endif
  End Function getoffx

  Function getoffz(zsize)
    Use par_mod
    Implicit None

    Integer:: getoffz, zsize
    If (zsize == lz0) Then
       getoffz = nzb
    Else If (zsize == lk0) Then
       getoffz = 0
    Else
       Stop "bad z size"
    Endif
  End Function getoffz

  Function getoffv(vsize)
    Use par_mod
    Implicit None

    Integer:: getoffv, vsize
    If (vsize == lv0) Then
       getoffv = nvb
    Else If (vsize == ll0) Then
       getoffv = 0
    Else
       Stop "bad v size"
    Endif
  End Function getoffv

  Function getoffw(wsize)
    Use par_mod
    Implicit None

    Integer:: getoffw, wsize
    If (wsize == lw0) Then
       getoffw = nwb
    Else If (wsize == lm0) Then
       getoffw = 0
    Else
       Stop "bad w size"
    Endif
  End Function getoffw

  Subroutine killgene
    Integer, Save:: ncalls=0
    Real:: rnull = 0., var
    ncalls = ncalls + 1
    if (ncalls < 9) return
    var = 1./rnull
  End Subroutine killgene

  Subroutine diag_array_1d(name, ph)
    Use par_mod
    Implicit None

    Character(*)::  name
    Complex::       ph(:)
    Character*50::  filename
    Integer:: ii

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do ii=1, Size(ph)
       Write (tfile, "(i2,2e30.12)") ii, ph(ii)
    Enddo

    Close(tfile)
  End Subroutine diag_array_1d

  Subroutine diag_array_2d(name, ph)
    Use par_mod
    Implicit None

    Character(*)::  name
    Complex::       ph(:,:)
    Character*50::  filename
    Integer:: ii, jj

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do ii=1, Size(ph, 1)
       Write (tfile, "(i4,(2e30.12))") ii, ph(ii, :)
    Enddo

    Close(tfile)
  End Subroutine diag_array_2d

  Subroutine diag_2d_nb_real(name, ph)
    Use par_mod
    Implicit None

    Character(*)::  name
    Real::          ph(0:, 0:)
    Character*50::  filename
    Integer::       i,j, offx

    offx = getoffx(Size(ph,1))
    If (Size(ph,2) /= lj0) Stop "diag_xy y"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do j = lj1, lj2
          Write (tfile, "(2i2,e30.12)")&
               i, j, ph(offx+(i-li1), j-lj1)
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_2d_nb_real

  Subroutine diag_2d_nb_complex(name, ph)
    Use par_mod
    Implicit None

    Character(*)::  name
    Complex::       ph(0:, 0:)
    Character*50::  filename
    Integer::       i,j, offx

    offx = getoffx(Size(ph,1))
    If (Size(ph,2) /= lj0) Stop "diag_nob y"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do j = lj1, lj2
          Write (tfile, "(2i2,2e30.12)")&
               i, j, ph(offx+(i-li1), j-lj1)
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_2d_nb_complex

  Subroutine diag_vw(name, ph)
    Use par_mod
    Implicit None

    Character(*)::  name
    Complex::       ph(ll1:, 0:)
    Character*50::  filename
    Integer::   l

    If (Size(ph,1) /= ll0) Stop "diag_vw v"
    If (Size(ph,2) /= lm0) Stop "diag_vw w"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do l = ll1, ll2
       Write (tfile, "(i2,2e30.12)") l, Sum(ph(l,:))
    Enddo

    Close(tfile)
  End Subroutine diag_vw

  Subroutine diag_densmu(name, eflds)
    Use par_mod
    Implicit None

    Character(*)::  name
    Complex::       eflds(li1:, lj1:, lk1:, lm1:, ln1:)
    Character*50::  filename
    Integer::   i, j, k

    If (Size(eflds,1) /= li0) Stop "diag_densmu x"
    If (Size(eflds,2) /= lj0) Stop "diag_densmu y"
    If (Size(eflds,3) /= lk0) Stop "diag_densmu z"
    If (Size(eflds,4) /= lm0) Stop "diag_densmu m"
    If (Size(eflds,5) /= 1) Stop "diag_densmu spec"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do j = lj1, lj2
          Write (tfile, "(2i3,2e30.12)") i, j, Sum(eflds(i, j, :, :, ln1))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_densmu

  Subroutine diag_3d_onlyb(name, phii)
    Use par_mod
    Implicit None

    Character(*):: name
    Complex::      phii(li1:, lj1:, lbz:)
    Character*50:: filename
    Integer::      i, j, k

    If (Size(phii,1) /= li0) Stop "diag_3d x"
    If (Size(phii,2) /= lj0) Stop "diag_3d y"
    If (Size(phii,3) /= lz0) Stop "diag_3d z"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do j = lj1, lj2
          Write (tfile, "(i3,i3,2e30.12)") i, j, Sum(phii(i, j, lk1:lk2))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_3d_onlyb

  Subroutine diag_3d_nb_real(name, phii)
    Use par_mod
    Implicit None

    Character(*):: name
    Real::         phii(0:, 0:, 0:)
    Character*50:: filename
    Integer::      i, j, k, offx, offz

    offx = getoffx(Size(phii,1))
    If (Size(phii,2) /= lj0) Stop "diag_3d_nb y"
    offz = getoffz(Size(phii,3))

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do k = lj1, lj2
          Write (tfile, "(i3,i3,e30.12)") i, j,&
               Sum(phii(offx+(i-li1), j-lj1, offz:offz+lk0-1))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_3d_nb_real

  Subroutine diag_3d_nb_complex(name, phii)
    Use par_mod
    Implicit None

    Character(*):: name
    Complex::      phii(0:, 0:, 0:)
    Character*50:: filename
    Integer::      i, j, k, offx, offz

    offx = getoffx(Size(phii,1))
    If (Size(phii,2) /= lj0) Stop "diag_3d_nb y"
    offz = getoffz(Size(phii,3))

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do j = lj1, lj2
          Write (tfile, "(i3,i3,2e30.12)") i, j,&
               Sum(phii(offx+(i-li1), j-lj1, offz:offz+lk0-1))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_3d_nb_complex

  Subroutine diag_6d_nb_real(name, ff_)
    Use par_mod
    Implicit None

    Character(*):: name
    Real::         ff_(0:,0:,0:,0:,0:,0:)
    Character*50:: filename
    Integer::      i, j, offx, offz, offv, offw

    offx = getoffx(Size(ff_,1))
    If (Size(ff_,2) /= lj0) Stop "diag_6d_nb_real y"
    offz = getoffz(Size(ff_,3))
    offv = getoffz(Size(ff_,4))
    offw = getoffw(Size(ff_,5))
    If (Size(ff_,6) /= 1) Stop "diag_6d_nb spec"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do j = lj1, lj2
          Write (tfile, "(2i3,2e30.12)") i, j,&
               Sum(ff_(offx+(i-li1), j-lj1, offz:offz+lk0-1,&
               offv:offv+ll0-1, offw:offw+lm0-1, 0))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_6d_nb_real

  Subroutine diag_6d_nb_complex(name, ff_)
    Use par_mod
    Implicit None

    Character(*):: name
    Complex::      ff_(0:,0:,0:,0:,0:,0:)
    Character*50:: filename
    Integer::      i, j, k, l, offx, offz, offv, offw

    offx = getoffx(Size(ff_,1))
    If (Size(ff_,2) /= lj0) Stop "diag_6d_nb_complex y"
    offz = getoffz(Size(ff_,3))
    offv = getoffv(Size(ff_,4))
    offw = getoffw(Size(ff_,5))
    If (Size(ff_,6) /= 1) Stop "diag_6d_nb_complex spec"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = li1, li2
       Do j = lj1, lj2
          Write (tfile, "(2i3,2e30.12)") i, j,&
               Sum(ff_(offx+(i-li1), j-lj1, offz:offz+lk0-1,&
               offv:offv+ll0-1, offw:offw+lm0-1, 0))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_6d_nb_complex

  Subroutine diag_fac(name, fac)
    Use par_mod
    Implicit None

    Character(*)::  name
    Real::      fac(:,:,:)
    Character*50::  filename
    Integer::   j, k

    If (Size(fac,1) /= li0) Stop "diag_fac: 1"
    If (Size(fac,2) /= lj0) Stop "diag_fac: 2"
    If (Size(fac,3) /= lk0) Stop "diag_fac: 3"

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do j = 1, lj0
       Do k = 1, lk0
          Write (tfile, "(i3,i3,2e30.12)") j + lj0*my_pey, k,&
               Sum(fac(:, j, k))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_fac


  Subroutine diag_klm(name, klm)
    Use par_mod
    Implicit None

    Character(*)::  name
    Real::      klm(lk1:lk2,ll1:ll2,lm1:lm2)
    Character*50::  filename
    Integer::   l, m

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do l = ll1, ll2
       Do m=lm1,lm2
          Write (tfile, *) l, m, Sum(klm(:,l,m))
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_klm

  Subroutine diag_real(name, vex_re)
    Use par_mod
    Implicit None

    Character(*):: name
    Complex::      vex_re(0:, 0:)
    Character*50:: filename
    Integer::      i, j

    Write(filename, "(a,i1,a)") name, mype, ".dat"
    Open(tfile, file=filename)

    Do i = 0, Ubound(vex_re,2)
       Do j=0, Ubound(vex_re,1)
          Write (tfile, "(2i3,2e30.12)") i, j, vex_re(j,i)
       Enddo
    Enddo

    Close(tfile)
  End Subroutine diag_real

End Module wrdebug
