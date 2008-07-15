*=======================================================================
*     Copyright (C) 1999,2001,2002
*     Associated Universities, Inc. Washington DC, USA.
*
*     This library is free software; you can redistribute it and/or
*     modify it under the terms of the GNU Library General Public
*     License as published by the Free Software Foundation; either
*     version 2 of the License, or (at your option) any later version.
*
*     This library is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Library General Public License for more details.
*
*     You should have received a copy of the GNU Library General Public
*     License along with this library; if not, write to the Free
*     Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
*     MA 02139, USA.
*
*     Correspondence concerning AIPS++ should be addressed as follows:
*            Internet email: aips2-request@nrao.edu.
*            Postal address: AIPS++ Project Office
*                            National Radio Astronomy Observatory
*                            520 Edgemont Road
*                            Charlottesville, VA 22903-2475 USA
*
*     $Id: fpbfridft.f,v 1.3 2005/11/10 18:12:00 cvsmgr Exp $
*-----------------------------------------------------------------------
C
C Compute the convolution function as IdealEij*fij, where fij has the
C first order approximation of the pointing errors.
C

      double precision function scaling(sigma,support,du,dv)
      implicit none
      double precision sigma,du,dv
      double precision sum,u,v
      integer i,j,support
      double precision pi
      data pi/3.14159265358979323846/

      sum = 0.0
      do i=-support,support
         u = i*du
         do j=-support,support
            v = j*dv
            sum = sum+exp(-(u**2+v**2)*pi**2/sigma)
         enddo
      enddo
      scaling = sum
      return
      end

      complex  function  Eij(griduvw,sigma,area,loc,
     $     raoff1, decoff1, raoff2, decoff2, 
     $     IdealEij, convsize, wconvsize,dovp)
      implicit none
      integer convsize,wconvsize,loc(3),dovp
      double precision raoff1, raoff2, decoff1, decoff2,
     $     scale(3),lambda,sigma,area
      double precision griduvw(3)
      complex IdealEij(convsize/2-1, convsize/2-1, wconvsize)
      complex Kij
      double precision Rij,Eij0, u,v
      double precision pi
      complex iota
      data pi/3.14159265358979323846/
      data iota/(0,1)/

      if (dovp .eq. 1) then
         u = griduvw(1)
         v = griduvw(2)

         Rij=((raoff1-raoff2)**2 + (decoff1-decoff2)**2)*sigma/2
         Kij=pi*iota*(u*(raoff1+raoff2) + v*(decoff1+decoff2))
         Eij0 = (u**2 + v**2)*pi**2/sigma

C         Eij = IdealEij(loc(1),loc(2),1)*(1.0-Rij)*(1.0-Kij)
         Eij = exp(-Eij0-Rij+Kij)/area
      else
         Eij = exp(-Eij0)/area
      endif

c      write(*,*) Eij0,Rij,Kij,exp(-Eij0-Rij+Kij)/area

      return 
      end
C
C Grid a number of visibility records
C
      subroutine pbgridft (uvw, dphase, values, nvispol, nvischan,
     $     gvalues,gnvispol,gnvischan,
     $     dopsf, flag, rflag, weight, nrow, rownum,
     $     scale, offset, grid, nx, ny, npol, nchan, freq, c,
     $     support, convsize, sampling, wconvsize, convfunc, 
     $     chanmap, polmap, sumwt,
     $     ant1, ant2, nant, scanno, sigma,raoff, decoff)

      implicit none
      integer nx, ny, npol, nchan, nvispol, nvischan, nrow
      integer gnvispol, gnvischan
      complex values(nvispol, nvischan, nrow)
      complex gvalues(gnvispol, gnvischan,nrow)
      complex grid(nx, ny, npol, nchan)
      double precision uvw(3, nrow), freq(nvischan), c, scale(3),
     $     offset(3), griduvw(3), lambda

      integer nant, scanno, ant1(nrow), ant2(nrow)
      real raoff(nant), decoff(nant)
      double precision sigma,area

      double precision dphase(nrow), uvdist
      complex phasor
      integer flag(nvispol, nvischan, nrow)
      integer rflag(nrow)
      real weight(nvischan, nrow), phase
      double precision sumwt(npol, nchan)
      integer rownum
      integer support(*), rsupport
      integer chanmap(nchan), polmap(npol)
      integer dopsf

      complex nvalue

      integer convsize, sampling, wconvsize
      complex convfunc(convsize/2-1, convsize/2-1, wconvsize)
      complex cwt,dcwt1,dcwt2

      real norm
      real wt

      logical pbogridft
      complex Eij
      external cppEij
      double precision scaling

      real pos(3)
      integer loc(3), off(3), iloc(3)
      integer rbeg, rend
      integer ix, iy, ipol, ichan
      integer apol, achan, irow
      double precision pi
      data pi/3.14159265358979323846/
      integer mysupport
      double precision mysigma, ra1,ra2,dec1,dec2

      irow=rownum

      if(irow.ge.0) then
         rbeg=irow+1
         rend=irow+1
      else 
         rbeg=1
         rend=nrow
      end if

      mysupport=10
      area = scaling(sigma,mysupport,1/scale(1),1/scale(2))
      do irow=rbeg, rend
         if(rflag(irow).eq.0) then 
            do ichan=1, nvischan
               achan=chanmap(ichan)+1

               lambda = c/freq(ichan)

               if((achan.ge.1).and.(achan.le.nchan).and.
     $              (weight(ichan,irow).gt.0.0)) then
c$$$                  call pbswproj(uvw(1,irow), dphase(irow), freq(ichan), c, 
c$$$     $                 scale, offset, sampling, pos, loc, off, phasor)
                  call pbsgridft(uvw(1,irow), dphase(irow), freq(ichan), 
     $                 c, scale, offset, sampling, pos, loc, off, 
     $                 phasor)
                  iloc(3)=max(1, min(wconvsize, loc(3)))
                  rsupport=support(iloc(3))
                  rsupport=mysupport
C                  if (pbogridft(nx, ny, wconvsize, loc, rsupport)) then
                  if (pbogridft(nx, ny, loc, rsupport)) then
                     do ipol=1, nvispol
                        apol=polmap(ipol)+1
                        if((flag(ipol,ichan,irow).ne.1).and.
     $                       (apol.ge.1).and.(apol.le.npol)) then
C     If we are making a PSF then we don't want to phase
C     rotate but we do want to reproject uvw
                           if(dopsf.eq.1) then
                              nvalue=cmplx(weight(ichan,irow))
                           else
                              nvalue=weight(ichan,irow)*
     $                             (values(ipol,ichan,irow)*phasor)
                           end if
C     norm will be the value we would get for the peak
C     at the phase center. We will want to normalize 
C     the final image by this term.
                           norm=0.0
                           do iy=-rsupport,rsupport
                              iloc(2)=1+abs(iy*sampling+off(2))
                              do ix=-rsupport,rsupport
                                 iloc(1)=1+abs(ix*sampling+off(1))
                                 griduvw(1) = ((loc(1)-offset(1)
     $                                +ix*sampling-2)
     $                                /scale(1))
                                 griduvw(2) = ((loc(2)-offset(2)
     $                                +iy*sampling-2)
     $                                /scale(2))

c$$$                                 cwt = Eij(griduvw,sigma,area,iloc,
c$$$     $                                raoff(ant1(irow)+1), 
c$$$     $                                decoff(ant1(irow)+1),
c$$$     $                                raoff(ant2(irow)+1), 
c$$$     $                                decoff(ant2(irow)+1),
c$$$     $                                convfunc,convsize,wconvsize,0)
                                 ra1 = raoff(ant1(irow)+1)
                                 ra2 = raoff(ant2(irow)+1)
                                 dec1= decoff(ant1(irow)+1)
                                 dec2= decoff(ant2(irow)+1)
                                 call cppEij(griduvw,mysigma,area,iloc,
     $                                ra1,dec1,ra2,dec2,
     $                                convsize,wconvsize,1,
     $                                cwt,dcwt1,dcwt2)
C                                 write(*,*) griduvw(1),griduvw(2),
C     $                                real(cwt),
C     $                                atan2(imag(cwt),real(cwt)),
C     $                                uvw(1,irow)/lambda,
C     $                                uvw(2,irow)/lambda

c$$$                                 if(uvw(3,irow).gt.0.0) then
c$$$                                    cwt=conjg(convfunc(iloc(1),
c$$$     $                                   iloc(2), iloc(3)))
c$$$                                 else
c$$$                                    cwt=convfunc(iloc(1),
c$$$     $                                   iloc(2), iloc(3))
c$$$                                 end if
                                 grid(loc(1)+ix,
     $                                loc(2)+iy,apol,achan)=
     $                                grid(loc(1)+ix,
     $                                loc(2)+iy,apol,achan)+
     $                                nvalue*cwt
                                 norm=norm+real(cwt)
                              end do
C                              write(*,*)
                           end do
C                           stop
                           sumwt(apol,achan)=sumwt(apol,achan)+
     $                          weight(ichan,irow)*norm
                        end if
                     end do
                  else
C     write(*,*) uvw(3,irow), pos(1), pos(2), pos(3),
C     $                 loc(1), loc(2), loc(3)
                  end if
               end if
            end do
         end if
      end do
      return
      end
C     
C Degrid a number of visibility records
C
      subroutine dpbgridft (uvw, dphase, values, nvispol, nvischan,
     $     gazvalues, gelvalues, gnvispol, gnvischan,doconj,
     $     flag, rflag, nrow, rownum, scale, offset, grid, 
     $     nx, ny, npol, nchan, freq, c, support, convsize, sampling, 
     $     wconvsize, convfunc, chanmap, polmap, ant1, ant2, nant, 
     $     scanno, sigma, raoff, decoff,area,dograd)
C ,Eij)

      implicit none
      integer nx, ny, npol, nchan, nvispol, nvischan, nrow,doconj,dograd
      complex values(nvispol, nvischan, nrow)
      integer gnvispol, gnvischan
      complex gazvalues(gnvispol, gnvischan, nrow)
      complex gelvalues(gnvispol, gnvischan, nrow)
      complex grid(nx, ny, npol, nchan)
      double precision uvw(3, nrow), freq(nvischan), c, scale(3),
     $     offset(3), griduvw(3), lambda
      double precision dphase(nrow), uvdist

      integer nant, ant1(nrow), ant2(nrow),scanno
      real raoff(nant), decoff(nant)
      double precision sigma

      complex phasor
      integer flag(nvispol, nvischan, nrow)
      integer rflag(nrow)
      integer rownum
      integer support(*), rsupport
      integer chanmap(*), polmap(*)

      complex nvalue,ngazvalue,ngelvalue

      integer convsize, wconvsize, sampling
      complex convfunc(convsize/2-1, convsize/2-1, wconvsize),
     $     cwt,dcwt1,dcwt2

      real norm, phase

      logical pbogridft
      complex Eij
      external cppEij
      double precision scaling
      real pos(3)
      integer loc(3), off(3), iloc(3)
      integer rbeg, rend
      integer ix, iy, ipol, ichan
      integer apol, achan, irow
      real wt, wtx, wty
      double precision pi,area
      data pi/3.14159265358979323846/
      integer mysupport
      double precision mysigma,ra1,ra2,dec1,dec2
      integer rttt,ttt

      irow=rownum

      if(irow.ge.0) then
         rbeg=irow+1
         rend=irow+1
      else 
         rbeg=1
         rend=nrow
      end if
C
      sampling=1
      mysupport=10
      mysupport=support(1)
      mysigma=sigma*4
c      area=scaling(mysigma,mysupport,1.0/scale(1),1.0/scale(2))

c      write(*,*) area

      do irow=rbeg, rend
         rttt = rflag(irow)
         ttt  = flag(1,1,irow)
         if(rflag(irow).eq.0) then
         do ichan=1, nvischan
            achan=chanmap(ichan)+1
            
            lambda = c/freq(ichan)

            if((achan.ge.1).and.(achan.le.nchan)) then
c$$$               call pbswproj(uvw(1,irow), dphase(irow), freq(ichan), c,
c$$$     $              scale, offset, sampling, pos, loc, off, phasor)
               call pbsgridft(uvw(1,irow), dphase(irow), freq(ichan), c,
     $              scale, offset, sampling, pos, loc, off, phasor)
               iloc(3)=max(1, min(wconvsize, loc(3)))
               rsupport=support(iloc(3))
               rsupport=mysupport
C               if (owproj(nx, ny, wconvsize, loc, rsupport)) then
               if (pbogridft(nx, ny, loc, rsupport)) then

                  do ipol=1, nvispol
                     apol=polmap(ipol)+1
                     if((flag(ipol,ichan,irow).ne.1).and.
     $                    (apol.ge.1).and.(apol.le.npol)) then
                        
                        nvalue=0.0
                        ngazvalue=0.0
                        ngelvalue=0.0
                        do iy=-rsupport,rsupport
                           iloc(2)=1+abs(iy*sampling+off(2))
                           do ix=-rsupport,rsupport
                              iloc(1)=1+abs(ix*sampling+off(1))
                              griduvw(1) = (loc(1)-offset(1)+ix-1)
     $                             /scale(1)-uvw(1,irow)/lambda
                              griduvw(2) = (loc(2)-offset(2)+iy-1)
     $                             /scale(2)-uvw(2,irow)/lambda
                              ra1 = raoff(ant1(irow)+1)
                              ra2 = raoff(ant2(irow)+1)
                              dec1= decoff(ant1(irow)+1)
                              dec2= decoff(ant2(irow)+1)
                              call cppEij(griduvw,mysigma,area,iloc,
     $                             ra1,dec1,ra2,dec2,
     $                             convsize,wconvsize,dograd,
     $                             cwt,dcwt1,dcwt2)
c                              write(*,*) cwt,ant1(irow)+1,ant2(irow)+1
C                              nvalue=nvalue+conjg(cwt)*
                              nvalue=nvalue+(cwt)*
     $                             grid(loc(1)+ix,loc(2)+iy,apol,achan)
                              if ((doconj .eq. 1) .and. 
     $                             (dograd .eq. 1)) then
                                 dcwt1 = conjg(dcwt1)
                                 dcwt2 = conjg(dcwt2)
                              endif
                              if (dograd .eq. 1) then
                                 ngazvalue=ngazvalue+(dcwt1)*
     $                                (grid(loc(1)+ix,loc(2)+iy,
     $                                apol,achan))
                                 ngelvalue=ngelvalue+(dcwt2)*
     $                                (grid(loc(1)+ix,loc(2)+iy,
     $                                apol,achan))
                              endif
                           end do
c                           if (irow .eq. 20) write(*,*)
                        end do
c                        if (irow .eq. 20) write(*,*)
                        values(ipol,ichan,irow)=nvalue*conjg(phasor)
                        if (dograd .eq. 1) then
                           gazvalues(ipol,ichan,irow)=ngazvalue*
     $                          conjg(phasor)
                           gelvalues(ipol,ichan,irow)=ngelvalue*
     $                          conjg(phasor)
                        endif
                     end if
                  end do
c                  if (irow .eq.20) stop
               end if
            end if
         end do
         end if
      end do
      return
      end
C
C Calculate gridded coordinates and the phasor needed for
C phase rotation. 
C
      subroutine pbswproj (uvw, dphase, freq, c, scale, offset, 
     $     sampling, pos, loc, off, phasor)
      implicit none
      integer loc(3), off(3), sampling
      double precision uvw(3), freq, c, scale(3), offset(3)
      real pos(3)
      double precision dphase, phase
      complex phasor
      integer idim
      double precision pi
      data pi/3.14159265358979323846/

C      pos(3)=(scale(3)*uvw(3)*freq/c)*(scale(3)*uvw(3)*freq/c)
C     $     +offset(3)+1.0
C      pos(3)=(scale(3)*uvw(3)*freq/c)+offset(3)+1.0
      pos(3)=sqrt(abs(scale(3)*uvw(3)*freq/c))+offset(3)+1.0
      loc(3)=nint(pos(3))
      off(3)=0

      do idim=1,2
         pos(idim)=scale(idim)*uvw(idim)*freq/c+
     $        (offset(idim)+1.0)
         loc(idim)=nint(pos(idim))
         off(idim)=nint((loc(idim)-pos(idim))*sampling)
      end do

      phase=-2.0D0*pi*dphase*freq/c
      phasor=cmplx(cos(phase), sin(phase))
      return 
      end

c$$$      logical function owproj (nx, ny, nw, loc, support)
c$$$      implicit none
c$$$      integer nx, ny, nw, loc(3), support
c$$$      owproj=(support.gt.0).and.
c$$$     $     (loc(1)-support.ge.1).and.(loc(1)+support.le.nx).and.
c$$$     $     (loc(2)-support.ge.1).and.(loc(2)+support.le.ny).and.
c$$$     $     (loc(3).ge.1).and.(loc(3).le.nw)
c$$$      return
c$$$      end
C
C Calculate gridded coordinates and the phasor needed for
C phase rotation.
C
      subroutine pbsgridft (uvw, dphase, freq, c, scale, offset,
     $     sampling, pos, loc, off, phasor)
      implicit none
      integer sampling
      integer loc(2), off(2)
      double precision uvw(3), freq, c, scale(2), offset(2)
      real pos(2)
      double precision dphase, phase
      complex phasor
      integer idim
      double precision pi
      data pi/3.14159265358979323846/

      do idim=1,2
         pos(idim)=scale(idim)*uvw(idim)*freq/c+(offset(idim)+1.0)
         loc(idim)=nint(pos(idim))
         off(idim)=nint((loc(idim)-pos(idim))*sampling)
      end do

      phase=-2.0D0*pi*dphase*freq/c
      phasor=cmplx(cos(phase), sin(phase))
      return 
      end
C
C Is this on the grid?
C
      logical function pbogridft (nx, ny, loc, support)
      implicit none
      integer nx, ny, loc(2), support
      pbogridft=(loc(1)-support.ge.1).and.(loc(1)+support.le.nx).and.
     $        (loc(2)-support.ge.1).and.(loc(2)+support.le.ny)
      return
      end
