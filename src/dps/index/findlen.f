      subroutine findlen(val,np,npkmin,nmultmin,idebug,fit,len)
      implicit none
C
C---- Parameters
C
      integer nbins
      parameter (nbins=10)      
      integer nopoints
      parameter (nopoints=2000)
      integer nopeaks
      parameter (nopeaks=2000)
C
C---- scalar arguments
C
      real fit,len
      integer np,npkmin,nmultmin,idebug
C
C---- ARRAY ARGUMENTS
      real val(*)

C     local arrays
      integer ipk(nopeaks),nbin(nbins),irej(nopoints),nrej(nbins)
      real xpk(nopeaks),valpk(nopeaks),rmsbin(nbins),sumbin(nbins), 
     +     sumsqbin(nbins),avgbin(nbins),rmsdbin(nbins)
C     local scalars
      integer i,j,k,nporg,npk,npkmod,nsplit,ipt,itrue,nstart,
     +        npbin,ibin,npstart
      real vali,sum,sumsq,avg,rms,rmsd,pkmax,threshorg,threshpk,scal,
     +     valmax,sumx,lenx
      logical last, found
C
C---- find mean and rms overall and for 10 different bins so that
C     the threshold can vary across the FFT (set lower for higher terms).
C
      npbin = np/(nbins+1)
C
C---- Initialise
C
      len = 0.0
      fit = 0.0
      sum = 0.0
      sumsq = 0.0

      if (np.le.1) return
C
      do 2 i = 1,nbins
        nbin(i) = 0
        nrej(i) = 0
        sumbin(i) = 0.0
        sumsqbin(i) = 0.0
 2    continue
C
C---- Find overall mean,rms, exclusing origin
C     
      do 10 i = 2,np
        vali = val(i)
        sum = sum + vali
        sumsq = sumsq + vali*vali
 10   continue
C
      avg = sum/(np-1)
      rmsd = sqrt(sumsq/(np-1) - avg*avg)
      rms = sqrt(sumsq/(np-1))
C
C---- flag outliers (ie peaks) to get a better estimate of mean and rms.
C
      do 14 i = 1,np
        irej(i) = 0
        if (abs(val(i)-avg).gt.5*rmsd) irej(i) = 1
 14   continue
C
C---- find mean and rms for each bin and overall, rejecting outliers (peaks).
C
      sum = 0.0
      sumsq = 0.0
      do 16 i = 1,np
        ibin = 1 + i/npbin
        if (ibin.gt.nbins) ibin = nbins
        if (irej(i).gt.0) then
          nrej(ibin) = nrej(ibin) + 1
          goto 16
        end if
        vali = val(i)
        sum = sum + vali
        sumsq = sumsq + vali*vali
        nbin(ibin) = nbin(ibin) + 1
        sumbin(ibin) = sumbin(ibin) + vali
        sumsqbin(ibin) = sumsqbin(ibin) + vali*vali
 16   continue
      do 18 i = 1,nbins
        if (nbin(i).eq.0) goto 18
        avgbin(i) = sumbin(i)/nbin(i)
        rmsdbin(i) = sqrt(sumsqbin(i)/nbin(i) -
     +                       avgbin(i)*avgbin(i))
 18   continue
C
C---- Set threshold for finding extent of origin peak and find
C     corresponding pixel number
C
      threshorg = avgbin(1) + 1.5*rmsdbin(1)
      do 20 I = 1,np
        if (val(i).le.threshorg) then
          nporg = i - 1
          goto 22
        end if
 20   continue
C
C---- Has not found extent of origin peak, return (len=fit=0)
C
      return
C
C---- Find maximum peak beyond this pixel
C
 22   pkmax = 0.0
      do 24 i = nporg+1,np
        if (val(i).gt.pkmax) then
          pkmax = val(i)
        end if
 24   continue
C
      if (idebug.gt.0) write(6,110)avg,rmsd,rms,nporg+1,pkmax,val(1),
     + (avgbin(i),rmsdbin(i),nrej(i),i=1,nbins)
 110  FORMAT(/,1X,'  Mean',F10.2,'  Rmsd',F8.2,
     +       '  Rms',F8.2,/,1X,'Origin peak boundary at pixel',I4,
     +       1X,'Max non-origin peak',F6.0,
     +       '  origin',F6.0,/,1X,'Average, rmsd and nrej in bins',/,
     +       (1x,5(F6.1,F6.1,I5)))
C
C---- Set threshold for finding peaks. Several (close to origin) peaks should
C     be very similar in height to the biggest peak (which may or may not be 
C     the first peak.
C     Need to modify threshold as a function of pixel, as rms gets lower
C     at high pixel values. Could use a variable rms to do this.
C
      scal = 0.75
      i = nporg
      npk = 0
C
 30   threshpk = avg + scal*(pkmax-avg)
      npkmod = 0
      if (idebug.gt.0) write(6,114) scal, i+1
 114  format(//,1X,'Looking for peaks with scal=',F4.2,
     +              ' starting at pixel',I5)
C
C---- Find all peaks above this threshold
C
 32    i = i + 1
       ibin = 1 + i/npbin
       if (ibin.gt.nbins) ibin = nbins
       threshpk = avgbin(ibin) + scal*(pkmax-avg)
       threshpk = max(threshpk, (avgbin(ibin)+2.5*rmsdbin(ibin))) 
       if (val(i).gt.threshpk) then
         npk = npk + 1
C         type *,'found peak npk at pixel',npk,i
 34      j = i + 1
C
C---- Find max value for this peak
C
         if (j.le.np) then
           if (val(j).gt.val(i)) then
             i = i + 1
             goto 34
           else
             ipk(npk) = i
             xpk(npk) = float(i)
             valpk(npk) = val(i)
             goto 32
           end if
         else
           ipk(npk) = i
         end if
       else
         if (i.lt.np) goto 32
       end if

       if (idebug.gt.0) write(6,120) scal, npk,
     +       (ipk(i),xpk(i),valpk(i),i=1,npk)
 120   format(1x,'scal=',F5.3,
     +    '  number of peaks',I5,/,1x,
     +    'values of pixel, x and peak height',
     +    /,(1x,4(i5,f7.1,f6.0,4x)))
C
C---- If FFT is very noisy, it is possible that it will not find even
C     the strongest peak
C
      if (npk.eq.0) goto 48
C
C---- Some peaks can be noisy, ie have more than one maximum. Merge these now
C     using the width of the origin peak as a marker.
C
      last = .false.
      i = 0
 40   i = i + 1
      if (i.lt.npk) then
        if (ipk(i+1)-ipk(i).lt.nporg) then
C
C---- Duplicate found, now find how many peaks are involved
C
          nsplit = 2
          do 42 j = i+1,npk-1
            if (ipk(j+1)-ipk(i).lt.nporg) then
              nsplit = nsplit + 1
            else
              goto 44
            end if
 42       continue
C
C---- Now find average of all points included, weighted by peak value
C
 44       sum = 0.0
          sumx = 0.0
          valmax = 0
          ibin = 1 + ipk(i)/npbin 
          if (ibin.gt.nbins) ibin = nbins
          avg = avgbin(ibin)
          do 46 j = 1,nsplit
            ipt = i + j -1
            valmax = max(valmax,valpk(ipt))
            sum = sum + valpk(ipt) - avg
            sumx = sumx + (valpk(ipt) - avg)*xpk(ipt)
 46       continue
          npkmod = npkmod + 1
          valpk(npkmod) = valmax
          if (sum.gt.0.0) then
            xpk(npkmod) = sumx/sum
          else
            xpk(npkmod) = xpk(i)
          end if
          ipk(npkmod) = ipk(i+nsplit-1)
          write(6,130) nsplit,i,valmax,xpk(npkmod),ipk(npkmod),avg
 130      format(1x,'Merged',I3,' peaks starting at number',I3,'  max'
     +         ' peak height',F6.1,' weighted mean x',F6.1,
     +         ' pixel',I5,' avg',F6.0)
          i = i + nsplit - 1
C
C---- If the very last peak has been included, decrement npkmod or
C     it will be counted twice
          if (i.eq.npk) last = .true.
        else
          npkmod = npkmod + 1
          valpk(npkmod) = valpk(i)
          xpk(npkmod) = xpk(i)
          ipk(npkmod) = ipk(i)
        end if
C
C---- If not yet the last but one peak, got to next peak
C
        if (i.lt.npk-1) goto 40
      end if
C
C---- All peaks tested for multiples. Now Add in the last peak if not
C     already included as a merged peak
C
      if (.not.last) then
        npkmod = npkmod + 1
        valpk(npkmod) = valpk(npk)
        xpk(npkmod) = xpk(npk)
        ipk(npkmod) = ipk(npk)
      end if
      write(6,140) npkmod,
     +       (ipk(i),xpk(i),valpk(i),i=1,npkmod)
 140  format(1x,'After merging multiple peaks',/,1X,
     +    '  number of peaks',I5,/,1x,
     +    'values of pixel,x and peak height',
     +    /,(1x,4(I5,f7.1,f6.0,4x)))
C
C---- If a sufficient number of peaks found, proceed, otherwise lower
C     threshold, and search outwards
C     from the last peak found (plus 2x*origin width).
C
 48   if (npkmod.lt.npkmin) then
        if (scal.gt.0.51) then
          scal = scal - 0.25
        else
          scal = scal - 0.1
        end if
        npk = npkmod
C
C---- If scal is ge 0.25, repeat search from beginning
C     otherwise start at last found point. This avoids
C     picking up noise peaks near origin when scal is very low.
C
        if (scal.gt.0.249) then
          npk = 0
          npkmod = 0
          i = nporg
        else
          if (npkmod.gt.0) then
            i = ipk(npkmod) + 2*nporg
          else
            i = nporg
          end if
        end if
        if (scal.ge.0) goto 30
      end if
C
C---- Sufficient peaks now found. Now have to find the primitive repeat length
C     Starting with the shortest, search list for a peak which has at least 
C     nmultmin multiples
C
      do 50 i = 1,npkmod-nmultmin
        len = xpk(i)
        nstart = i + 1
C
C---- Loop over multiplicities
C
        do 52 j = 1,nmultmin
           lenx = (j+1)*len
           found = .false.
C
C---- Loop over unused peaks looking for this multiplicity
C
           do 54 k = nstart,npkmod
C             type *,'i,j,k,lenx,xpk(k)',i,j,k,lenx,xpk(k),orgwidth
             found = (abs(lenx-xpk(k)).lt.(nporg+j*0.5))
             if (found) then
               write(6,150) i,j+1,k
 150         format(1x,'For peak',I3,' found multiple',I3,
     +              ' at peak',I3)
               nstart = k + 1
               goto 52
             end if
 54        continue
           if (.not.found) goto 50
 52     continue
        if (found) then
          itrue = i
          goto 60
        end if
 50   continue
C
C---- NOT FOUND
C
      if (idebug.gt.0) write(6,152)
 152  format(1x,'** FAILURE **, no suitable peaks found')
      len = 0.0
      return
 60   len = xpk(itrue)
      if (idebug.gt.1) write(6,154) itrue,nmultmin,len
 154  format(//,1X,'** SUCCESS **Peak',I3,' has the required',I3,
     +          ' multiples, len=',F6.1)
      return
      end
