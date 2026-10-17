#!/bin/csh -f
cd mosflm
set count = 0
foreach I ( `grep -l 'c     PARAMETER (IXWDTH=8192)' *` ) 
if (-e $I) then
@ count = $count + 1
echo First pass for file $I
sed s/"c     PARAMETER (IXWDTH=8192)"/"      PARAMETER (IXWDTH=8192)"/ $I > tmp1
sed s/"      PARAMETER (IXWDTH=12288)"/"c     PARAMETER (IXWDTH=12288)"/ tmp1 > tmp2

echo Second pass
sed s/"c     PARAMETER (IYLENGTH=4096)"/"      PARAMETER (IYLENGTH=4096)"/ tmp2 > tmp1
sed s/"      PARAMETER (IYLENGTH=6144)"/"c      PARAMETER (IYLENGTH=6144)"/ tmp1 > $I
echo "******Finished*******"
endif
end
if ( $count == 0 ) then
	echo no files were changed
else
	echo $count files were changed
endif

echo mosflm cannot now deal with unbinned Q315 images!
/bin/rm -f tmp1
/bin/rm -f tmp2
cd ..
