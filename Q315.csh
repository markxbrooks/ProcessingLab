#!/bin/csh -f
cd mosflm
set count = 0
foreach I ( `grep -l '      PARAMETER (IXWDTH=8192)' *` ) 
if (-e $I) then
@ count = $count + 1
echo First pass for file $I
sed s/"      PARAMETER (IXWDTH=8192)"/"c     PARAMETER (IXWDTH=8192)"/ $I > tmp1
sed s/"c     PARAMETER (IXWDTH=12288)"/"      PARAMETER (IXWDTH=12288)"/ tmp1 > tmp2

echo Second pass
sed s/"      PARAMETER (IYLENGTH=4096)"/"c     PARAMETER (IYLENGTH=4096)"/ tmp2 > tmp1
sed s/"c     PARAMETER (IYLENGTH=6144)"/"       PARAMETER (IYLENGTH=6144)"/ tmp1 > $I
echo "******Finished*******"
endif
end
if ( $count == 0 ) then
	echo no files were changed\; mosflm can already deal with unbinned Q315 images!
else
	echo $count files were changed\; mosflm can now deal with unbinned Q315 images!
endif
/bin/rm -f tmp1
/bin/rm -f tmp2
cd ..
