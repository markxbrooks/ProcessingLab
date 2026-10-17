#!/bin/csh -f
cd mosflm
foreach I ( `grep -l 'PARAMETER (BORDER = 3)' *f` )
echo First pass
sed s/"PARAMETER (BORDER = 3)"/"PARAMETER (BORDER = 4)"/ $I > tmp1

echo Second pass
sed s/"PARAMETER (MAX_PIXEL = 600)"/"PARAMETER (MAX_PIXEL = 800)"/ tmp1 > tmp2

echo Third pass
sed s/"$           PAR_TITLE = 1, PAR_FONT = 2, PAR_MENU = 0"/"$           PAR_TITLE = 1, PAR_FONT = 2, PAR_MENU = 0)"/ tmp2 > tmp1

echo Fourth pass
sed s/"PARAMETER (MAX_MEN_ITMS = 19, MAX_MEN_NAME = 19, MEN_FONT = 2,"/"PARAMETER (MAX_MEN_ITMS = 19, MAX_MEN_NAME = 19, MEN_FONT = 4,"/ tmp2 > tmp1

echo Fifth pass
sed s/"$           NOT_TITLE = +1, NOT_FONT = 1, NOT_MENU = 0)"/"$           NOT_TITLE = +1, NOT_FONT = 2, NOT_MENU = 0)"/ tmp1 > tmp2

echo Sixth pass
sed s/"PARAMETER (MAX_PMN_ITMS = 0, MAX_PMN_NAME = 12, PMN_FONT = 2,"/"PARAMETER (MAX_PMN_ITMS = 0, MAX_PMN_NAME = 12, PMN_FONT = 4,"/ tmp2 > tmp1

echo Seventh pass
sed s/"PARAMETER (PPB_FONT = 2, PPB_COLR=3)"/"PARAMETER (PPB_FONT = 4, PPB_COLR=3)"/ tmp1 > tmp2

echo Eighth pass
sed s/"PARAMETER (BUSY_FONT = 2)"/"PARAMETER (BUSY_FONT = 3)"/ tmp2 > $I

echo "******Finished $I *******"
end
echo mosflm display is now big!
/bin/rm -f tmp1
/bin/rm -f tmp2
cd ..
