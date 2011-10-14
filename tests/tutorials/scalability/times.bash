for i in $*
do 
rm tmp tmp2
grep -B33 'category: WFN_OPT' $i | sed  s/'Total CPU time for category: WFN_OPT'//g | awk '{print $2}' > tmp
tr -s '\n' ' ' <tmp> tmp2
echo $i `cat tmp2`
done 

