for i in *.xvg;do
    mv "$i" "$(echo $i | tr a-z A-Z)";
done
for i in *.XVG;do
    mv "$i" "${i//TABLE/table}"
done
for i in *.XVG;do    
    mv "$i" "${i//XVG/xvg}"
done
