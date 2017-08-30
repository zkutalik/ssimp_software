
for file in *.pdf; do \
echo $file;\
convert -trim -density 900 -quality 100 $file `echo $file|cut -f1 -d'.'`.jpg;\
convert `echo $file|cut -f1 -d'.'`.jpg -resize 30% `echo $file|cut -f1 -d'.'`.jpg;\
done
