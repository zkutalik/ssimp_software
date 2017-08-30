
for file in *.pdf; do \
echo $file;\
convert -trim -density 900 -quality 100 $file `echo $file|cut -f1 -d'.'`.jpg;\
done

