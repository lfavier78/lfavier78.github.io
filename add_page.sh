#!/bin/bash

for file in *.html
do
 echo $file
 cp $file save/$file
 cat $file | sed 's/<td align=\"center\"><a href=\"files.html\">Files<\/a><\/td>/<td align=\"center\"><a href=\"files.html\">Files<\/a><\/td>\n<tr>\n<td align=\"center\"><a href=\"reading.html\">Reading Seminar<\/a><\/td>\n<\/tr>/g' > tmp.txt
 mv tmp.txt $file
done

#cat cv.html | sed 's/<td align=\"center\"><a href=\"research.html\">Research<\/a><\/td>/<td align=\"center\"><a href=\"research.html\">Research<\/a><\/td>\n<tr>\n<td align=\"center\"><a href=\"projects.html\">Projects<\/a><\/td>\n<\/tr>/g' > tmp.txt
