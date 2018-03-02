#!/bin/sh

for file in $(find ./ -iname operms.2D.*.rds); do 
	out=$(du -sh $file); 
	size=$(echo $out | cut -f1 -d " "); 
	if [ "$size" == "0" ]; then 
		echo $out; 
	fi;  
done
