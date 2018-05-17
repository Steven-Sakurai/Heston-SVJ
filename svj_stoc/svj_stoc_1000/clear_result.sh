i="1"
dirName="session"
while [ $i -lt 16 ]
do
	name=$dirName$i
	rm $name/conv.csv
	rm $name/est.csv
	i=$[$i+1]
done