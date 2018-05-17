echo '' > ../log
dirName="session"
i="1"
while [ $i -lt 16 ]
do
    name=$dirName$i
    cd $name
    echo $name
    cat test.Rout >> ../log
    cd ..
    i=$[$i+1]
done