mess=$1
code=$2

echo message ${mess}
echo code ${code}
echo

if [ $# -eq 2 ]
then
    git add ${code}
    git commit -m "${mess}"
else
    git commit ./ -m "${mess}"
fi

#no more pushing origin
#git push origin local
