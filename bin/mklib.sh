code=$1

rm -f lib${code}.so  *.o *Dict.* 

cflag="$(root-config --auxcflags | tr -s ' ' '\n' | grep c++1)  -I$(root-config --incdir) -I$(pwd)"

echo $cflag

echo 
comm="g++ $cflag -fPIC -c ${code}.cxx -o ${code}.o"
echo $comm
eval $comm || exit 1

rootver=$(root-config --version | awk -F\. '{print $1}')
if [ ${rootver} == 5 ]
then
    rcopt="-c"
elif [ ${rootver} == 6 ]
then 
    rcopt=""
else
    echo unknown root version
    exit 1
fi

comm="rootcint -f ${code}Dict.cxx ${rcopt} $cflag  ${code}.h ${code}LinkDef.h"
echo $comm
eval $comm  || exit 1

comm="g++ $cflag  -fPIC -c ${code}Dict.cxx -o ${code}Dict.o"
echo $comm
eval $comm || exit 1


comm="g++ -shared -O3 -Wall -Werror $cflag $(root-config --libs) ${code}.o ${code}Dict.o -o lib${code}.so"
echo $comm
eval $comm
    
rm -f  ${code}.o ${code}Dict.h ${code}Dict.cxx ${code}Dict.o
#${code}Dict_rdict.pcm #this is needed!

echo
echo "=========================     ${code} done!     ========================="
echo

#ln -s $(pwd)/lib${code}.so ../lib/
