echo
echo
echo "              Remember pull first!!!"
echo "              Remember to merge xlu <-> master"
echo "              Remember to push origin master"
echo

echo "===================================== " Branch...
git branch

echo "===================================== " Status...
git status

exit

for bb in master 
# xlu
do
    echo
    echo
    echo "**********************************************"  In ${bb} "**********************************************"
    git checkout ${bb}
    
    echo "===================================== " Branch...
    git branch
    
    echo "===================================== " Status...
    git status
done


