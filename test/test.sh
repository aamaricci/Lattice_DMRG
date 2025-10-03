#Code to run all tests in test/bin/
#N.B. test should end by .x

set -e

if [ ! -d "bin" ]
then
    echo "\e[31m ERROR \e[0m"
    echo " There is no *bin* directory"
    echo " Try  'make all' before testing"
    return 1
fi

CHECK_LIB=$(pkg-config --libs dmrg)
if [ -z ${CHECK_LIB} ]
then
    echo "\e[31m ERROR \e[0m"
    echo " EDIpack not loaded"
    return 1
fi

WITH_MPI=$(pkg-config --variable=mpi dmrg)

cd bin/
HERE=`pwd`
echo $HERE



#Clean the list_dir input list to avoid repetitions.
awk '!seen[$0]++' list_dir > tmp
mv tmp list_dir


while read DIR; do
    if [ -d $DIR ]; then
	if ls $DIR/*.x  1> /dev/null 2>&1; then
	    echo "TESTING $DIR"
	    cd $DIR
	    pwd
	    if [ -f DONE.out ];then
    		echo "Test $DIR has already been passed. Skip"
    	    else
		for exe in *.x
		do
		    echo "Running $exe:"
		    if [ -z ${WITH_MPI} ]
    		    then
    			echo "./$exe "
    			./$exe 
    			echo ""
    			echo ""
    			echo "leaving $DIR..."
    			echo ""
    			echo ""
    			touch DONE.out
    		    else
    			echo "mpiexec -np 2 ./$exe "
    			mpiexec -np 2 ./$exe  < /dev/null
    			echo ""
    			echo ""
    			echo "leaving $DIR..."
    			echo ""
    			echo ""
    			touch DONE.out
    		    fi
		done
	    fi
	    cd $HERE
	fi
    fi
done<list_dir
    

