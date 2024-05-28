if test -z "$1"
then
    ls *.inp
    echo "Select the INPUT"
    read INPUT
else
    INPUT=$1
fi

if test -z "$2"
then
    ls *.inp
    echo "Select the INPUT"
    read INPUT2
else
    INPUT2=$2
fi
OUTT=${INPUT%'.inp'}.out

./CorrelationFunctionGeneral.e <<EOF > $OUTT &
$INPUT
$INPUT2
EOF
sleep 1


   
