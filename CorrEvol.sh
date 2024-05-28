dirstart='/media/lorenzo/Dati/Desktop/Linux-File/Open-Quantum-System/Redfield/RedfieldLinear/GeneralDynamic/'
dirnow=`pwd`
cp "${dirstart}CorrelationFunctionGeneral.e" "${dirstart}CORR.sh" $dirnow

if test -z "$1"
then
    ls *.inp
    echo 'Select the input file'
    read inp
else
    inp="$1"
fi

if test -z "$2"
then
    ls *.inp
    echo 'Select the input file for time evolution'
    read timeinp
else
    timeinp="$2"
fi

if test -z "$3"
then
  bash CORR.sh $inp $timeinp
else
    dir="$3"
    echo "Dir found! $dir"
    mkdir $dir
    cp CorrelationFunctionGeneral.e CORR.sh $dir
    listfile=`cat $inp`
    echo $listfile
    for i in $listfile
    do
      echo $i
      if test -z "$i"
      then
        echo "File $i does not exist!"
      else
        cp $i $dir
      fi
    done
    cp $inp $timeinp $dir
    cd $dir
    bash CORR.sh $inp $timeinp
    cd ..
fi

