for p in 1 2 4; do
    for t in 1 2 4 8; do
        mpisubmit.pl -p $p \
            -W 00:10 \
            --stdout /dev/null \
            --stderr /dev/null \
            -t $t \
            ./main_mpi 128 1
        
        mpisubmit.pl -p $p \
            -W 00:10 \
            --stdout /dev/null \
            --stderr /dev/null \
            -t $t \
            ./main_mpi 128
        
        mpisubmit.pl -p $p \
            -W 00:10 \
            --stdout /dev/null \
            --stderr /dev/null \
            -t $t \
            ./main_mpi 256 1
        
        mpisubmit.pl -p $p \
            -W 00:10 \
            --stdout /dev/null \
            --stderr /dev/null \
            -t $t \
            ./main_mpi 256
    done
done