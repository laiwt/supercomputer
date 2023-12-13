for t in 1 2 4 8 16; do
    bsub -n 1 \
        -q normal \
        -W 00:10 \
        -o /dev/null \
        -e /dev/null \
        -R "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" \
        OMP_NUM_THREADS=$t \
        ./main 128 1

    bsub -n 1 \
        -q normal \
        -W 00:10 \
        -o /dev/null \
        -e /dev/null \
        -R "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" \
        OMP_NUM_THREADS=$t \
        ./main 128
done

for t in 1 4 8 16 32; do
    bsub -n 1 \
        -q normal \
        -W 00:10 \
        -o /dev/null \
        -e /dev/null \
        -R "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" \
        OMP_NUM_THREADS=$t \
        ./main 256 1

    bsub -n 1 \
        -q normal \
        -W 00:10 \
        -o /dev/null \
        -e /dev/null \
        -R "affinity[core(10,same=socket,exclusive=(socket,alljobs)):membind=localonly:distribute=pack(socket=1)]" \
        OMP_NUM_THREADS=$t \
        ./main 256
done