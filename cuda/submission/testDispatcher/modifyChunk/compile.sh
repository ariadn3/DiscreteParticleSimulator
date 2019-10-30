/usr/local/cuda/bin/nvcc -o withfMad source/*.cu -fmad=true
/usr/local/cuda/bin/nvcc -o withoutfMad source/*.cu -fmad=false
