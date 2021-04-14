wget https://github.com/LLNL/sundials/releases/download/v5.7.0/sundials-5.7.0.tar.gz
tar -zxvf sundials-5.7.0.tar.gz
mkdir build-sundials && cd build-sundials
cmake ../sundials-5.7.0 -DCMAKE_INSTALL_PREFIX=/usr/local/sundials \
    -DENABLE_CUDA=ON -DCMAKE_CUDA_COMPILER=nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=61 -DSUNDIALS_INDEX_SIZE=32 \
    -DENABLE_KLU=ON -DKLU_INCLUDE_DIR=/usr/local/suitesparse/include \
    -DKLU_LIBRARY_DIR=/usr/local/suitesparse/lib -DEXAMPLES_ENABLE_C=OFF \
    -DEXAMPLES_ENABLE_CUDA=OFF -DEXAMPLES_INSTALL=OFF
make -j4 && sudo make install
# make a copy to home as a cache
sudo cp -r /usr/local/sundials $HOME/sundials
