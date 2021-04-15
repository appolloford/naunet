cd $HOME/naunet_example0
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials
cmake --build build

cd $HOME/naunet_example1
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials
cmake --build build

cd $HOME/naunet_example2
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
      -DCMAKE_CUDA_ARCHITECTURES=61
cmake --build build

cd $HOME/naunet_example3
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DBOOST_ROOT=/usr/local
cmake --build build
