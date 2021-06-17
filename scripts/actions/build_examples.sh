cd $HOME/naunet_example0
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
      -DMAKE_PYTHON=ON
      # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
cmake --build build
# cd build && make install

cd $HOME/naunet_example1
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
      -DMAKE_PYTHON=ON
      # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
cmake --build build
# cd build && make install

cd $HOME/naunet_example2
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
      -DCMAKE_CUDA_ARCHITECTURES=61 \
      -DMAKE_PYTHON=ON
      # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
cmake --build build
# cd build && make install

cd $HOME/naunet_example3
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
      -DBOOST_ROOT=/usr/local \
      -DMAKE_PYTHON=ON
      # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
cmake --build build
# cd build && make install
