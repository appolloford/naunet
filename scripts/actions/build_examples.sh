for idx in 0 1 4 5 8 9 12 13
do
  cd $HOME/naunet_example$idx
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
          -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
          -DMAKE_PYTHON=ON -DMAKE_TEST=ON
          # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
          -DMAKE_PYTHON=ON -DMAKE_TEST=ON
  fi
  cmake --build build
  # cd build && make install
done

for idx in 2 6 10
do
  cd $HOME/naunet_example$idx
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
          -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
          -DCMAKE_CUDA_ARCHITECTURES=61 \
          -DMAKE_PYTHON=ON -DMAKE_TEST=ON
          # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
    cmake --build build
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    echo "No CUDA driver! Skip! "
  fi
  # cd build && make install
done

for idx in 3 7 11 14
do
  cd $HOME/naunet_example$idx
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
          -DBOOST_ROOT=/usr/local \
          -DMAKE_PYTHON=ON -DMAKE_TEST=ON
          # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
          -DMAKE_PYTHON=ON -DMAKE_TEST=ON
  fi
  cmake --build build
  # cd build && make install
done 
