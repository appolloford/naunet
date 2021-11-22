cd $HOME/naunet_example0
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
        -DMAKE_PYTHON=ON
        # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
elif [[ "$OSTYPE" == "darwin"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DMAKE_PYTHON=ON
fi
cmake --build build
# cd build && make install

cd $HOME/naunet_example1
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
        -DMAKE_PYTHON=ON
        # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
elif [[ "$OSTYPE" == "darwin"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DMAKE_PYTHON=ON
fi
cmake --build build
# cd build && make install

cd $HOME/naunet_example2
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
        -DCMAKE_CUDA_ARCHITECTURES=61 \
        -DMAKE_PYTHON=ON
        # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
  cmake --build build
elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo "No CUDA driver! Skip! "
fi
# cd build && make install

cd $HOME/naunet_example3
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DBOOST_ROOT=/usr/local \
        -DMAKE_PYTHON=ON
        # -DCMAKE_INSTALL_PREFIX=./ # not work on github workflow
elif [[ "$OSTYPE" == "darwin"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DMAKE_PYTHON=ON
fi
cmake --build build
# cd build && make install

cd $HOME/naunet_example4
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
        -DMAKE_PYTHON=ON
elif [[ "$OSTYPE" == "darwin"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DMAKE_PYTHON=ON
fi
cmake --build build

cd $HOME/naunet_example5
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
        -DMAKE_PYTHON=ON
elif [[ "$OSTYPE" == "darwin"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DMAKE_PYTHON=ON
fi
cmake --build build

cd $HOME/naunet_example6
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials \
        -DCMAKE_CUDA_ARCHITECTURES=61 \
        -DMAKE_PYTHON=ON
  cmake --build build  
elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo "No CUDA driver! Skip! "
fi

cd $HOME/naunet_example7
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release \
        -DBOOST_ROOT=/usr/local \
        -DMAKE_PYTHON=ON
elif [[ "$OSTYPE" == "darwin"* ]]; then
  cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DMAKE_PYTHON=ON
fi
cmake --build build
