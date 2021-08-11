cd $HOME/naunet_example0/build
ctest -R "single|pymodule" -j4 --output-on-failure

cd $HOME/naunet_example1/build
ctest -R "single|pymodule" -j4 --output-on-failure

cd $HOME/naunet_example2/build
ctest -R "single|pymodule" -j4 --output-on-failure

cd $HOME/naunet_example3/build
ctest -R "single|pymodule" -j4 --output-on-failure
