wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v5.9.0.tar.gz 
tar -zxvf v5.9.0.tar.gz && cd SuiteSparse-5.9.0
make library -j4
sudo make install INSTALL=/usr/local/suitesparse
# make a copy to home as a cache
sudo cp -r /usr/local/suitesparse $HOME/suitesparse
