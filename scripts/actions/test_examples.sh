for idx in {0..18}
do
  if [[ "2 6 10 14" =~ (^|[[:space:]])"$idx"($|[[:space:]]) ]]; then
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      cd $HOME/naunet_example$idx/build
      ctest -R "single|pymodule" -j4 -V --output-on-failure
    fi
  else
    cd $HOME/naunet_example$idx/build
    ctest -R "single|pymodule" -j4 -V --output-on-failure
  fi
done
