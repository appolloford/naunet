for idx in {0..14}
do
  if [[ "2 6 10" =~ (^|[[:space:]])"$idx"($|[[:space:]]) ]]; then
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      cd $HOME/naunet_example$idx/build
      ctest -R "single|pymodule" -j4 --output-on-failure
    fi
  elif [[ "7" =~ (^|[[:space:]])"$idx"($|[[:space:]]) ]]; then
    # TODO: not sure why pymodule test freezes
    # ctest -R "single|pymodule" -j4 --output-on-failure
    cd $HOME/naunet_example$idx/build
    ctest -R "single" -j4 --output-on-failure
  else
    cd $HOME/naunet_example$idx/build
    ctest -R "single|pymodule" -j4 --output-on-failure
  fi
done
