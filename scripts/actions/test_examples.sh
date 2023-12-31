for idx in {0..18}
do
  if [[ "2 6 10 14" =~ (^|[[:space:]])"$idx"($|[[:space:]]) ]]; then
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      ctest -R "single|pymodule" -j4 -V --output-on-failure \
            --test-dir $HOME/naunet_example$idx/build
    fi
  else
    ctest -R "single|pymodule" -j4 -V --output-on-failure \
          --test-dir $HOME/naunet_example$idx/build
  fi
done
