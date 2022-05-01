cd $HOME
for idx in {0..18}
do
  naunet example --select=$idx --path=naunet_example$idx
done
