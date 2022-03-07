cd $HOME
for idx in {0..14}
do
  naunet example --select=$idx --dest=naunet_example$idx
done
