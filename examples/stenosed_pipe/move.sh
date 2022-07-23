cd input
for file in *; do
  mv "$file" ../input1/"${file}.csv"
done
