for i in {0..504}; do
     file_name="${i}_hits.csv"
        [ ! -f $file_name ] && echo "$file_name"
done
