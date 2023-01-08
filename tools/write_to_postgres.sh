#!/bin/bash
# ./write_to_postgres.sh path_to_output_folder database postgres_password

echo "Have you remembered to run \`CREATE EXTENSION postgis;\` in the data base and clear the data base?"
echo "If no type \"n\" then ENTER, otherwise press any key"
read x
if [ "$x" = "n" ]; then
    exit
fi
cd $1"/kmls"
echo "Copying files to PostgreSQL"
for FILE in *; do 
    echo $FILE; 
    ogr2ogr -f "PostgreSQL" PG:"host=localhost user=postgres dbname=$2 password=$3" $FILE; 
done
cd -
