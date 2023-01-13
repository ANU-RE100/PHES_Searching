#!/bin/bash
# ./write_to_postgres.sh path_to_output_folder database postgres_password

# Check if database exists
if PGPASSWORD=$3 psql -U postgres -lqt | cut -d \| -f 1 | grep -qw $2; then
    PGPASSWORD=$3 psql -U postgres -d $2 -c "DROP SCHEMA public CASCADE; CREATE SCHEMA public;"
else
    PGPASSWORD=$3 createdb $2 -U postgres
fi

PGPASSWORD=$3 psql -U postgres -d $2 -c "CREATE EXTENSION postgis;"

# Copy files
cd $1"/kmls"
echo "Copying files to PostgreSQL"
for FILE in *; do 
    echo $FILE; 
    ogr2ogr -f "PostgreSQL" PG:"host=localhost user=postgres dbname=$2 password=$3" $FILE;
done
cd -
