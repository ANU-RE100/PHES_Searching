#!/bin/bash
# ./write_to_postgres.sh path_to_output_folder postgres_password re100_username PHES_type [country]

if [ $# -eq 4 ]; then
    database="global_"${4,,}
elif [ $# -eq 5 ]; then
    database=${5,,}"_"${4,,}
else
    echo "Wrong number of inputs, input format is:"
    echo "./write_to_postgres.sh path_to_output_folder postgres_password re100_username PHES_type [country]"
    exit
fi

echo "Remember to connect to Global Protect"
ssh -M -S ~/my-ctrl-socket -fNT -L 5432:localhost:5432 $3@re100.anu.edu.au

# Check if database exists
if PGPASSWORD=$2 psql -U postgres --host localhost --port 5432 -lqt | cut -d \| -f 1 | grep -qw $database; then
    PGPASSWORD=$2 psql -U postgres --host localhost --port 5432 -d $database -c "DROP SCHEMA public CASCADE; CREATE SCHEMA public;"
else
    PGPASSWORD=$2 createdb $database -U postgres --host localhost --port 5432
fi

PGPASSWORD=$2 psql -U postgres --host localhost --port 5432 -d $database -c "CREATE EXTENSION postgis;"

# Copy files
cd $1"/kmls"
echo "Copying files to PostgreSQL"
for FILE in *; do 
    echo $FILE; 
    ogr2ogr -f "PostgreSQL" PG:"host=localhost user=postgres dbname=$database password=$2" $FILE;
done
cd -

ssh -S ~/my-ctrl-socket -O exit $3@re100.anu.edu.au