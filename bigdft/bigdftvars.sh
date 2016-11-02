echo "While sourcing this script thethe following environment variable will be set:"
for file in futile_environment.sh bigdft_environment.sh
do
sh $file
eval `sh $file`
done
echo '...done.'
