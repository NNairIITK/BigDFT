echo "While sourcing this script the following environment variable will be set:"
for file in futile_environment.sh bigdft_environment.sh
do
sh $file
eval `sh $file`
done
$(return >/dev/null 2>&1)
if [ "$?" -eq "0" ]
then
    echo "...done."
else
    echo "ERROR: This script is NOT sourced! Execute 'source" $0 "' instead"
fi