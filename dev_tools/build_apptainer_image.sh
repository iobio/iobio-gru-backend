gru_version_arg=${1}

if [ -z ${gru_version_arg} ]
then
    gru_version=$(git describe --tags)
else
    gru_version=$gru_version_arg
fi

echo "Building version $gru_version"

apptainer build gru_${gru_version}.sif docker-daemon://iobio/iobio-gru-backend:${gru_version}
