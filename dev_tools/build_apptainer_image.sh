gru_version=${1}

apptainer build gru_${gru_version}.sif docker-daemon://iobio/iobio-gru-backend:${gru_version}
