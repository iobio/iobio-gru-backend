gru_version=${1}

apptainer build gru_${gru_version}.sif docker-daemon://quay.io/iobio/iobio-gru-backend:${gru_version}
