#!/bin/bash

old_version=${1}
new_version=${2}

function check_dir {
        if [ ! -d "${1}" ]; then
                echo "${1} does not exist"
                exit 1
        fi
}

old_dir=data/gru_data_${old_version}
check_dir ${old_dir}

update_dir=updates/updates_${old_version}_to_${new_version}
check_dir ${update_dir}

new_dir=data/gru_data_${new_version}
mkdir -p ${new_dir}

cp -alf ${old_dir}/* ${new_dir}/

chmod -R +w ${new_dir}

while read file; do
        rm "${new_dir}/${file}"
done < ${update_dir}/FILES_TO_DELETE

cp -alf ${update_dir}/* ${new_dir}/

rm -f "${new_dir}/FILES_TO_DELETE"

chmod -R -w ${new_dir}

echo "Results:"

diff -qr ${old_dir} ${new_dir}
