#!/bin/bash

git clone https://github.com/atks/vt
cd vt
git checkout db574945014233a797eeafedad6e5f7e1ad2297a
make -j4
strip vt
