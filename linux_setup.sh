#!/bin/sh
eval "$(conda shell.bash hook)"
#conda activate <env-name>
conda activate engens-conda-env

echo $(python --version)

pip install -i https://test.pypi.org/simple/ engens
git clone https://github.com/hsidky/hde.git

cp hde.patch ./hde/hde.patch
cd hde
git apply hde.patch
cd ..

pip install ./hde

