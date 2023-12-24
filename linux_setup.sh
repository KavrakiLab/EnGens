#!/bin/sh
eval "$(conda shell.bash hook)"
#conda activate <env-name>
conda activate latest

echo $(python --version)

pip install -i https://test.pypi.org/simple/ engens
git clone https://github.com/hsidky/hde.git

cp hde.patch ./hde/hde.patch
cd hde
git apply hde.patch
cd ..

pip install ./hde

#wget https://yanglab.nankai.edu.cn/mTM-align/version/mTM-align.tar.bz2
# yanglab temporarily unavailable!
cp dependencies/mTM-align.tar.bz2 ./mTM-align.tar.bz2
tar -xvf mTM-align.tar.bz2
cp mTM-align/src/mTM-align ${CONDA_PREFIX}/bin/mTM-align
rm mTM-align.tar.bz2

pypatch apply ./pdbfixer.patch pdbfixer
