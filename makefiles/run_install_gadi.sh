#!/bin/sh

cd ..
installpath='/scratch/cr78/apps/ga-aem/2.0.beta'

mkdir -p $installpath
mkdir -p $installpath/bin/gadi/intel

cp -pru docs $installpath
cp -pru examples $installpath
cp -pru python $installpath
cp -pru matlab $installpath
cp -pru bin/gadi/intel $installpath/bin/gadi

chmod -R go+rx $installpath/bin/gadi/intel/*.exe
chmod -R go+rx $installpath/matlab/bin/linux/*.so
chmod -R go+rx $installpath/python/gatdaem1d/*.so


