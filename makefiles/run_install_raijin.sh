#!/bin/sh

cd ..
installpath='/short/public/rcb547/apps/ga-aem/2.0.beta'

mkdir -p $installpath
mkdir -p $installpath/bin/raijin/intel

cp -pru docs $installpath
cp -pru examples $installpath
cp -pru python $installpath
cp -pru matlab $installpath
cp -pru bin/raijin/intel $installpath/bin/raijin

chmod -R go+rx $installpath/bin/raijin/intel/*.exe
chmod -R go+rx $installpath/matlab/bin/linux/*.so
chmod -R go+rx $installpath/python/gatdaem1d/*.so


