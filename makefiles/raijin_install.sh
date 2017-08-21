#!/bin/sh

installpath='/short/public/rcb547/apps/ga-aem/dev'
cd ..

mkdir -p $installpath
mkdir -p $installpath/bin/raijin/intel

cp -pru docs $installpath/docs
cp -pru examples $installpath/examples
cp -pru python $installpath/python
cp -pru matlab $installpath/matlab
cp -pru bin/raijin/intel $installpath/bin/raijin

chmod -R go+rx $installpath/bin/raijin/intel/*.exe
chmod -R go+rx $installpath/matlab/bin/linux/*.so
chmod -R go+rx $installpath/python/gatdaem1d/*.so


