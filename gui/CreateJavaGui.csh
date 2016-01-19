#!/bin/csh

make CFLAGS="-shared -fPIC -O2" FFLAGS="-shared -fPIC -ffixed-line-length-132" javalib -f HathorMakeFile

make clean -f HathorMakeFile

gfortran dummy.f -o dummy.exe
set gfortranlib=`ldd dummy.exe | grep libgfortran | awk '{print $1}'`
rm -f dummy.exe
 
set dir=`pwd`
if ( $?LHAPDF ) then
   set lhapdf=$LHAPDF   
else
   set lhapdf=$HOME/local/lhapdf
endif

set lhalibs=`$lhapdf/bin/lhapdf-config --libdir`
set lhapath=`$lhapdf/bin/lhapdf-config --datadir`
set ldpreload=libLHAPDF.so:$gfortranlib


rm -f xhathor
cat > xhathor <<EOF
#!/bin/csh
if !(\$?HATHORPATH) then 
   setenv HATHORPATH $dir
endif
if !(\$?LHAPATH) then 
   setenv LHAPATH $lhapath
endif
setenv LD_LIBRARY_PATH $lhalibs/:\$HATHORPATH
setenv LD_PRELOAD $ldpreload
java -jar \$HATHORPATH/JHathor.jar 
EOF

chmod 700 xhathor

