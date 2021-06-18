ECHO Setting up QGIS DEV ENV x64_VS2019
 
$env:PYTHONPATH=''
$env:VS17COMNTOOLS="C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\Common7\Tools" 
$env:OSGEO4W_ROOT = 'c:\Program Files\QGIS 2.18\'

cmd.exe /c ${env:OSGEO4W_ROOT}\bin\o4w_env.bat
 
$env:QMAKESPEC="win64-msvc2019"
$env:PATH="${env:OSGEO4W_ROOT}\bin\;${env:OSGEO4W_ROOT}\apps\qgis\bin;${env:PATH}"
 
$env:INCLUDE="${env:INCLUDE};${env:OSGEO4W_ROOT}\include;${env:OSGEO4W_ROOT}\apps\qgis\include"
$env:LIB="${env:LIB};${env:OSGEO4W_ROOT}\lib;${env:OSGEO4W_ROOT}\apps\qgis\lib"
 
$env:PATH="${env:OSGEO4W_ROOT}\bin;${env:SYSTEMROOT}\System32;${env:SYSTEMROOT};${env:SYSTEMROOT}\System32\wbem;${env:PATH}"

