import sys
import zipfile          # file zip managment.
tempFileZip='archivoCalibration.zip'

try:
# Se Lee el zip descargado
	archivo_zip = zipfile.ZipFile(tempFileZip, 'r')
	print '0'
	if(len(sys.argv) >1):
		print sys.argv[1]
	sys.exit(0)
except:
	print '1'
	if(len(sys.argv) >1):
		print sys.argv[1]
	sys.exit(1)
