modes = [\
'hillshade',\
'slope',\
'aspect',\
'color-relief',\
'TRI',\
'TPI',\
'roughness']

gtiff_deflate = '\
-of GTiff \
-co COMPRESS=DEFLATE \
-co BLOCKXSIZE=512 \
-co BLOCKYSIZE=512 \
-oo NUM_THREADS=ALL_CPUS \
-co NUM_THREADS=ALL_CPUS \
--config GDAL_CACHEMAX 500'

color_table = '\
0 0 255 0/n\
64 252 141 89/n\
128 255,255,191/n\
192 145 207 96/n\
255 255 0 0'

color_table_file = '/tmp/temp1/color_table.txt'
with open(color_table_file, 'w+') as ctf:
    ctf.append(color_table)

for mode in modes:
    if mode == 'hillshade':
        mode_param = '\
        -z 3 \
        -s 1 \
    	-az 315 \
    	-alt 45 \
    	-multidirectional'
    elif mode == 'slope':
        mode_param = '\
        -s 1'
    elif mode == 'color-relief':
        mode_param = '\
        {} \
        -nearest_color_entry \
        '.format(color_table)

    cmd = '{} {} {} {} \
	-compute_edges \
	-alg Horn \
	-b 1 \
    {} \
    '.format(gdaldem, mode, input_dem, output_map, gtiff_deflate)
