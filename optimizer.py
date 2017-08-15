import sys
import os
import glob
import time
import shutil
import random
import numpy
import operator
from optparse import OptionParser
from osgeo import gdal

platform = sys.platform
if sys.platform.startswith('win'):
    platform = 'win'
    projmon = '<job dir> <job,aoi> <username>'
elif sys.platform.startswith('linux'):
    platform = 'linux'
    projmon = ''
else:
    print '\nCan not determine system platform.'
    sys.exit()

usage = 'usage: optimizer_dev <input> <output> {} [options]\n \
\n <input> <output> can take one of the following forms:\n \
\n  <input file> <output dir> \
\n  <input dir> <output dir> \
\n  <space delimited text file listing input dirs/files and output dirs>'.format(projmon)

parser = OptionParser(usage=usage)
parser.add_option("--it",dest='input_type',default='tif',\
help='Default=tif - Input any GDAL readable file extension. \
e.g. tif,pix,dem,bil,adf,vrt,ecw,sid,img... www.gdal.org/formats_list.html')
parser.add_option("--dt",dest='datatype',default='1',\
help='Default=1 - Dataset type. 1 = 8 bit imagery (tif-slab), \
2 = high bit depth terrain or bathy, 3 = cache tile (tif-cache)')
parser.add_option("--ts",dest='xy_tilesize',default='0,0',\
help='Default=None - Comma delimited pixel and line dimensions of output tile.')
parser.add_option("--pb",dest='pixel_buffer',type='int',default=0,\
help='Default=0 - Buffer distance in pixels for output tile.')
parser.add_option("--cmp",dest='compression',type='int',default=85,\
help='Default=85 - JPEG compression quality for 8 bit imagery.')
parser.add_option("--emp",dest='empty',type='float',default=0.7,\
help='Default=0.7 - Maximum mean value for empty output tiles, losers are skipped. \
Set to -32767 with --dt=2.')
parser.add_option("--rz",dest='randomize',action='store_true',default=False,\
help='Default=False - Randomize the order in which input tiles are processed. \
Not available with --vrt.')
parser.add_option("--vrt",dest='vrt',action='store_true',default=False,\
help='Default=False - Assemble all tiles in folder to a VRT for re-tiling. \
Only functions while input path is a dir or textfile list of dirs.')
parser.add_option("--vrt_proj",dest='vrt_proj',default='',\
help='Default=None - VRT SRS code for re-projecting. Use wm for Web Mecator. e.g. --vrt_proj EPSG:4326')
parser.add_option("--nc",dest='nocopy',action='store_true',default=False,\
help='Default=False - No copying of input to local dir. Set to True with --vrt.')
parser.add_option("--id",dest='script_id',type='int',default=1,\
help='Default=1 - ID of the script. Windows only.')
parser.add_option("--ns",dest='n_scripts',default=1,type="int",\
help="Default=1 - Number of scripts to launch. Windows only.")
parser.add_option("--opts",dest='gdal_options',default='',\
help='Default=None - gdal_translate options in double quotes. \
e.g. "-b 1 -b 2 -b 3 -mask 4" to create RGB image with external mask from internal band 4.')
parser.add_option("--upload",dest='gcs',default='',\
help='Default=None - e.g. existing_bucket_guid/new_folder Upload optimized data to Google Cloud Stroage bucket. \
Windows platform must have Google Cloud SDK installed and authenticated.')

(options,args) = parser.parse_args()

if len(args) == 0:
    print '\n'
    parser.print_help()
    sys.exit()

inp = args[0]
s_id = options.script_id
ns = options.n_scripts
it = options.input_type
compression = options.compression
x = int(options.xy_tilesize.split(',')[0])
y = int(options.xy_tilesize.split(',')[1])
pb = options.pixel_buffer
randomize = options.randomize
vrt = options.vrt
nocopy = options.nocopy
empty = options.empty
opts = options.gdal_options
datatype = options.datatype
t_epsg = options.vrt_proj
gcs = options.gcs

if platform == 'win':
    import i3functions_osgeo4w as i
    from i3functions import spatial as s
    from i3functions import common as c

    gdaladdo = i.gdaladdo
    gdal_translate = i.gdal_translate
    gdalbuildvrt = i.gdalbuildvrt
    gdalwarp = i.gdalwarp

    if len(args) != 5 and not inp.endswith('.txt'):
        print '\n'
        parser.print_help()
        sys.exit()
		
	elif len(args) != 4 and inp.endswith('.txt'):
        print '\n'
        parser.print_help()
        sys.exit()

    stime = i.printStartTime()

    if inp.endswith('.txt'):
        pass
        workdir = args[1]
        name = args[2]
        username = args[3]
    else:
        outdir = args[1]
        workdir = args[2]
        name = args[3]
        username = args[4]

    if not ',' in name:
        print '\nProject name and AOI name should be specified as comma separated values.\n'
        sys.exit()
    else:
        jobno = name.split(',')[0]
        name = name.replace(',','+')

    fullname = c.getFullName(username)
    if not fullname:
        print 'Nice try. "{}" user name does not exist. Please use a valid user name.\n'.format(username)
        sys.exit()
    else:
        print '\nUser:',fullname

    if ns > 1:
        for sid in range(2,ns+1):
            launch_cmd = '%s\\optimizer_dev.bat '.format(i.scripts) \
            + ' '.join(['"{}"'.format(a) \
            if not '"' in a else a for a in sys.argv[1:]]) \
            + ' --id {} --ns 1'.format(sid)
            os.system(r'start /SEPARATE {}'.format(launch_cmd))
            print launch_cmd,'\n'

    activef,activef_h = c.writeLockFile(s_id,__file__,workdir,jobno,fullname)

    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    pointer_txt  = os.path.join(os.path.dirname(c.scripts),'ProjMon','{}.txt'\
    .format(name))
    pointer_h = open(pointer_txt,'w')
    workdir = c.getUNCPath(workdir)
    pointer_h.write('%s,%s,0,%s'%(workdir,0,fullname))
    pointer_h.close()

elif platform == 'linux':
    gdaladdo = '/bin/gdaladdo'
    gdal_translate = '/bin/gdal_translate'
    gdalbuildvrt = '/bin/gdalbuildvrt'
    gdalwarp = '/bin/gdalwarp'

    if (args) < 2 and not os.path.isfile(inp.endswith('.txt')):
        print '\n'
        parser.print_help()
        sys.exit()

    outdir = args[1]


def removeFile(infile):
    infiles = glob.glob(os.path.splitext(infile)[0]+'*.*')
    for item in infiles:
        try:
            os.remove(item)
        except Exception,e:
            'Could not remove {}'.format(item)
            print '\n',e

def getTileInfo(infile,temp,xpxsz,ypxsz,pxbuff):
    ds = gdal.Open(infile)
    width = ds.RasterXSize
    height = ds.RasterYSize
    bands = ds.RasterCount

    if width < xpxsz or width == xpxsz:
        xpxsz = width
        print '\nAdjusting user defined x to %s' %(width)
    if xpxsz == 0:
        num_xtile = 1
        last_x = width
    if xpxsz != 0:
        num_xtile = width/xpxsz
        remain_x = operator.mod(width,xpxsz)
        if remain_x > 0 and width > xpxsz:
            num_xtile = num_xtile + 1
            last_x = remain_x
        elif remain_x == width:
            last_x = width
        elif remain_x == 0 and width > xpxsz:
            last_x = xpxsz
        else:
            last_x = width

    if height < ypxsz or height == ypxsz:
        ypxsz = height
        print '\nAdjusting user defined y to %s' %(height)
    if ypxsz == 0:
        num_ytile = 1
        last_y = height
    if ypxsz != 0:
        num_ytile = height/ypxsz
        remain_y = operator.mod(height,ypxsz)
        if remain_y > 0 and height > ypxsz:
            num_ytile = num_ytile + 1
            last_y = remain_y
        elif remain_y == height:
            last_y = height
        elif remain_y == 0 and height > ypxsz:
            last_y = ypxsz
        else:
            last_y = height
    return width,height,num_xtile,num_ytile,last_x,last_y,xpxsz,ypxsz,bands

def getTileName(infile,num_xtile,num_ytile,col_id,row_id):
    if num_xtile > 1 and num_ytile > 1:
        if col_id < 10:
            col = "C0" + str(col_id)
        else:
            col = "C" + str(col_id)
        if row_id < 10:
            row = "R0" + str(row_id)
        else:
            row = "R" + str(row_id)
        outimage = os.path.basename(infile.split('.')[0]) + "_" + col + row
    elif num_xtile == 1 and num_ytile > 1:
        col = "C01"
        if row_id < 10:
            row = "R0" + str(row_id)
        else:
            row = "R" + str(row_id)
        outimage = os.path.basename(infile.split('.')[0]) + "_" + col + row
    elif num_xtile > 1 and num_ytile == 1:
        if col_id < 10:
            col = "C0" + str(col_id)
        else:
            col = "C" + str(col_id)
        row = "R01"
        outimage = os.path.basename(infile.split('.')[0]) + "_" + col + row
    else:
        outimage = os.path.basename(infile.split('.')[0])
    return outimage

def getTileCoords(xpxsz,ypxsz,pxbuff,col_id,row_id,num_xtile,num_ytile,last_x,last_y):
    if num_xtile == 1:
        xstart = 0
        xsize = last_x
    elif num_xtile > 1 and col_id == 1:
        xstart = 0
        xsize = xpxsz + pxbuff
    elif num_xtile > 1 and col_id == num_xtile:
        xstart = ((col_id-1)*xpxsz-pxbuff)
        xsize = last_x + pxbuff
    else:
        xstart = ((col_id-1)*xpxsz-pxbuff)
        xsize = xpxsz + pxbuff + pxbuff
    if num_ytile == 1:
        ystart = 0
        ysize = last_y
    elif num_ytile > 1 and row_id == 1:
        ystart = 0
        ysize = ypxsz + pxbuff
    elif num_ytile > 1 and row_id == num_ytile:
        ystart = ((row_id-1)*ypxsz-pxbuff)
        ysize = last_y + pxbuff
    else:
        ystart = ((row_id-1)*ypxsz-pxbuff)
        ysize = ypxsz + pxbuff + pxbuff
    return xstart,ystart,xsize,ysize

def checkMean(infile,xstart,ystart,xsize,ysize,maxval,num_xtile,num_ytile,datatype):
    if datatype == '2':
        maxval = -32767
    valid = False
    try:
        if num_xtile != 1 and num_ytile != 1:
            print '\nCheck mean...\n{}'.format(infile)
            ds = gdal.Open(infile)
            band = ds.GetRasterBand(1)
            print 'Srcwin: {} {} {} {}'.format(xstart,ystart,xsize,ysize)
            subset = band.ReadAsArray(xstart,ystart,xsize,ysize,1000,1000)
            stat = numpy.mean(subset)
            print 'Mean: {} > {}'.format(stat,maxval)
            if stat > maxval:
                valid = True
            ds = None
        else:
            valid = True
    except Exception,e:
        print '\n',e
    return valid

def validate(filename,datatype):
    class ValidateCloudOptimizedGeoTIFFException(Exception):
        pass

    try:
        # Check local GDAL version
        if int(gdal.VersionInfo('VERSION_NUM')) < 2020000:
            raise ValidateCloudOptimizedGeoTIFFException(
                "GDAL 2.2 or above required")

        gdal.PushErrorHandler()
        ds = gdal.Open(filename)
        gdal.PopErrorHandler()

        main_band = ds.GetRasterBand(1)
        block_size = main_band.GetBlockSize()

        if ds is None:
            raise ValidateCloudOptimizedGeoTIFFException(
                "Invalid file : %s" % gdal.GetLastErrorMsg())

        # Check Tiled GTiff
        if ds.GetDriver().ShortName != 'GTiff':
            raise ValidateCloudOptimizedGeoTIFFException(
                "The file is not a GeoTIFF")

        if block_size[0] == main_band.XSize:
            raise ValidateCloudOptimizedGeoTIFFException(
                "Full resolution image is not tiled")

        # Check IFD structure
        ifd_offset = [int(main_band.GetMetadataItem('IFD_OFFSET', 'TIFF'))]
        if ifd_offset[0] != 8:
            raise ValidateCloudOptimizedGeoTIFFException(
                "The offset of the main IFD should be 8. It is %d instead" %
                ifd_offset[0])

        if datatype != '3':

            # Check overviews
            ovr_count = main_band.GetOverviewCount()
            if filename + '.ovr' in ds.GetFileList():
                raise ValidateCloudOptimizedGeoTIFFException(
                    "Overviews should be internal")

            if main_band.XSize >= 512 or main_band.YSize >= 512:
                if ovr_count == 0:
                    raise ValidateCloudOptimizedGeoTIFFException(
                        "The file should have overviews")

            for i in range(ovr_count):
                # Check that overviews are by descending sizes
                ovr_band = ds.GetRasterBand(1).GetOverview(i)
                if i == 0:
                    if ovr_band.XSize > main_band.XSize or \
                       ovr_band.YSize > main_band.YSize:
                            raise ValidateCloudOptimizedGeoTIFFException(
                               "First overview has larger dimension than main band")
                else:
                    prev_ovr_band = ds.GetRasterBand(1).GetOverview(i-1)
                    if ovr_band.XSize > prev_ovr_band.XSize or \
                       ovr_band.YSize > prev_ovr_band.YSize:
                            raise ValidateCloudOptimizedGeoTIFFException(
                               "Overview of index %d has larger dimension than "
                               "overview of index %d" % (i, i-1))

                block_size = ovr_band.GetBlockSize()
                if block_size[0] == ovr_band.XSize:
                    raise ValidateCloudOptimizedGeoTIFFException(
                        "Overview of index %d is not tiled" % i)

                # Check that the IFD of descending overviews are sorted by increasing
                # offsets
                ifd_offset.append(int(ovr_band.GetMetadataItem('IFD_OFFSET', 'TIFF')))
                if ifd_offset[-1] < ifd_offset[-2]:
                    if i == 0:
                        raise ValidateCloudOptimizedGeoTIFFException(
                            "The offset of the IFD for overview of index %d is %d, "
                            "whereas it should be greater than the one of the main "
                            "image, which is at byte %d" %
                            (i, ifd_offset[-1], ifd_offset[-2]))
                    else:
                        raise ValidateCloudOptimizedGeoTIFFException(
                            "The offset of the IFD for overview of index %d is %d, "
                            "whereas it should be greater than the one of index %d, "
                            "which is at byte %d" %
                            (i, ifd_offset[-1], i-1, ifd_offset[-2]))

            # Check that the imagery starts by the smallest overview and ends with
            # the main resolution dataset
            data_offset = [int(main_band.GetMetadataItem('BLOCK_OFFSET_0_0', 'TIFF'))]
            for i in range(ovr_count):
                ovr_band = ds.GetRasterBand(1).GetOverview(i)
                data_offset.append(int(
                    ovr_band.GetMetadataItem('BLOCK_OFFSET_0_0', 'TIFF')))

            if data_offset[-1] < ifd_offset[-1]:
                if ovr_count > 0:
                    raise ValidateCloudOptimizedGeoTIFFException(
                        "The offset of the first block of the smallest overview "
                        "should be after its IFD")
                else:
                    raise ValidateCloudOptimizedGeoTIFFException(
                        "The offset of the first block of the image should "
                        "be after its IFD")
            for i in range(len(data_offset)-2, 0, -1):
                if data_offset[i] < data_offset[i+1]:
                    raise ValidateCloudOptimizedGeoTIFFException(
                        "The offset of the first block of overview of index %d should "
                        "be after the one of the overview of index %d" %
                        (i-1, i))
            if len(data_offset) >= 2 and data_offset[0] < data_offset[1]:
                raise ValidateCloudOptimizedGeoTIFFException(
                    "The offset of the first block of the main resolution image "
                    "should be after the one of the overview of index %d" %
                    (ovr_count - 1))

        msg = '\n{}: is a valid cloud optimized GeoTIFF!'.format(os.path.basename(filename))
        return True, msg
    except ValidateCloudOptimizedGeoTIFFException as e:
        msg = '\n{}: is NOT a valid cloud optimized GeoTIFF!\n{}\n'.format(os.path.basename(filename),str(e))
        return False, msg

def assembleVRT(folder,mos_vrt):
    vrt_cmd = '{} {} {}'.format(gdalbuildvrt,mos_vrt,folder)
    try:
        print '\n',vrt_cmd
        os.system(vrt_cmd)
    except Exception,e:
        print '\n',e

def warpVRT(vrt_file,warp_file,t_epsg,datatype):
    if datatype == '2':
        resample = 'near'
    else:
        resample = 'cubic'
    if t_epsg == 'wm':
        wkt = 'PROJCS["WGS_1984_Web_Mercator",GEOGCS["GCS_WGS_1984_Major_Auxiliary_Sphere",DATUM["WGS_1984_Major_Auxiliary_Sphere",SPHEROID["WGS_1984_Major_Auxiliary_Sphere",6378137.0,0.0]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Mercator_1SP"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0],EXTENSION["PROJ4","+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"],AUTHORITY["EPSG","900913"]]'
        wm_wkt = os.path.join(temp2,'wm.wkt')
        with open(wm_wkt, 'a+') as wf:
            wf.write('{}'.format(wkt))
        t_epsg = wm_wkt
    warp_cmd = '{} -of VRT -r {} -t_srs {} {} {}'.format(gdalwarp,resample,t_epsg,vrt_file,warp_file)
    try:
        print '\n',warp_cmd
        os.system(warp_cmd)
    except Exception,e:
        print '\n',e

def optimizeRaster(infile,outdir,outimage,temp_dir1,temp_dir2,num_xtile,num_ytile,xstart,ystart,xsize,ysize,compression,opts,bands,datatype):
    bn = os.path.splitext(os.path.basename(outimage))[0]
    wrk_tif = os.path.join(temp_dir2, bn +'.tif')
    if datatype != '3':
        deflate = os.path.join(temp_dir1, bn +'_deflate.tif')
        deflate_cmd = '{} {} -co BIGTIFF=YES -co TILED=YES -co COMPRESS=DEFLATE -oo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS --config GDAL_CACHEMAX 500 -srcwin {} {} {} {} "{}" "{}"'.format(gdal_translate,opts,xstart,ystart,xsize,ysize,infile,deflate)
        if datatype == '1':
            if bands > 3:
                photometric = 'MINISBLACK'
            elif bands == 1:
                photometric = ''
            else:
                photometric = 'YCBCR'
            gdaladdo_cmd = '{} -r average -oo NUM_THREADS=ALL_CPUS {} 2 4 8 16 32'.format(gdaladdo,deflate)
            gdal_param = '-of GTiff -co TILED=YES -co COMPRESS=JPEG -co JPEG_QUALITY={0} -co BLOCKXSIZE=512 -co BLOCKYSIZE=512 -co PHOTOMETRIC={1} -co COPY_SRC_OVERVIEWS=YES -oo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS --config JPEG_QUALITY_OVERVIEW {0} --config GDAL_TIFF_OVR_BLOCKSIZE 512 --config GDAL_CACHEMAX 500'.format(compression,photometric)
        elif datatype == '2':
            gdaladdo_cmd = '{} -r nearest -oo NUM_THREADS=ALL_CPUS {} 2 4 8 16 32'.format(gdaladdo,deflate)
            gdal_param = '-of GTiff -stats -co TILED=YES -co COMPRESS=LZW -co PREDICTOR=2 -co COPY_SRC_OVERVIEWS=YES -co BLOCKXSIZE=512 -co BLOCKYSIZE=512 -oo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS --config GDAL_TIFF_OVR_BLOCKSIZE 512 --config GDAL_CACHEMAX 500'
        translate_cmd = '{} {} "{}" "{}"'.format(gdal_translate,gdal_param,deflate,wrk_tif)
    else:
        deflate_cmd = ''
        gdaladdo_cmd = ''
        gdal_param = '-of GTiff -co TILED=YES -co COMPRESS=JPEG -co JPEG_QUALITY={} -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 -co PHOTOMETRIC=YCBCR -oo NUM_THREADS=ALL_CPUS -co NUM_THREADS=ALL_CPUS --config GDAL_CACHEMAX 500'.format(compression)
        translate_cmd = '{} {} "{}" "{}"'.format(gdal_translate,gdal_param,infile,wrk_tif)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    proc_tile = False
    if platform == 'win':
        if not num_xtile == 1 and not num_ytile == 1:
            wrk_tile = i.createWrkfile(os.path.join(outdir, outimage) + '.tif')
            if wrk_tile:
                proc_tile = True
        else:
            if not os.path.isfile(os.path.join(outdir, outimage) + '.tif'):
                if not os.path.isfile(os.path.join(outdir, outimage) + '.wrk'):
                    proc_tile = True
    elif platform == 'linux':
        if not os.path.isfile(os.path.join(outdir, outimage) + '.tif'):
            proc_tile = True
        else:
            print '\nOutput already exists. Skipping.'
    if proc_tile:
        try:
            if datatype != '3':
                print '\nCreating cloud optimized GeoTIFF...'
                print '\n',deflate_cmd
                os.system(deflate_cmd)
                print '\n',gdaladdo_cmd
                os.system(gdaladdo_cmd)
                print '\n',gdal_param
                os.system(translate_cmd)
            else:
                print '\nCreating cloud optimized GeoTIFF...'
                print '\n',gdal_param
                os.system(translate_cmd)
            valid, valid_msg = validate(wrk_tif,datatype)
            if not valid:
                print valid_msg,'\nHalting further processing! Check gdal_translate parameters!'
                sys.exit()
            else:
                print valid_msg
            if platform == 'win':
                i.copyFile(wrk_tif,outdir,exclude='/xf *.xml *.vrt')
                os.remove(wrk_tile)
            elif platform == 'linux':
                cp_cmd = 'cp {} {}'.format(wrk_tif.replace('.tif','.*'),outdir)
                print '\nCopying to destination...{}'.format(os.path.join(outdir,wrk_tif))
                os.system(cp_cmd)
        except Exception,e:
            print '\n',e

def uploadGCS(outdir, bucket_folder):
    upload_cmd = 'gsutil -m cp -r {}\*.tif* gs://{}/'.format(outdir, bucket_folder)
    try:
        print '\nUploading data to Google Cloud Storage...\n\n',upload_cmd
        os.system(upload_cmd)
    except Exception,e:
        print '\n',e

################################################################################

if os.path.isfile(inp) and not inp.endswith('.txt'):
    count = 0
    infolder = os.path.split(inp)[0]

    process = False
    if platform == 'win':
        comp_name,no_processors,temp1,temp2 = i.createTempdir('optimizer',s_id)
        stime_tile = time.time()
        wrk_file = os.path.join(workdir,'{}.wrk'.format(os.path.splitext(os.path.basename(inp))[0]))
        wrk_message = c.createWrkFileDirectly(wrk_file,id=s_id,count=1,total=1,indir=infolder,outdir=outdir,workdir=workdir,stime_tile=stime_tile)
        if wrk_message == 'Created':
            process = True
        activef,activef_h = c.writeLockFile(s_id,__file__,workdir,jobno,fullname)
    elif platform == 'linux':
        if not os.path.isfile(os.path.join(os.path.join(outdir,'{}.tif'.format(os.path.splitext(os.path.basename(inp))[0])))):
            process = True
            temp1 = '/tmp/temp1'
            temp2 = '/tmp/temp2'
            if not os.path.exists(temp1):
                os.makedirs(temp1)
            if not os.path.exists(temp2):
                os.makedirs(temp2)

    it = inp.split('.')[1]
    temp_inimage = os.path.join(temp1,os.path.basename(inp))

    if process:
        csv_table = ['Raster,Name']
        width,height,num_xtile,num_ytile,last_x,last_y,xpxsz,ypxsz,b = getTileInfo(inp,temp1,x,y,pb)
        if not nocopy:
            if it == 'adf':
                if platform == 'win':
                    i.copyFolder(infolder,temp1)
                elif platform == 'linux':
                    cp_cmd = 'cp {} {}'.format(os.path.join(infolder,'*'),temp1)
                    print 'Copying to scratch...{}'.format(temp1)
                    os.system(cp_cmd)
            elif it == 'vrt':
                temp_inimage = inp
                nocopy = True
            else:
                if platform == 'win':
                    i.copyFile(inp,temp1,exclude='/xf *.aux')
                elif platform == 'linux':
                    cp_cmd = 'cp {} {}'.format(inp,temp1)
                    print 'Copying to scratch...{}'.format(temp1)
                    os.system(cp_cmd)
        else:
            temp_inimage = inp

        for k in range(1,num_xtile+1):
            for j in range(1,num_ytile+1):
                outname = getTileName(temp_inimage,num_xtile,num_ytile,k,j)
                print '\nInput Image: {}\nImage Size: {}, {}\nOutput Image: {}.tif\nCols/Rows: {}x{}'\
                .format(os.path.basename(temp_inimage),width,height,outname,num_xtile,num_ytile)
                xstart,ystart,xsize,ysize = getTileCoords(xpxsz,ypxsz,pb,k,j,num_xtile,num_ytile,last_x,last_y)
                if checkMean(temp_inimage,xstart,ystart,xsize,ysize,empty,num_xtile,num_ytile,datatype):
                    optimizeRaster(temp_inimage,outdir,outname,temp1,temp2,num_xtile,num_ytile,xstart,ystart,xsize,ysize,compression,opts,b,datatype)
                    if k == 1 and j == 1:
                        if platform == 'win':
                            try:
                                image_size_gb = s.getUncompressedSize(temp_outimage)
                            except Exception,e:
                                print '\n',e
                                image_size_gb = 1
                    if gcs != '':
                        csv_table.append('/vsicurl/http://storage.googleapis.com/{}/{},{}'.format(gcs,os.path.basename(temp_outimage),os.path.basename(temp_outimage).replace('.tif','')))
                else:
                    image_size_gb = 1
                    print 'Output tile contains only NoData and will be skipped.'
                j+=1
                continue
            k+=1
            continue

        if gcs != '':
            uploadGCS(outdir,gcs)

        removeFile(os.path.join(temp1,'*'))

        if platform == 'win':
            pointer_h = open(pointer_txt,'w')
            pointer_workdir = c.getUNCPath(workdir)
            pointer_h.write('%s,%s,%s,%s'%(pointer_workdir,1,image_size_gb,fullname))
            pointer_h.close()
            ptime_tile = time.time() - stime_tile
            c.createDoneFile(wrk_file,id=s_id,ptime_tile=ptime_tile,wrk_file_h = '',tile_size_gb=image_size_gb,jobno=jobno)

    if gcs != '':
        print 'Writing table...'
        logfile = os.path.join(temp2,'_table.txt')
        with open(logfile, 'a+') as lf:
            for item in csv_table:
                lf.write('{}\n'.format(item))

        if platform == 'win':
            i.copyFile(logfile,outdir,exclude='/xf *.tif* *.vrt')
        elif platform == 'linux':
            cp_cmd = 'cp {} {}'.format(logfile,outdir)
            print 'Copying table...'
            os.system(cp_cmd)

    removeFile(os.path.join(temp2,'*'))

elif os.path.isdir(inp):
    count = 0
    infolder = os.path.split(inp)[0]

    if platform == 'win':
        comp_name,no_processors,temp1,temp2 = i.createTempdir('optimizer',s_id)
    elif platform == 'linux':
        temp1 = '/tmp/temp1'
        temp2 = '/tmp/temp2'
        if not os.path.exists(temp1):
            os.makedirs(temp1)
        if not os.path.exists(temp2):
            os.makedirs(temp2)

    pattern = os.path.join(inp,'*.{}'.format(it))
    inimages = []
    if vrt:
        names = glob.glob(pattern)
        wrk_vrt = os.path.join(temp1, os.path.basename(names[0]).replace('.{}'.format(it),'.vrt'))
        inimages.append(wrk_vrt)
        nocopy = True
    else:
        inimages = glob.glob(pattern)
        if randomize:
            random.shuffle(inimages)
    len_inimages = len(inimages)

    if platform == 'win':
        pointer_txt  = os.path.join(os.path.dirname(c.scripts),'ProjMon','%s.txt'%name)
        pointer_h = open(pointer_txt,'w')
        workdir = c.getUNCPath(workdir)
        pointer_h.write('%s,%s,0,%s'%(workdir,len_inimages,fullname))
        pointer_h.close()

    for image in inimages:
        count +=1
        process = False
        if platform == 'win':
            stime_tile = time.time()
            wrk_file = os.path.join(workdir,'%s.wrk' %os.path.splitext(os.path.basename(image))[0])
            wrk_message = c.createWrkFileDirectly(wrk_file,id=s_id,count=count,total=len_inimages,indir=infolder,outdir=outdir,workdir=workdir,stime_tile=stime_tile)
            if wrk_message == 'Created':
                process = True
            activef,activef_h = c.writeLockFile(s_id,__file__,workdir,jobno,fullname)
        elif platform == 'linux':
            if not os.path.isfile(os.path.join(os.path.join(outdir,'{}*.tif'.format(os.path.splitext(os.path.basename(inp))[0])))):
                process = True

        temp_inimage = os.path.join(temp1,os.path.basename(image))

        if process:
            csv_table = ['Raster,Name']
            if not nocopy or not vrt:
                if platform == 'win':
                    i.copyFile(image,temp1,exclude='/xf *.aux')
                elif platform == 'linux':
                    cp_cmd = 'cp {} {}'.format(image,temp1)
                    print 'Copying to scratch...{}'.format(temp1)
                    os.system(cp_cmd)
            else:
                temp_inimage = image
            if vrt:
                assembleVRT(pattern,temp_inimage)
                if t_epsg != '':
                    wrk_warp = os.path.join(temp2,os.path.basename(temp_inimage))
                    warpVRT(temp_inimage,wrk_warp,t_epsg,datatype)
                    temp_inimage = wrk_warp
            width,height,num_xtile,num_ytile,last_x,last_y,xpxsz,ypxsz,b = getTileInfo(temp_inimage,temp1,x,y,pb)
            for k in range(1,num_xtile+1):
                for j in range(1,num_ytile+1):
                    outname = getTileName(temp_inimage,num_xtile,num_ytile,k,j)
                    print '\nInput Image: {}\nImage Size: {}, {}\nOutput Image: {}.tif\nCols/Rows: {}x{}'.format(os.path.basename(temp_inimage),width,height,outname,num_xtile,num_ytile)
                    xstart,ystart,xsize,ysize = getTileCoords(xpxsz,ypxsz,pb,k,j,num_xtile,num_ytile,last_x,last_y)
                    if checkMean(temp_inimage,xstart,ystart,xsize,ysize,empty,num_xtile,num_ytile,datatype):
						optimizeRaster(temp_inimage,outdir,outname,temp1,temp2,num_xtile,num_ytile,xstart,ystart,xsize,ysize,compression,opts,b,datatype)
                        if k == 1 and j == 1:
                            if platform == 'win':
                                try:
                                    image_size_gb = s.getUncompressedSize(temp_outimage)
                                except Exception,e:
                                    print '\n',e
                                    image_size_gb = 1
                        if gcs != '':
                            csv_table.append('/vsicurl/http://storage.googleapis.com/{}/{},{}'.format(gcs,os.path.basename(temp_outimage),os.path.basename(temp_outimage).replace('.tif','')))
                    else:
                        image_size_gb = 1
                        print 'Output tile contains only NoData and will be skipped.'
                    j+=1
                    continue
                k+=1
                continue

            if gcs != '':
                uploadGCS(outdir,gcs)

            removeFile(os.path.join(temp1,'*.*'))

            if platform == 'win':
                pointer_h = open(pointer_txt,'w')
                pointer_workdir = c.getUNCPath(workdir)
                pointer_h.write('%s,%s,%s,%s'%(pointer_workdir,1,image_size_gb,fullname))
                pointer_h.close()
                ptime_tile = time.time() - stime_tile
                c.createDoneFile(wrk_file,id=s_id,ptime_tile=ptime_tile,wrk_file_h = '',tile_size_gb=image_size_gb,jobno=jobno)

        if gcs != '':
            print 'Writing table...'
            logfile = os.path.join(temp2,'_table.txt')
            with open(logfile, 'a+') as lf:
                for item in csv_table:
                    lf.write('{}\n'.format(item))

            if platform == 'win':
                i.copyFile(logfile,outdir,exclude='/xf *.tif* *.vrt')
            elif platform == 'linux':
                cp_cmd = 'cp {} {}'.format(logfile,outdir)
                print 'Copying table...'
                os.system(cp_cmd)

        removeFile(os.path.join(temp2,'*.*'))

elif os.path.isfile(inp) and inp.endswith('.txt'):
    count = 0
    infolder = os.path.split(inp)[0]

    if platform == 'win':
        comp_name,no_processors,temp1,temp2 = i.createTempdir('optimizer',s_id)
    elif platform == 'linux':
        temp1 = '/tmp/temp1'
        temp2 = '/tmp/temp2'
        if not os.path.exists(temp1):
            os.makedirs(temp1)
        if not os.path.exists(temp2):
            os.makedirs(temp2)

    lines = open(inp).readlines()
    for line in lines:
        line = line.strip()
        if line:
            indir = line.split()[0]
            outdir = line.split()[1]
            inimages = []
            if os.path.isfile(indir):
                inimages.append(indir)
                total = 1
            elif os.path.isdir(indir):
                pattern = os.path.join(indir,'*.{}'.format(it))
                if vrt:
                    names = glob.glob(pattern)
                    wrk_vrt = os.path.join(temp1, os.path.basename(names[0]).replace('.{}'.format(it),'.vrt'))
                    inimages.append(wrk_vrt)
                    nocopy = True
                else:
                    inimages = glob.glob(pattern)
                    if randomize:
                        random.shuffle(inimages)

            print '\nInput: {}\nOutput: {}\n'.format(indir,outdir)

            len_inimages = len(inimages)

            if platform == 'win':
                pointer_txt  = os.path.join(os.path.dirname(c.scripts),'ProjMon','%s.txt'%name)
                pointer_h = open(pointer_txt,'w')
                workdir = c.getUNCPath(workdir)
                pointer_h.write('%s,%s,0,%s'%(workdir,len_inimages,fullname))
                pointer_h.close()

            for image in inimages:
                count +=1
                process = False
                if platform == 'win':
                    stime_tile = time.time()
                    wrk_file = os.path.join(workdir,'%s.wrk' %os.path.splitext(os.path.basename(image))[0])
                    wrk_message = c.createWrkFileDirectly(wrk_file,id=s_id,count=count,total=len_inimages,indir=indir,outdir=outdir,workdir=workdir,stime_tile=stime_tile)
                    if wrk_message == 'Created':
                        process = True
                    activef,activef_h = c.writeLockFile(s_id,__file__,workdir,jobno,fullname)
                elif platform == 'linux':
                    if not os.path.isfile(os.path.join(os.path.join(outdir,'{}*.tif'.format(os.path.splitext(os.path.basename(inp))[0])))):
                        process = True

                temp_inimage = os.path.join(temp1,os.path.basename(image))

                if process:
                    csv_table = ['Raster,Name']
                    if not nocopy or not vrt:
                        if platform == 'win':
                            i.copyFile(image,temp1,exclude='/xf *.aux')
                        elif platform == 'linux':
                            cp_cmd = 'cp {} {}'.format(image,temp1)
                            print 'Copying to scratch...{}'.format(temp1)
                            os.system(cp_cmd)
                    else:
                        temp_inimage = image
                    if vrt:
                        assembleVRT(pattern,temp_inimage)
                        if t_epsg != '':
                            wrk_warp = os.path.join(temp2,os.path.basename(temp_inimage))
                            warpVRT(temp_inimage,wrk_warp,t_epsg,datatype)
                            temp_inimage = wrk_warp
                    width,height,num_xtile,num_ytile,last_x,last_y,xpxsz,ypxsz,b = getTileInfo(temp_inimage,temp1,x,y,pb)
                    for k in range(1,num_xtile+1):
                        for j in range(1,num_ytile+1):
                            outname = getTileName(temp_inimage,num_xtile,num_ytile,k,j)
                            print '\nInput Image: {}\nImage Size: {}, {}\nOutput Image: {}.tif\nCols/Rows: {}x{}'.format(os.path.basename(temp_inimage),width,height,outname,num_xtile,num_ytile)
                            xstart,ystart,xsize,ysize = getTileCoords(xpxsz,ypxsz,pb,k,j,num_xtile,num_ytile,last_x,last_y)
                            if checkMean(temp_inimage,xstart,ystart,xsize,ysize,empty,num_xtile,num_ytile,datatype):
                                optimizeRaster(temp_inimage,outdir,outname,temp1,temp2,num_xtile,num_ytile,xstart,ystart,xsize,ysize,compression,opts,b,datatype)
                                if k == 1 and j == 1:
                                    if platform == 'win':
                                        try:
                                            image_size_gb = s.getUncompressedSize(temp_outimage)
                                        except Exception,e:
                                            print '\n',e
                                            image_size_gb = 1
                                if gcs != '':
                                    csv_table.append('/vsicurl/http://storage.googleapis.com/{}/{},{}'.format(gcs,os.path.basename(temp_outimage),os.path.basename(temp_outimage).replace('.tif','')))
                            else:
                                image_size_gb = 1
                                print 'Output tile contains only NoData and will be skipped.'
                            j+=1
                            continue
                        k+=1
                        continue

                    if gcs != '':
                        uploadGCS(outdir,gcs)

                    removeFile(os.path.join(temp1,'*.*'))

                    if platform == 'win':
                        pointer_h = open(pointer_txt,'w')
                        pointer_workdir = c.getUNCPath(workdir)
                        pointer_h.write('%s,%s,%s,%s'%(pointer_workdir,1,image_size_gb,fullname))
                        pointer_h.close()
                        ptime_tile = time.time() - stime_tile
                        c.createDoneFile(wrk_file,id=s_id,ptime_tile=ptime_tile,wrk_file_h = '',tile_size_gb=image_size_gb,jobno=jobno)

                if gcs != '':
                    print 'Writing table...'
                    logfile = os.path.join(temp2,'_table.txt')
                    with open(logfile, 'a+') as lf:
                        for item in csv_table:
                            lf.write('{}\n'.format(item))

                    if platform == 'win':
                        i.copyFile(logfile,outdir,exclude='/xf *.tif* *.vrt')
                    elif platform == 'linux':
                        cp_cmd = 'cp {} {}'.format(logfile,outdir)
                        print 'Copying table...'
                        os.system(cp_cmd)

                removeFile(os.path.join(temp2,'*.*'))

if platform == 'win':
    c.printProcessingTime(stime)
    c.removeLockFile(activef,activef_h)
