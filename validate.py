def validate(filename, datatype):
	if int(gdal.VersionInfo('VERSION_NUM')) < 2020000:
        raise ValidateCloudOptimizedGeoTIFFException(
            "GDAL 2.2 or above required")

	gdal.PushErrorHandler()
	ds = gdal.Open(filename)
	gdal.PopErrorHandler()
		
	main_band = ds.GetRasterBand(1)
	block_size = main_band.GetBlockSize()
	ovr_count = main_band.GetOverviewCount()
	ifd_offset = [int(main_band.GetMetadataItem('IFD_OFFSET', 'TIFF'))]
	
	if ds is None:
		raise ValidateCloudOptimizedGeoTIFFException(
			"Invalid file : %s" % gdal.GetLastErrorMsg())
			
	if ds.GetDriver().ShortName != 'GTiff':
		raise ValidateCloudOptimizedGeoTIFFException(
			"The file is not a GeoTIFF")
	
	if block_size[0] == main_band.XSize:
		raise ValidateCloudOptimizedGeoTIFFException(
			"Full resolution image is not tiled")
	
	if ifd_offset[0] != 8:
		raise ValidateCloudOptimizedGeoTIFFException(
			"The offset of the main IFD should be 8. It is %d instead" %
			ifd_offset[0])
	
	if datatype != '3':
		if block_size[0] != '512':
			raise ValidateCloudOptimizedGeoTIFFException(
				"Full resolution image is not tiled with 512 block size")
				
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

			ovr_block_size = ovr_band.GetBlockSize()
			if ovr_block_size[0] == ovr_band.XSize:
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
	else:
		if block_size[0] != '256':
			raise ValidateCloudOptimizedGeoTIFFException(
				"The tif-cache is not tiled with 256 block size")
				
	ds = None
