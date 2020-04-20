from urlparse import urlparse
import cioppy
import geopandas as gp
import pandas as pd
import os 
import gdal
import osr
from shapely.wkt import loads
from shapely.geometry import box

gdal.UseExceptions()
import sys

sys.path.append('/opt/OTB/lib/python')
sys.path.append('/opt/OTB/lib/libfftw3.so.3')
os.environ['OTB_APPLICATION_PATH'] = '/opt/OTB/lib/otb/applications'
os.environ['LD_LIBRARY_PATH'] = '/opt/OTB/lib'
os.environ['ITK_AUTOLOAD_PATH'] = '/opt/OTB/lib/otb/applications'
os.environ['GDAL_DATA'] = '/opt/anaconda/share/gdal/'
import otbApplication


def get_vsi_url(input_reference, username, api_key):
    
    enclosure = cioppy.Cioppy().search(input_reference, 
                                       [],
                                       'enclosure',
                                       'GeoTime', 
                                       creds='{}:{}'.format(username, api_key))[0]['enclosure']   
    
    parsed_url = urlparse(enclosure)

    vsi_url = '/vsicurl/{}://{}:{}@{}{}'.format(urlparse(enclosure).scheme,
                        username, 
                        api_key,
                        urlparse(enclosure).netloc,
                       urlparse(enclosure).path)
    
    return vsi_url 


def get_product_metadata(input_references, username, api_key):
    
    temp_searches = []

    for index, reference in enumerate(input_references):

        search_temp = gp.GeoDataFrame(cioppy.Cioppy().search(end_point=reference,
                                              params=[],
                                               output_fields='self,track,enclosure,identifier,wkt,startdate,enddate,platform,cc', 
                                               model='EOP',
                                                            creds='{}:{}'.format(username, api_key)))


        temp_searches.append(search_temp)

    search = gp.GeoDataFrame(pd.concat(temp_searches, ignore_index=True)) 

    search['geometry'] = search['wkt'].apply(loads)
    search['cc'] = pd.to_numeric(search['cc'])
    search['startdate'] = pd.to_datetime(search['startdate'])
    search['enddate'] = pd.to_datetime(search['enddate'])
    
    return search


def get_image_wkt(product):
    
    src = gdal.Open(product)
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()

    max_x = ulx + (src.RasterXSize * xres)
    min_y = uly + (src.RasterYSize * yres)
    min_x = ulx 
    max_y = uly

    source = osr.SpatialReference()
    source.ImportFromWkt(src.GetProjection())

    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(source, target)

    result_wkt = box(transform.TransformPoint(min_x, min_y)[0],
                     transform.TransformPoint(min_x, min_y)[1],
                     transform.TransformPoint(max_x, max_y)[0],
                     transform.TransformPoint(max_x, max_y)[1]).wkt
    
    return result_wkt


def cog(input_tif, output_tif):
    
    translate_options = gdal.TranslateOptions(gdal.ParseCommandLine('-co TILED=YES ' \
                                                                    '-co COPY_SRC_OVERVIEWS=YES ' \
                                                                    ' -co COMPRESS=LZW'))

    ds = gdal.Open(input_tif, gdal.OF_READONLY)

    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
    ds.BuildOverviews('NEAREST', [2,4,8,16,32])
    
    ds = None

    ds = gdal.Open(input_tif)
    gdal.Translate(output_tif,
                   ds, 
                   options=translate_options)
    ds = None

    os.remove('{}.ovr'.format(input_tif))
    os.remove(input_tif)
    
def pre_proces(input_references, bbox, ndvi_threshold, bsi_threshold, output_name, username, api_key):
    
    
    x_min, y_min, x_max, y_max = bbox
    
    # Create a VRT with all products provided 
    vsi_list = [get_vsi_url(input_reference, 
                            username,
                            api_key) for input_reference in input_references]

    vrt_options = gdal.BuildVRTOptions(VRTNodata=-20)

    vrt = 'my.vrt'

    if os.path.exists(vrt):
        os.remove(vrt)
    
    gdal.UseExceptions()
    ds_vrt = gdal.BuildVRT(vrt, vsi_list, options=vrt_options)
    
    while ds_vrt is None:
        ds_vrt = gdal.Open(vrt,  gdal.OF_UPDATE)
    
    ds_vrt.FlushCache()
    
    ds_vrt = None
    
    
    # Update the VRT bands description
    src_ds = gdal.Open(get_vsi_url(input_references[0], 
                                   username,
                                   api_key))

    ds_vrt = gdal.Open(vrt,  gdal.OF_UPDATE)

    for band in range(ds_vrt.RasterCount):

        band += 1
        print band, src_ds.GetRasterBand(band).GetDescription()

        src_band = ds_vrt.GetRasterBand(band)
        src_band.SetDescription(src_ds.GetRasterBand(band).GetDescription())  

    ds_vrt.FlushCache()
    
    ds_vrt = None
    src_ds = None
    
    # clip the VRT with the AOI provided
    clipped_vrt = 'clipped.vrt'
    
    if os.path.exists(clipped_vrt):
        os.remove(clipped_vrt)

    gdal.Translate(clipped_vrt,
                   vrt,
                   projWin=[x_min, y_max, x_max, y_min],
                   projWinSRS='EPSG:4326',
                   format='VRT')
    
    # Build the expressions
    band_2 = 'im1b2'
    band_4 = 'im1b4'
    band_8 = 'im1b8'
    band_11 = 'im1b11'
    scl = 'im1b13'
    
    
    ndvi_expression = '({1}-{0})/({1}+{0})'.format(band_4, band_8)
    bsi_expression = '(({3}/{4}+{1}/{4})-({2}/{4}+{0}/{4}))/(({3}/{4}+{1}/{4})+({2}/{4}+{0}/{4}))'.format(band_2, band_4, band_8, band_11, '10000')
    invalid_expression = '{0} == 0 || {0} == 1 || {0} == 3 || {0} == 8 || {0} == 9'.format(scl) 


    # if ndvi  >= 0.3 and bsi <= 0 then IS vegetation, otherwise no vegetation
    mask_expression = '{0} ? 128 : {1} >= {2} && {3} <= {4} ? 10 : 20'.format(invalid_expression,
                                                                              ndvi_expression, 
                                                                              ndvi_threshold,
                                                                              bsi_expression,
                                                                              bsi_threshold)
    # NDVI MASK & BSI MASK added-----
    bsi_mask_expression='{0} ? 128 :{1} <= {2}? 1:0'.format(invalid_expression,bsi_expression,bsi_threshold)
    ndvi_mask_expression='{0} ? 128 :{1} >= {2}? 1:0'.format(invalid_expression,ndvi_expression,ndvi_threshold)

    
    
    expressions = []

    invalid_data = 128
    expressions.append(mask_expression)
    
    expressions.append('{} ? {} : {}'.format(invalid_expression, invalid_data, bsi_expression))
    expressions.append(bsi_mask_expression)
    
    expressions.append('{} ? {} : {}'.format(invalid_expression, invalid_data, ndvi_expression))
    expressions.append(ndvi_mask_expression)
    
    expressions.append('{} ? {} : 0'.format(invalid_expression, invalid_data))

    
    
    
    
    
    
    band_names = ['vegetation_mask',
                  'bsi_values',
                  'bsi_mask',
                  'ndvi_values',
                  'ndvi_mask',
                  'cloud_mask']

    metadata = dict()
    metadata['B02'] = 'im1b2'
    metadata['B04'] = 'im1b4'
    metadata['B08'] = 'im1b8'
    metadata['B11'] = 'im1b11'
    metadata['SCL'] = 'im1b13'

    
    for index, expression in enumerate(expressions):
        # Apply the BandMathX OTB operator
        BandMathX = otbApplication.Registry.CreateApplication('BandMathX')

        BandMathX.SetParameterStringList('il', [clipped_vrt])
        BandMathX.SetParameterString('out', '{}_{}'.format(index, output_name))
        BandMathX.SetParameterString('exp', expression)
        BandMathX.SetParameterOutputImagePixelType('out', otbApplication.ImagePixelType_double)

        BandMathX.ExecuteAndWriteOutput()
        
        ds_temp = gdal.Open('{}_{}.tif'.format(index, output_name), gdal.OF_UPDATE)

        for band_index in range(ds_temp.RasterCount):


            metadata['BAND_EXPRESSION'] = '({})'.format(expression)

            src_band = ds_temp.GetRasterBand(band_index+1)
            src_band.SetMetadata(metadata)
            src_band.SetDescription(band_names[index])  
            src_band.SetRasterColorInterpretation(gdal.GCI_Undefined)

        ds_temp.FlushCache()
    
        ds_temp = None
    
    os.remove(vrt)
    os.remove(clipped_vrt)
    
def get_mosaic_expressionx(no_data, image_count):
    nodata_expr = '(' + ' && '.join(['im{}b1 == 128'.format(n) for n in range(1, image_count + 1)]) + ') ? {} : '.format(no_data)
    
    c = []  
    for n in range(1, image_count + 1): 
        c.append('(im{0}b1 >= -1 && im{0}b1 <= 1) ? im{0}b1'.format(n))
    
    return nodata_expr + ' : '.join(c) + ' : 255'



def get_mask_expressionx(no_data, image_count, classes):
    
    nodata_expr = '(' + ' && '.join(['im{}b1 == 128'.format(n) for n in range(1, image_count + 1)]) + ') ? {} : '.format(no_data)
    
    c = []  
    for k in classes: 
        
        c.append('(' + ' || '.join(['im{0}b1 == {1}'.format(n, k) for n in range(1, image_count + 1)]) + ') ? {}'.format(k))
    
    return nodata_expr + ' : '.join(c) + ' : 255'