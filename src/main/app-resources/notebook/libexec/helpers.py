from urlparse import urlparse
import cioppy
import geopandas as gp
import pandas as pd
import os 
import gdal
import osr
from shapely.wkt import loads
from shapely.geometry import box

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