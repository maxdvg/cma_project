import gpxpy
import shapely
from pyproj import Transformer
import geopandas as gpd
from tqdm import tqdm
from os.path import isfile

"""
A couple of utility functions
"""

def gpx_point_to_shapely(gpx_point: gpxpy.gpx.GPXTrackPoint):
      """
      Converts GPX point `gpx_point` to a shapely point, stripping additional data like
      time and elevation, and changing to the swiss coordinate system from WGS.
      """
      transformer = Transformer.from_crs("EPSG:4326", "EPSG:2056", always_xy=True)
      swiss_point = transformer.transform(gpx_point.longitude, gpx_point.latitude)
      return shapely.Point(swiss_point)

def osrm_response_point_to_shapely(resp_point: list[float]):
      """
      Converts an OSRM response point, like [8.3221, 46.242] into a shapely point in
      ESPG:2056 format
      """
      transformer = Transformer.from_crs("EPSG:4326", "EPSG:2056", always_xy=True)
      swiss_point = transformer.transform(resp_point[0], resp_point[1])
      return shapely.Point(swiss_point)

def gpx_route_to_gdf(gpx_route: gpxpy.gpx.GPXTrackSegment, cache: str = None):
      """
      Converts a GPX route into a geodataframe
      """
      if cache is not None and isfile(cache):
            gpx_gdf = gpd.read_file(cache)
      else:
            shapely_points = []
            for point in tqdm(gpx_route.points,
                              desc="Converting GPX to GDF",
                              unit="Points"):
                  shapely_points.append(gpx_point_to_shapely(point))
            gpx_gdf = gpd.GeoDataFrame(geometry=shapely_points, crs='epsg:2056')
            if cache is not None:
                  gpx_gdf.to_file(cache, driver='GPKG', layer='GPX_ROUTE')
      return gpx_gdf
