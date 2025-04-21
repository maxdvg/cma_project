"""
Examining route preferences for recreational hiking and cycling
The following code was written without any AI 
"""

import gpxpy
import sys
import subprocess
import geopandas as gpd
import os.path
import numpy as np
import shapely
import matplotlib.pyplot as plt
import typing
import requests
import yaml

from config import Config
from tqdm import tqdm
from pyproj import Transformer
from utils import gpx_point_to_shapely, osrm_response_point_to_shapely, gpx_route_to_gdf


def osrm_alternate_routes_query(start_point, end_point, profile: str, alternatives=3):
    osrm_base = f"http://localhost:5000/route/v1/{profile}/"
    def to_str(point):
        point = point["geometry"]
        return f"{point.x},{point.y}"
    coords = f"{to_str(start_point)};{to_str(end_point)}"
    alternates_req = f"?alternatives={alternatives}"
    geojson_req = "geometries=geojson"
    detailed_overview = "overview=full"
    return osrm_base + coords + alternates_req + "&" + geojson_req + "&" + detailed_overview


def get_osrm_routes(start_point, end_point, profile: str):
    """
    :requires: OSRM server must be running at localhost:5000
    :returns: a list of GeoDataFrames containing the routes as a series of shapely 
    points in EPSG:2056
    """
    query = osrm_alternate_routes_query(start_point, end_point, profile)
    osrm_response = requests.get(query).json()
    routes = []
    if osrm_response['code'] == 'Ok':
        for unparsed_route in osrm_response['routes']:
            route = []
            wgs_path_coords = unparsed_route['geometry']['coordinates']
            for cord_pair in wgs_path_coords:
                route.append(osrm_response_point_to_shapely(cord_pair))
            routes.append(gpd.GeoDataFrame(geometry=route, crs='epsg:2056'))
    return routes

def construct_bounding_box(gpx: gpxpy.gpx.GPX, output_CRS = "EPSG:2056"):
    """
    Constructs a (maximally tight) bounding box around a GPX track
    :param gpx: parsed GPX file which is to be bounded
    :returns: the bounding box represented as a tuple of the format 
    (minx, miny, maxx, maxy)
    """
    bbox_east: float = sys.float_info.min
    bbox_west: float = sys.float_info.max
    bbox_north: float = sys.float_info.min
    bbox_south: float = sys.float_info.max
    
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                if bbox_east < point.longitude:
                    bbox_east = point.longitude
                if bbox_west > point.longitude:
                    bbox_west = point.longitude
                if bbox_north < point.latitude:
                    bbox_north = point.latitude
                if bbox_south > point.latitude:
                    bbox_south = point.latitude

    transformer = Transformer.from_crs("EPSG:4326", output_CRS, always_xy=True)
    bbox_west, bbox_south = transformer.transform(bbox_west, bbox_south)
    bbox_east, bbox_north = transformer.transform(bbox_east, bbox_north)
    return (bbox_west, bbox_south, bbox_east, bbox_north)


def clip_gpkg(bbox, 
              gpkg,
              output_file,
              buffer = 10000,
              ):
    """
    Clips a GeoPackage file `gpkg` to the confines of a bounding box `bbox` with a buffer
    of `buffer` meters. 
    :requires: ogr2ogr must be installed and executable from the working directory of this project
    :param bbox: bounding box, as a tuple, (minx, miny, maxx, maxy)
    :param gpkg: Path to GeoPackage file to be clipped
    :param output_file: path to file to which the created clipped gpkg should be saved
    :param buffer: buffer in meters around the bounding box which should also be included in the generated clipped geopackage
    :side-effects: if `output_file` does not exist, creates said output file and populates it with the clipped GPKG, 
    otherwise if `output_file` exists, prompt user whether they wish to overwrite with clipped contents or not
    :returns: `output_file` filename, see "side-effects" above
    """
    minx = bbox[0] - buffer
    miny = bbox[1] - buffer
    maxx = bbox[2] + buffer
    maxy = bbox[3] + buffer
    cmd = f"ogr2ogr -clipdst {minx} {miny} {maxx} {maxy} {output_file} {gpkg}"
    print(f"cmd: {cmd}")
    overwrite = "y"
    if os.path.isfile(output_file):
        overwrite = input(f"{output_file} already exists, overwrite it (y/n)?: ")
        if overwrite == "y":
            print(f"Running '{cmd}'")
            subprocess.call(cmd, shell=True)
    else:
        print(f"Running '{cmd}'")
        subprocess.call(cmd, shell=True)
    return output_file 

def find_closest_path_to_point(point: shapely.geometry.point.Point,
                               paths_geopackage: gpd.geodataframe.GeoDataFrame):
    """
    Given a point and a collection of paths, find the path in the collection 
    which is closest to said point
    :param point: A shapely point for which the nearest path in `paths_geopackage`
    should be found
    :param paths_geopackage: geopackage containing paths (LineStrings) to which distances
    to `point` should be computed
    :returns: A tuple (path_idx, path) where
            - path_idx is the index of the nearest path in `paths_geopackage` to `point` 
            in `paths_geopackage`
            - path is the entire entry of the path in the `paths_geopackage`
    """
    distances: np.ndarray = paths_geopackage.distance(point).values
    path_idx = np.argmin(distances)
    return path_idx, paths_geopackage.iloc[path_idx]

def find_paths_taken(gpx_route_gdf: gpd.geodataframe.GeoDataFrame,
                     all_paths: gpd.geodataframe.GeoDataFrame,
                     cache: typing.Optional[str] = None,
                     path_cutoff: int = 10):
    """
    Given a geodataframe of shapely points `gpx_route_gdf` representing some track
    and a geodataframe of paths `all_paths` near
    said GPX track, find the paths in the geodataframe which were most likely taken
    during the recording of the GPX track

    Naive implementation, simply assumes that the user took the tracks closest
    to the GPX points. A more advanced implementation might use hidden markov models.

    :param gpx_route_gdf: GPX route as a GeoDataFrame of Shapely points, as returned by
    gpx_route_to_gdf from utils.py
    :param all_paths: The database of all paths to which the gpx route could be matched
    :param cache: The name of the file from which / to which the matched paths should
    be saved. If the file exists, loads from said file. If the file does not exist,
    stores result of calculating (paths taken) as GPKG
    :param path_cutoff: number of GPX points which must match with a path for that 
    path to be considered part of the route
    :returns: geodataframe of all of the paths in `all_paths` which match with at least
    `path_cutoff` paths in `gpx_route_gdf`
    """
    if cache is None or not os.path.isfile(cache):
        taken_path_idxs = {}
        for _ , p in tqdm(gpx_route_gdf.iterrows(),
                          desc="Finding closest paths to GPX route",
                          unit="GPX points"):
            point: shapely.geometry.point.Point = p["geometry"]
            closest_path_idx, _ = find_closest_path_to_point(point,
                                                             all_paths)
            if closest_path_idx not in taken_path_idxs:
                taken_path_idxs[closest_path_idx] = 1
            else:
                taken_path_idxs[closest_path_idx] += 1
        path_idxs_that_meet_threshold = []
        for p_idx, num_times_taken in taken_path_idxs.items():
            if num_times_taken >= path_cutoff:
                path_idxs_that_meet_threshold.append(p_idx)
        taken_paths = [all_paths.iloc[p_idx] for p_idx in path_idxs_that_meet_threshold]
        taken_paths = gpd.GeoDataFrame(taken_paths)
        taken_paths = taken_paths.set_crs('epsg:2056')
        if cache is not None:
            taken_paths.to_file(cache, driver='GPKG', layer='Paths') 
    else:
        taken_paths = gpd.read_file(cache, layer='Paths')    

    return taken_paths

def path_statistics(path_gpkg: str):
    """
    For a given `path_gpkg` which is a geopackage containing all of the path
    segments for a given path with all the data from SwissTLM3D, calculate
    some statistics regarding the path segments
    :param path_gpkg: Filename of a geopackage containing the SwissTLM3D path
    segments which are to be analyzed
    :returns: A nested dictionary mapping from interesting attribute types to
    a dictionary mapping from that attribute's subtypes to path lengths of said
    subtype
    """
    taken_paths = gpd.read_file(path_gpkg)
    # Attributes of interest from SwissTLM3D
    interesting_attributes = ["objektart",
                              "wanderwege",
                              "befahrbarkeit",
                              "richtungsgetrennt",
                              "belagsart",
                              "eigentuemer"]
    stats_dict = {}
    for attribute in interesting_attributes:
        stats_dict[attribute] = {}
    for i in range(len(taken_paths)):
        sub_path = taken_paths.iloc[i]
        sub_path_len = sub_path.geometry.length
        for attribute in interesting_attributes:
            attr_val = sub_path[attribute]
            attr_stats_dict = stats_dict[attribute]
            if attr_val not in attr_stats_dict:
                attr_stats_dict[attr_val] = sub_path_len
            else:
                attr_stats_dict[attr_val] += sub_path_len
    # Replace None type key with "No Value Given" 
    # hopefully avoid downstream runtime errors
    for attribute, attr_stats_dict in stats_dict.items():
        if None in attr_stats_dict and "N/A" not in attr_stats_dict:
            attr_stats_dict["N/A"] = attr_stats_dict[None]
            del attr_stats_dict[None]
        elif None in attr_stats_dict and "N/A" in attr_stats_dict:
            raise RuntimeError("Whoever made the dataset you're working with is an idiot")

    return stats_dict

if __name__ == "__main__":
    """
    sys.argv[1] must be a config yaml file as layed specified by
    config.py. See config.yaml for an example of such a config

    General Warning: All of the code working with Geopandas in this
    project assumes that the indices remain untouched. Should the
    indices of any Geopandas DataFrame be modified during or between
    the execution of any of the function provided, no guarantees can
    be made regarding behaviour. Most likely, things would break.
    """

    # Load in the config file
    with open(sys.argv[1], "r") as f:
        config_data = yaml.safe_load(f)
    config = Config.model_validate(config_data)

    # Existing files
    swisstlm = config.swisstlm3d.swiss_tlm_gpkg
    personal_gps_track = config.gpx_tracks.gpx_track_file
    # Files to be created
    clipped_swisstlm = config.swisstlm3d.clip_gpkg
    taken_paths_gpkg = config.paths.taken_paths_gpkg


    """
    GPX Track route matching and display
    """
    # Load in GPX track and extract the data
    # assumes data is in format as created when exporting
    # from Garmin Connect to GPX
    with open(personal_gps_track, 'r') as fluntern_gpx_file:
        gpx = gpxpy.parse(fluntern_gpx_file)
        route = gpx.tracks[0].segments[0]

    # Create a bounding box to construct a clipped gpkg of SwissTLM3D
    # otherwise a lot of operations take an extremely long time due to
    # how giant SWISSTLM3D is
    bbox_gpx = construct_bounding_box(gpx)
    print(f"bounding box of hike: {bbox_gpx}")
    clipped_gpkg = clip_gpkg(bbox_gpx,
                             swisstlm,
                             clipped_swisstlm,
                             buffer=8000)

    # Extract streets and paths from the newly created SwissTLM3D
    streets_layer = config.swisstlm3d.network_layer
    geopkg = gpd.read_file(clipped_gpkg, layer=streets_layer)

    route_gdf = gpx_route_to_gdf(route, cache=config.gpx_tracks.gdf_gpx_path)
    # Find all the paths taken on the GPX route
    paths = find_paths_taken(route_gdf, geopkg, cache=taken_paths_gpkg)

    """
    Do some basic plotting, make sure everything functioning reasonably
    """

    # plot paths taken
    fig, ax = plt.subplots()
    geopkg.plot(ax=ax, color="yellow")
    paths.plot(ax=ax, color="blue")
    plt.savefig(f"{config.name.photo_loc}_path.png")

    """
    Finding and plotting alternate routes
    """
    wgs_route = route_gdf.to_crs("epsg:4326")
    # Example: find alternate routes from start to end point
    osrm_routes = get_osrm_routes(wgs_route.iloc[0], wgs_route.iloc[len(wgs_route) - 1], config.osrm_server_config.profile)
    # Get the alternate routes and save them for sanity check
    for or_idx, osrm_route in enumerate(osrm_routes):
        alternate_route_filename = f"{config.paths.alternate_routes_prefix}_{or_idx}.gpkg"
        fig, ax = plt.subplots()
        geopkg.plot(ax=ax, color="black")
        osrm_route.plot(ax=ax, color="yellow")
        paths = find_paths_taken(osrm_route, geopkg, cache=alternate_route_filename, path_cutoff=3)
        paths.plot(ax=ax, color="blue")
        plt.savefig(f"{config.name.photo_loc}_{config.paths.alternate_routes_prefix}_{or_idx}.png")



