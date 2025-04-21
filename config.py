from pydantic import BaseModel
from enum import Enum

"""
Configuration file for running "gpx_route_matching.py" and "fig_gen.py"
Attributes marked with "Exists" must exist before gpx_route_matching.py
is called, those marked with "Created" will be created/cached and
need not exist before the program is called
"""

# Not actually used by the programs but it's nice to know where
# they are in any case and it could be used with e.g. 
# subprocess.Popen() if desired 
class OSRMServerConfig(BaseModel):
    server_startup_filename: str  # Exists
    pbf_filename: str # Exists
    profile: str

class SwissTLM3D(BaseModel):
    swiss_tlm_gpkg: str # Exists
    clip_gpkg: str # Created
    network_layer: str # Created

class GPXTracks(BaseModel):
    gpx_track_file: str # Exists
    gdf_gpx_path: str # Created

class Paths(BaseModel):
    taken_paths_gpkg: str # Created
    alternate_routes_prefix: str # created

class Name(BaseModel):
    name: str
    photo_loc: str

class Config(BaseModel):
    osrm_server_config: OSRMServerConfig
    swisstlm3d: SwissTLM3D
    gpx_tracks: GPXTracks
    paths: Paths
    name: Name
